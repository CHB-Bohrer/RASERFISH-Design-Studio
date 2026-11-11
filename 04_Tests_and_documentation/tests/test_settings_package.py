"""Portable settings roundtrips and hostile/incomplete input bundles."""
import copy
import gzip
import io
import json
import shutil
import subprocess
import unittest
import warnings
import zipfile
from unittest.mock import patch

from compact_outputs import SCHEMA
from design_engine import Policy
from sequence_tools import make_bins, validate_imagers
from studio import BASE, bootstrap, export_settings, import_settings, normalize_config_inputs


def archive_bytes(entries):
    output=io.BytesIO()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',UserWarning)
        with zipfile.ZipFile(output,'w',zipfile.ZIP_DEFLATED) as archive:
            for name,data in entries:
                archive.writestr(name,data)
    return output.getvalue()


def rewrite_archive(data,change):
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        entries={name:archive.read(name) for name in archive.namelist()}
    change(entries)
    return archive_bytes(entries.items())


class SettingsPackageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        boot=bootstrap()
        cls.defaults=copy.deepcopy(boot['default_config'])
        cls.defaults['references']=copy.deepcopy(boot['references'])
        cls.package=export_settings(cls.defaults)

    def expected_config(self,config):
        expected=normalize_config_inputs(config)
        expected['policy']=vars(Policy.from_dict(expected.get('policy',{})))
        expected['imagers']=validate_imagers(expected['imagers'])
        return expected

    def test_verified_default_configuration_roundtrips_exactly(self):
        restored=import_settings(self.package)
        self.assertEqual(restored,self.expected_config(self.defaults))
        self.assertEqual(restored['policy']['bound_strand'],'+')
        self.assertEqual(restored['policy']['tm_min'],45)
        self.assertEqual(restored['policy']['tm_max'],55)
        self.assertEqual(restored['bin_width'],20000)
        self.assertEqual(list(restored['bounds'].values()),[[189019,209018]])
        bins,_=make_bins(restored['references'],restored['bin_width'],restored['tail_policy'],
                         restored['assignments'],restored['bounds'])
        self.assertEqual([b['imager'] for b in bins],['DP1'])

    def test_save_does_not_mutate_the_live_configuration(self):
        config=copy.deepcopy(self.defaults);before=copy.deepcopy(config)
        export_settings(config)
        self.assertEqual(config,before)

    def test_manifest_is_compact_and_full_context_is_compressed_once(self):
        with zipfile.ZipFile(io.BytesIO(self.package)) as archive:
            names=archive.namelist();manifest=json.loads(archive.read('run_config.json'))
            self.assertEqual(manifest['schema'],SCHEMA)
            self.assertIn('inputs/target_region.fa.gz',names)
            self.assertIn('inputs/reference_context.fa.gz',names)
            self.assertNotIn(self.defaults['references'][0]['sequence'],json.dumps(manifest))
            self.assertFalse(any(name.endswith(('.nhr','.nin','.nsq','.ndb')) for name in names))
            self.assertLess(len(archive.read('run_config.json')),100000)

    def test_background_and_rna_fasta_inputs_survive_roundtrip(self):
        config=copy.deepcopy(self.defaults)
        config['background_text']='>construct1 sequence used in vector\nACGTGCTAGCTAGCTAGCAT\n>construct2\nTTAGCGATCGATCGTA'
        config['transcript_text']='>rna1 annotated test RNA\nacguugcauugccaug\n>rna2\nUUUACGUACGU'
        expected=self.expected_config(config)
        restored=import_settings(export_settings(config))
        self.assertEqual(restored,expected)
        self.assertEqual(restored['transcripts'][0]['sequence'],'ACGTTGCATTGCCATG')
        self.assertEqual(restored['transcripts'][0]['description'],'rna1 annotated test RNA')
        self.assertNotIn('transcript_text',restored)
        self.assertNotIn('background_text',restored)

    def test_existing_additional_record_metadata_and_database_paths_survive(self):
        config=copy.deepcopy(self.defaults)
        config['background']=[dict(id='vector',sequence='ACGT'*20,source='laboratory vector',
                                   features=[dict(start=1,end=8,type='promoter')])]
        config['transcripts']=[dict(id='RNA_X',sequence='TGCA'*20,description='Supplied transcript',
                                   provenance={'accession':'synthetic_fixture'})]
        config['databases']={'human':'/a path/genome','rna_human':'/another path/rna'}
        restored=import_settings(export_settings(config))
        self.assertEqual(restored,self.expected_config(config))

    def test_missing_manifest_is_rejected(self):
        with self.assertRaises((OSError,ValueError)):
            import_settings(archive_bytes([('inputs/placeholder','nothing')]))

    def test_manifest_without_required_sequence_asset_is_rejected(self):
        def change(entries):entries.pop('inputs/reference_context.fa.gz')
        with self.assertRaises((OSError,ValueError)):
            import_settings(rewrite_archive(self.package,change))

    def test_manifest_without_required_annotation_asset_is_rejected(self):
        def change(entries):entries.pop('inputs/annotations.json.gz')
        with self.assertRaises((OSError,ValueError)):
            import_settings(rewrite_archive(self.package,change))

    def test_changed_input_asset_fails_hash_verification(self):
        def change(entries):entries['inputs/reference_context.fa.gz']=gzip.compress(b'>context_1\nACGT\n')
        with self.assertRaisesRegex(ValueError,'SHA-256'):
            import_settings(rewrite_archive(self.package,change))

    def test_manifest_undeclared_sequence_record_is_rejected(self):
        def change(entries):
            manifest=json.loads(entries['run_config.json'])
            manifest['configuration']['references'][0]['sequence']['$sequence']='missing'
            entries['run_config.json']=json.dumps(manifest)
        with self.assertRaisesRegex(ValueError,'Missing manifest sequence'):
            import_settings(rewrite_archive(self.package,change))

    def test_archive_traversal_absolute_and_windows_paths_are_rejected(self):
        for filename in ('../escape','inputs/../../escape','/absolute','inputs\\..\\escape'):
            with self.subTest(filename=filename),self.assertRaisesRegex(ValueError,'Unsafe'):
                import_settings(archive_bytes([(filename,b'bad')]))

    def test_duplicate_archive_member_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'duplicate'):
            import_settings(archive_bytes([('run_config.json',b'{}'),('run_config.json',b'{}')]))

    def test_manifest_asset_cannot_escape_extracted_directory(self):
        def change(entries):
            manifest=json.loads(entries['run_config.json'])
            manifest['assets']['../outside']={'sha256':'0'*64,'format':'fasta.gz'}
            entries['run_config.json']=json.dumps(manifest)
        with self.assertRaisesRegex(ValueError,'safe relative path'):
            import_settings(rewrite_archive(self.package,change))

    def test_oversize_upload_is_rejected_before_zip_read(self):
        with self.assertRaisesRegex(ValueError,'25 MB'):
            import_settings(b'x'*25_000_001)

    def test_too_many_archive_files_are_rejected(self):
        with self.assertRaisesRegex(ValueError,'too many files'):
            import_settings(archive_bytes([(f'input-{i}',b'x') for i in range(101)]))

    def test_outer_zip_expansion_is_bounded(self):
        output=io.BytesIO()
        with zipfile.ZipFile(output,'w',zipfile.ZIP_DEFLATED) as archive:
            with archive.open('large.txt','w') as member:
                for _ in range(50):member.write(b'x'*1_000_000)
                member.write(b'x')
        self.assertLess(len(output.getvalue()),100000)
        with self.assertRaisesRegex(ValueError,'expanded-size'):
            import_settings(output.getvalue())

    def test_inner_gzip_expansion_is_bounded(self):
        # Lower the asset limit so the nested-compression boundary is exercised
        # without allocating a production-size decompression bomb in this test.
        with patch('compact_outputs.MAX_ASSET_BYTES',1024):
            with self.assertRaisesRegex(ValueError,'decompressed size'):
                import_settings(self.package)

    def test_legacy_plain_configuration_in_zip_remains_importable(self):
        restored=import_settings(archive_bytes([('run_config.json',json.dumps(self.defaults))]))
        self.assertEqual(restored,self.defaults)

    @unittest.skipUnless(shutil.which('node'),'Node.js is optional; run with Node available for UI settings regression')
    def test_actual_ui_restore_and_save_preserve_additional_screening_inputs(self):
        # Execute the actual two UI functions in a small control harness. The
        # browser layout is tested separately; this exercises preservation logic
        # for records, legacy textarea imports, edits, and subsequent resets.
        script=r'''
const fs=require('fs'),vm=require('vm'),assert=require('assert');
const source=fs.readFileSync(process.argv[1],'utf8');
const namesCode=source.slice(source.indexOf('// Readable run names'),source.indexOf('// End readable run names.'));
const configCode=source.slice(source.indexOf('function configuration()'),source.indexOf('\nfunction renderRunSummary'));
const restoreCode=source.slice(source.indexOf('function restoreConfiguration('),source.indexOf("\n$('saveConfig').onclick"));
assert(configCode.startsWith('function configuration()'));
assert(restoreCode.startsWith('function restoreConfiguration('));
const controls=new Map();
const get=id=>{if(!controls.has(id))controls.set(id,{value:'',checked:false,type:'text',classList:{toggle(){}}});return controls.get(id);};
const ctx={S:{bins:[]},B:{presets:{errfi1:{description:'ERRFI1'}},environment:{databases:{}}},
  $:get,copy:x=>JSON.parse(JSON.stringify(x)),resolveImagerEntry:x=>x,
  loadedRefs(){},makeBins(){},renderAllowed(){},updateAdaptiveControls(){},markStale(){},
  getBounds:()=>ctx.cfg.bounds,getPolicy:()=>ctx.cfg.policy};
vm.createContext(ctx);vm.runInContext(namesCode+'\n'+configCode+'\n'+restoreCode,ctx);
const base={references:[{id:'X',sequence:'ACGT',length:4}],imagers:[],policy:{preset_id:'errfi1'},bounds:{X:[1,4]}};
const roundtrip=extra=>{ctx.cfg={...base,...extra};vm.runInContext('restoreConfiguration(cfg)',ctx);return JSON.parse(JSON.stringify(vm.runInContext('configuration()',ctx)));};
let result=roundtrip({background:[{id:'vector',description:'vector original description',sequence:'ACGT',source:'original metadata'}],transcripts:[{id:'rna',sequence:'TGCA',source:'original RNA metadata'}]});
assert.strictEqual(result.background[0].source,'original metadata');
assert.strictEqual(result.transcripts[0].source,'original RNA metadata');
assert.strictEqual(get('constructFasta').value,'>vector original description\nACGT');
get('constructFasta').value='>edited\nAAAA';result=JSON.parse(JSON.stringify(vm.runInContext('configuration()',ctx)));
assert.strictEqual(result.background_text,'>edited\nAAAA');assert(!('background' in result));
result=roundtrip({background_text:'>legacy vector\nCCGG',transcript_text:'>legacy RNA\nAUGC'});
assert.strictEqual(result.background_text,'>legacy vector\nCCGG');
assert.strictEqual(result.transcript_text,'>legacy RNA\nAUGC');
assert(!('background' in result));assert(!('transcripts' in result));
result=roundtrip({});
assert(!result.background?.length);assert(!result.transcripts?.length);
assert(!result.background_text);assert(!result.transcript_text);
'''
        completed=subprocess.run([shutil.which('node'),'-e',script,str(BASE/'web/app.js')],capture_output=True,text=True)
        self.assertEqual(completed.returncode,0,completed.stdout+completed.stderr)


if __name__=='__main__':unittest.main()
