Packaging note (2026-11-11): This is a dated development record. See README.md for current setup and PACKAGE_NOTES.md for this source snapshot. Paths to private evidence have been omitted; those files are not included.

# RASER-FISH Design Studio 1.6 - current methods

This is a candidate designer: QC participates before and during selection. It is
not just the older audit viewer and not the full published OligoMiner pipeline.
Each run preserves exact inputs through a hashed manifest and compressed assets.
The v1.6 defaults were checked against saved ERRFI1 run 20261015_100112_62f902
(verification/errfi1_baseline.json). Scientific acceptance thresholds are unchanged;
search budgets/methods and output storage have changed.

## 1. Target sequence, layout and readout

Inputs are GRCh38 genomic coordinates, DNA/FASTA or annotated GenBank/SnapGene files.
The supplied reference is the top (+) strand; genomic mapping orientation is retained.
IUPAC characters retain their coordinates; arms containing ambiguity are excluded.
Independent references are not arbitrarily concatenated. Circular/integration
junctions require explicit sequence. Local targets are limited to 1 Mb per run.

Tracing bins are equal width, adjacent within one bounded span and nonoverlapping.
A remainder is reported/excluded or causes a strict-division stop. The 20 kb ERRFI1
windows are labeling widths, not optical precision. The existing UI still supports
one continuous bounded span per reference per run; five separated windows require
separate runs and later combined-pool QC. Importing a BED only adds a feature track.

UCSC sequence, Ensembl gene/transcript/exon context and compatible ENCODE H3K27ac
tracks remain supported. Cell line is supplied by the user. Optional fallback is
requested cell -> HCT116 -> IMR-90; a proxy is labeled with its actual source.
No network error is treated as absence of data. Peak scores and continuous signals
are different. Annotation/ChIP tracks are not measurements of probe accessibility,
transcript expression, fluorescence or purity in the user's cells.

The readout can be entered as the imaging oligo OR its primary docking site, with
an explicit convention. Both stored strings are 5' to 3'; they are reverse
complements. The full primary is docking + TT + arm + TT + docking. Its handles
are present before folding, off-target and pair checks. The diagram draws the
imager antiparallel for pairing, not as an alternative order sequence. Reusing an
imager does not identify separate bins; same-dye distinct imagers need distinct
readout handling. The inventory is not automatically proven orthogonal.

The bundled setup is the verified 400,000 bp forward GRCh38 reference
chr1:7,820,001-8,220,000. Local inclusive bounds 189019-209018 give one strict
20,000 bp region at genomic chr1:8,009,019-8,029,018. Readout DP1 uses imager
5'-CTAACGAATCGG-3' and primary docking 5'-CCGATTCGTTAG-3'. Its reference hash matches
the saved run. An original HCT116 proxy annotation is not a measurement in the
user's cells.

## 2. One trace-wide DNA-bound strand

A production run requires a valid + or - DNA-bound strand. Top (+) is the visible,
verified ERRFI1 default throughout the UI, backend and preset. Bottom (-) remains
an explicit choice; unselected, mixed and automatic modes remain invalid. There
is no per-gene switch. A top-targeting arm is the reverse complement
of the corresponding forward reference window. A bottom-targeting arm equals it.
Every exact repeat copy and final selected site must obey the same trace-wide rule.
Differently oriented fragments of a mapped chromosome are checked for consistency.
All windows in independent runs belonging to the same trace must use the same choice.
This is necessary sequence consistency, not proof of observed chromatid/homolog identity.

## 3. Candidate generation and retained temperature model

Current defaults: preferred 36-42 nt arms, one-base adaptive alternatives through
90 nt, GC 0-76% (0 means no additional lower-GC gate), at most five consecutive
identical bases, and exclusion of eight dinucleotide repeat units. Candidate start
stride is one base. Arms must satisfy strictly 45 < Tm < 55 C under the active model.

Default model: Biopython DNA_NN3, 390 mM monovalent Na, 2 nM EACH probe species and
0 nM target as a limiting trace-target assumption, entropy salt correction method
5, aqueous Tm rounded to 0.01 C before subtracting 50 x 0.65 C for formamide. This is
not measured in-cell Tm or a recovered Ellenberg command line. The NN3 25/25 nM and
previous NN4 2/0 nM comparisons remain available. The 45-55 C window does not change
hybridization temperature, 42 C aqueous folding, or 42/35 C DNA-off-target thresholds.

Up to two shortest locally compatible lengths per start per length band are retained.
Separate detailed-QC shortlists now retain up to 400 preferred and 200 extended
candidates per region, within the unchanged 5,000-candidate run limit. Fixed 1 kb
sectors relative to the region start balance these budgets, with spaced starts
within each sector. A 20 kb region has twenty sectors; a short final sector is
supported. Budgets redistribute when sectors lack candidates. Sectors do not alter
labeling-region boundaries or readouts. This is not every possible sequence;
larger shortlists can change final yield and runtime.
No longer arm inherits an earlier candidate's QC. Exact matches across the supplied
references are enumerated; incompatible cross-bin or opposite-strand copies reject
a candidate. Useful noncompeting copies inside the intended bin can be retained.

Actual cheap generation on the bundled ERRFI1 reference produced 600 unique
candidates in 85.467 seconds: 400 preferred and 200 extended, with 20+10 per sector.
All were Top (+), with strict Tm range 45.01-54.98 C. Pre-shortlist generation/filter
counts matched saved v1.5 counts and reference hashes were unchanged. This check
ran no BLAST, full folding or final selection: 600 candidates are not 600 full-QC
passes or a new final yield. See verification/errfi1_shortlist_validation.json.

## 4. Revised RNA policy

DEFAULT: transcript-sequence screening, `avoid_gene_rna=false`, `screen_rna=true`.
Loaded annotations remain visible and saved, but whole gene/transcript spans do
not automatically mask candidates. No inferred whole-span/unspliced local RNA is
generated by that annotation layer. Explicit user RNA inputs are searched locally;
full mode additionally requires the public annotated-human-RNA BLAST database.

RNA-sense references count complementary hits, not same-sense identity. Strong RNA
hits have identity >=90% AND (aligned length >=30 nt OR query coverage >=90%). The
maximum strong RNA-risk cluster count is zero by default. An exact 18-base RNA
match alone is review evidence/ranking input, not an automatic rejection in the
relaxed policy. Short matches are NOT established as harmless. The models do not
predict RNA-DNA mismatch Tm, expression or RNA-binding frequency.

STRICT OPTION: `avoid_gene_rna=true`, retaining the RNA search. Full annotated and
confirmed transcript/gene spans, including introns, constrain candidates on the
chosen strand; conflicting/unknown directions cause rejection, never strand switching.
Those spans also generate local potential unspliced-RNA references; available explicit
annotated spliced features may be included. Exact complementary RNA tracts >=18 nt
also count as a hard RNA risk under this mode. Raw/construct direction declarations
must be confirmed; missing mapped annotations do not count as an RNA-free result.

Strict mode is intentionally conservative. In default mode, unspliced, antisense,
novel and unprovided construct RNAs remain additional uncertainty. A public database
does not contain a custom reporter automatically. Supply explicit transcript FASTA
for sequence screening. No cell-specific expression or RNase cleanup is assumed.
RNA controls are independent of DNA on-target exemptions. Applying the active
ERRFI1 preset restores all policy controls, including Top (+) and the verified RNA
policy, while keeping target/region/readout layout. Full ERRFI1 restore also resets
the bundled reference, exact window, DP1 and standard output. Legacy comparison
presets are documentation/import history, not active starting presets. Imported
settings keep saved values until explicitly reset. Results record actual policy.

## 5. Sequence specificity and intended human sites

Small local references are searched for full-arm ungapped matches with up to three
substitutions and exact full-primary partial tracts >=18 nt, including handle/arm
junctions. Full mode also runs local blastn-short on arms AND complete primaries
against the installed human DNA and requested RNA databases. Commands, counts,
query/output hashes and heuristic limits are recorded. Full audit retains raw
outputs compressed; standard output discards them after parsing. Scientific BLAST
settings/hit-retention criteria are unchanged. A 100,000-subject limit remains finite;
BLAST is heuristic, not an exhaustive proof of hybridization absence.

Intended sites are exempted by bin, reference, coordinate and physical strand, not
by being human. Other bins remain distinct targets. The advanced human-DNA report-only
override preserves evidence but removes human-DNA gates/penalties and is labeled.
It does not exempt RNA, local cross-bin, folding or spacing failures. Mapped targets
must recover the expected perfect genomic arm alignment from the supplied database.

Returned exact-18-nt or sufficiently long/high-identity alignments are review evidence.
For DNA, the ACTUAL aligned query/target sequence is passed to the NN model with
supported mismatches, never simply assigned the on-target Tm. Unsupported/gapped
energetics remain unknown; exact-subsegment proxies are identified separately.
Overlapping HSPs at the same source/reference/strand form clusters. Alternative
contigs and transcript isoforms may still represent redundant physical biology;
counts are searched-coordinate clusters, not experimentally independent loci.

Retained default DNA gates: zero unintended clusters with estimated risk Tm >=42 C;
at most ten moderate/high clusters with risk Tm >=35 C; reject long/high-identity
alignments with unknown energetics: identity >=90% AND (aligned length >=30 nt
OR query coverage >=90%).
Risk Tm is the maximum available whole-alignment or exact-subsegment estimate.
These are project heuristics, not validated FISH failure boundaries or OligoMiner
LDA scores. Spatial search ranking changed; these acceptance comparisons did not.

## 6. Complete-primary folding and joint selection

Full mode requires ViennaRNA DNA Mathews 2004 parameters, 42 C, 0.390 M monovalent
salt, dangles=2, no GU wobble and no G-quadruplex mode. Complete primaries are folded
without truncation. Monomer MFE, ensemble free energy and base-unpaired probabilities
are calculated. Mean arm-unpaired probability <0.40 rejects; equality is allowed.
Self-pairs and candidate-to-selected full-primary pairs are checked before acceptance.
An exact complementary stem >11 nt or cofold MFE gain below -12 kcal/mol with
intermolecular pairs is a default exclusion. Stem 11 and gain exactly -12 are
allowed by these comparisons. All final selected pairs are rechecked.
These are aqueous heuristics: no validated formamide folding model, pseudoknot model,
whole-pool equilibrium or target competition is included. Energy gains are not
concentration-dependent dimer fractions; unpaired probability is not binding yield.

Default selection compares a short-only feasible baseline, additions to that
baseline, joint site-packing and spatial-balance greedy searches. The highest exact
intended-site count among attempted valid solutions is retained. Ties consider
the smallest largest uncovered gap, sum of regional largest gaps and sector
unevenness, before filled-bin/minimum-bin counts, oligo count and arm-base total.
The spatial pass rewards sparse sectors and gap filling. Bounded repair tries
QC-passing alternatives, removes at most two incompatible selected primaries and
refills under every gate: at most two rounds, 128 seed checks and 32 refills. Only
a strict improvement of the count-first objective is accepted, so repair cannot
reduce exact-site count. Final self-pair, pair and spacing checks remain mandatory.

Target 100 and minimum 50 are verified per-region reporting benchmarks in maximum-
site mode, not stopping caps. The default remains four BLAST threads. Boundaries
and QC are not relaxed to fill a region. The finite greedy/repair search has no
global-optimum certificate and does not guarantee more final probes.

Reports show 1 kb sector counts, largest uncovered arm-coordinate gap including
region edges, and median exact-site-center spacing. Sector counts measure exact
binding sites; repeated sites can make these differ from oligo counts. Gap metrics
are distinct from the enforced extended-footprint spacing gate.

Five-base minimum gaps apply between modeled full footprints, extended into TT/
handle bases only while contiguous exact matching can continue. All exact copies
and boundaries are checked. This does not model every partial/mismatched occupancy
state or physically guarantee five unpaired bases in every conformation.

## 7. Installation, outputs and source scope

Place the complete v1.6 folder beside installed RASERFISH_Design_Studio_v1_5.
If v1.6 has no local .qc_env, `bash start_studio.sh` automatically reuses its sibling
runtime and human DNA/RNA databases. An alternative is
`bash start_studio.sh --qc-root "/path/to/installed/Studio"`. The installed folder
must remain at its original location and have its reference manifest. The launcher
executes v1.6 source on port 8769 and saves new results in v1.6/design_runs.
Startup does not install, download, overwrite previous runs or stop old jobs.
A separate first-time install.sh route remains for users without a ready reusable
installation; it installs into v1.6 and records real engine/reference checks.

Standard output is the default. run_config.json holds a small reproducible manifest,
source metadata, mapping, effective policy, hashes and compressed-asset references.
Designed slices are saved in inputs/target_region.fa.gz. Full loaded context is
retained where needed by local exact-site QC; annotations are externalized once.
Keep the manifest and inputs/ together. Reload verifies hashes and restores the
original sequence context/coordinates. External databases are not copied into runs.

candidate_decisions.csv.gz records every shortlisted candidate and its rejection/
status information. design_result.json and the viewer contain selected candidates,
QC/spatial summaries and aggregate rejection counts without repeated full reference
sequences or alignment collections. Selected probes/sites, draft order, strand
JSON/CSV, generation counts and selection-search reports remain available.
Pair summaries count all selected self/pairs; flagged-pair CSV includes failures
and near-threshold warnings. The display-only margins (one stem base, 1 kcal/mol
cofold gain, 0.05 arm-unpaired probability) never change scientific acceptance.

Standard output discards raw BLAST after parsing. Full audit, selected BEFORE a
run, additionally saves compressed raw search/query files and detailed candidate,
alignment and selected-pair evidence under debug/. Both modes perform identical
scientific analysis. Temporary BLAST files and in-memory evidence can still be
large during computation. Full audit can remain substantially larger; a checksum
cannot recover an omitted raw artifact afterwards.

Settings export creates a portable ZIP with its manifest and compressed input
assets; import requires the complete ZIP. Legacy settings JSON with embedded input
sequence remains supported. Review database paths when moving settings between
installations. Imported background/RNA records retain their metadata.

RNA policy is recorded separately from completed analyses. Missing analyses are
not passes. The synthetic local tutorial is marked PREVIEW. Historical CFTR order/
input copies are unchanged. No newly full-QC-screened ERRFI1 order is supplied,
and there is no new final-yield claim. Prior 41-probe tables are baseline evidence;
final JSON/HTML absence is recorded separately and does not establish complete
prior export. Tests freeze retained scientific implementations and verify defaults,
reference hashes and settings roundtrips against the saved evidence.

Sources for inherited algorithms/settings are recorded in SOURCES.md and source
code. This release makes the explicit user-requested policy and packaging changes;
it is not an experimental validation or exact reproduction of all published
LoopTrace/OligoMiner components. LDA and whole-genome k-mer filtering remain absent.
