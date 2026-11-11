"""Offline export integration tests; no QC engine or new biological design run."""
import copy
import json
from pathlib import Path
import re
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import studio

class BindingViewerGeneration(unittest.TestCase):
    def setUp(self):
        self.result=studio.load_demo()
    def boot(self, html):
        return json.loads(re.search(r'window.STUDIO_BOOT=(.*?);</script>',html,re.S).group(1))
    def test_new_renderer_embedded(self):
        html=studio.standalone(self.result,True)
        self.assertIn('root.ImagerBinding = api',html)
        self.assertNotIn('src="/web/imager_binding.js"',html)
        self.assertIn('id="assignedBindingPreview"',html)
    def test_direct_guide_initial_step_three(self):
        html=studio.standalone(self.result,True,3)
        self.assertEqual(self.boot(html)['initial_step'],3)
        self.assertTrue(self.boot(html)['preview'])
    def test_bad_initial_step_rejected(self):
        with self.assertRaises(ValueError):studio.standalone(self.result,True,99)
    def test_result_not_mutated(self):
        previous=copy.deepcopy(self.result)
        studio.standalone(self.result,False)
        self.assertEqual(self.result,previous)
    def test_saved_data_and_version_preserved(self):
        boot=self.boot(studio.standalone(self.result,False))
        self.assertEqual(boot['version'],'1.6.1')
        self.assertEqual(boot['demo'],self.result)
        self.assertTrue(boot['offline'])
        self.assertFalse(boot['preview'])

if __name__=='__main__':unittest.main()
