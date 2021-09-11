import hashlib
import unittest
import os
import sys
testdir = os.path.dirname(__file__)
sys.path.insert(0, os.path.abspath(os.path.dirname(testdir)))
from covidbio.mash_utils import make_nr


class TestCOVID(unittest.TestCase):
    '''
    Input fasta has seven entries, output has six entries
    '''
    outfile = 'cysprot-nr.fa'
    
    def test_make_nr(self):
        infile = os.path.join(testdir, 'cysprot.fa')
        make_nr(infile, threshold=0.15)
        with open(self.outfile) as f:
            num = len([l for l in f if l.startswith(">")])
        self.assertEqual(num, 6, "Correct number of entries")
        self.remove_files()

    def remove_files(self):
        os.remove(self.outfile)

if __name__ == '__main__':
    unittest.main()
