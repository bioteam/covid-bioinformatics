#!/usr/bin/env python3

import argparse
import sys
import os
import glob
# import re
# import tempfile
import subprocess
from Bio import SearchIO

'''
BED Format (https://genome.ucsc.edu/FAQ/FAQformat.html#format1)

The first three required BED fields are:

chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the 
feature, however, the number in position format will be represented. For example, the first 100 bases of chromosome 1 are defined as 
chrom=1, chromStart=0, chromEnd=100, and span the bases numbered 0-99 in our software (not 0-100), but will represent the position 
notation chr1:1-100.

The 9 additional optional BED fields are:

name - Defines the name of the BED line. 
score - A score between 0 and 1000.
strand - Defines the strand. Either "." (=no strand) or "+" or "-".
thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). When there 
is no thick part, thickStart and thickEnd are usually set to the chromStart position.
thickEnd - The ending position at which the feature is drawn thickly (for example the stop codon in gene displays).
itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will 
determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors 
or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
blockCount - The number of blocks (exons) in the BED line.
blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. 
The number of items in this list should correspond to blockCount.
In BED files with block definitions, the first blockStart value must be 0, so that the first block begins at chromStart. Similarly, 
the final blockStart position plus the final blockSize value must equal chromEnd. Blocks may not overlap.

Here's an example of an annotation track, introduced by a header line:

track name=pairedReads description="Clone Paired Reads" useScore=1
chr22 1000 5000 cloneA 960 + 1000 5000 0 2 567,488, 0,3512
chr22 2000 6000 cloneB 900 - 2000 6000 0 2 433,399, 0,3601
'''

parser = argparse.ArgumentParser()
parser.add_argument('-verbose', default=False, action='store_true', help="Verbose")
parser.add_argument('-hmmdir', default=os.path.dirname(os.path.abspath(__file__)), help="HMM directory")
parser.add_argument('files', nargs='+', help='File names')
args = parser.parse_args()


def main():
    annotator = Annotate_With_Hmms(args.verbose, args.hmmdir, args.files)
    annotator.get_hmms()
    annotator.annotate()


class Annotate_With_Hmms:
    def __init__(self, verbose, hmmdir, files):
        self.verbose = verbose
        self.hmmdir = hmmdir
        self.files = files

    def annotate(self):
        for file in self.files:
            name = os.path.basename(file).split('.')[0]
            if self.verbose:
                print("Annotating {}".format(name))

    def get_hmms(self):
        self.hmms = [f for f in glob.glob(
            self.hmmdir + '/*nt.hmm', recursive=False)]
        if self.verbose:
            print("HMMs: {}".format(self.hmms))

    def read(self):
        for path in [f for f in self.files if os.path.isfile(f)]:
            seqs = []
            # Get basename without suffix, for example "S-aa"
            name = os.path.basename(path).split('.')[0]
            try:
                if self.verbose:
                    print("Reading Fasta file: {}".format(path))
                # Each one of these is a multiple fasta file
                for index, fa in enumerate(SeqIO.parse(path, "fasta")):
                    seqs.append(fa)
                self.seqs[name] = self.remove_dups(seqs)
            except (RuntimeError) as exception:
                print("Error parsing sequences in '" +
                      str(path) + "':" + str(exception))

    def make_align(self):
        for name in self.seqs:
            # Most aligners will reject a file with a single sequence so just copy
            if len(self.seqs[name]) == 1:
                cmd = ['cp', name + '.fa', name + '.fasta']
                subprocess.run(cmd, check=True)
            else:
                seqfile = tempfile.NamedTemporaryFile('w', delete=False)
                SeqIO.write(self.seqs[name], seqfile.name, 'fasta')
                if self.verbose:
                    print("Alignment input sequence file: {}".format(seqfile.name))
                # out_filename is used to redirect the STDOUT to file
                # when self.aligner is "mafft" and it requires redirect to file
                cmd, out_filename = self.make_align_cmd(seqfile.name, name)
                if self.verbose:
                    print("Alignment command is '{}'".format(cmd))
                try:
                    if out_filename:
                        with open(out_filename, "w") as f:
                            subprocess.run(cmd, check=True, stdout=f)
                    else:
                        subprocess.run(cmd, check=True)
                except (subprocess.CalledProcessError) as exception:
                    print("Error running '{}': ".format(
                        self.aligner) + str(exception))

            # Create additional Maf format alignment
            if self.maf:
                AlignIO.convert(name + '.fasta', 'fasta',
                                name + '.maf', 'maf', alphabet=None)

            self.alns[name] = name + '.fasta'

    def make_align_cmd(self, infile, name):
        '''
        '''
        if self.aligner == 'muscle':
            return [self.aligner, '-quiet', '-in', infile, '-out', name + '.fasta'], None
        elif self.aligner == 'mafft':
            return [self.aligner, '--auto', infile], name + '.fasta'
        elif self.aligner == 'clustalo':
            return [self.aligner, '-i', infile, '-o', name + '.fasta', '--outfmt=fasta'], None
        else:
            sys.exit("No command for aligner {}".format(self.aligner))

    def make_hmm(self):
        for name in self.alns:
            if self.verbose:
                print("{0} input file is '{1}'".format(
                    self.hmmbuild, self.alns[name]))
            # Either --amino or --dna
            opt = '--amino' if '-aa' in name else '--dna'
            try:
                subprocess.run(
                    [self.hmmbuild, opt, name + '.hmm', self.alns[name]])
            except (subprocess.CalledProcessError) as exception:
                print("Error running {}: ".format(
                    self.hmmbuild) + str(exception))

            self.hmms[name] = name + '.hmm'


if __name__ == "__main__":
    main()
