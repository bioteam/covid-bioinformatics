#!/usr/bin/env python3

import argparse
import sys
import os
import glob
import tempfile
import subprocess
from Bio import SearchIO
from Bio import SeqIO

'''
BED Format (https://m.ensembl.org/info/website/upload/bed.html)

Required fields
The first three fields in each feature line are required:

chrom - name of the chromosome or scaffold. Any valid seq_region_name can be used, with or without the 'chr' prefix.
chromStart - Start position of the feature in standard chromosomal coordinates (i.e. first base is 0).
chromEnd - End position of the feature in standard chromosomal coordinates

Optional fields
Nine additional fields are optional. Note that columns cannot be empty - lower-numbered fields must always be populated if higher-numbered ones are used.

name - Label to be displayed under the feature, if turned on in "Configure this page".
score - A score between 0 and 1000. See track lines, below, for ways to configure the display style of scored data.
strand - defined as + (forward) or - (reverse).
thickStart - coordinate at which to start drawing the feature as a solid rectangle
thickEnd - coordinate at which to stop drawing the feature as a solid rectangle
itemRgb - an RGB colour value (e.g. 0,0,255). Only used if there is a track line with the value of itemRgb set to "on" (case-insensitive).
blockCount - the number of sub-elements (e.g. exons) within the feature
blockSizes - the size of these sub-elements
blockStarts - the start coordinate of each sub-element

Track lines

Track definition lines can be used to configure the display further, e.g. by grouping features into separate tracks. 
Track lines should be placed at the beginning of the list of features they are to affect.

The track line consists of the word 'track' followed by space-separated key=value pairs - see the example below. 
Valid parameters used are:

name - unique name to identify this track when parsing the file
description - Label to be displayed under the track in Region in Detail
priority - integer defining the order in which to display tracks, if multiple tracks are defined.
color - as RGB, hex or X11 named color.
useScore - set to 1 to render the track in greyscale based on the values in the score column.
itemRgb - if set to 'on' (case-insensitive), the individual RGB values defined in tracks will be used.

track name="ItemRGBDemo" description="Item RGB demonstration" itemRgb="On"
chr7  127471196  127472363  Pos1  0  +  127471196  127472363  255,0,0
chr7  127472363  127473530  Pos2  0  +  127472363  127473530  255,0,0
chr7  127475864  127477031  Neg1  0  -  127475864  127477031  0,0,255
chr7  127477031  127478198  Neg2  0  -  127477031  127478198  0,0,255
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
    annotator.write()


class Annotate_With_Hmms:
    def __init__(self, verbose, hmmdir, files):
        self.verbose = verbose
        self.hmmdir = hmmdir
        self.files = files
        self.bed = dict()


    def annotate(self):
        for file in self.files:
            name = self.write_fasta(file)
            self.bed[name] = []
            self.bed[name].append(['track name="' + name + '"',
                'description="HMM-based annotation of COV sequence ' + name + '"',
                'itemRgb="on"'])
            for hmm in self.hmms:
                feat = os.path.basename(hmm).split('-')[0]
                if self.verbose:
                    print("Annotating {0} with {1}".format(name,hmm))
                hit = self.run_hmmsearch(name,hmm)
                # Create feature line
                featureline = []
                featureline.append(str(1))
                featureline.append(str(hit.hit_start - 1))
                featureline.append(str(hit.hit_end))
                featureline.append(feat)
                featureline.append(str(hit.bitscore))
                featureline.append('+')
                # If there is an ATG
                thickStart = str(hit.hit_start - 1) if 'ORF' in feat or feat in 'EMNS' else ''
                featureline.append(thickStart)
                thickEnd = str(hit.hit_start + 2) if 'ORF' in feat or feat in 'EMNS' else ''
                featureline.append(thickEnd)
                # Color features
                if feat in 'EMNS':
                    featureline.append('66,245,173')
                elif 'ORF' in feat:
                    featureline.append('66,123,245')
                else:
                    featureline.append('245,66,197')
                self.bed[name].append(featureline)


    def get_hmms(self):
        # Get all nucleotide HMMs
        self.hmms = [f for f in glob.glob(self.hmmdir + '/*nt.hmm', recursive=False)]
        if self.verbose:
            print("HMMs: {}".format(self.hmms))


    def write_fasta(self, file):
        name = os.path.basename(file).split('.')[0]
        gb = SeqIO.parse(file, 'gb')
        # Create fasta "reference" sequence
        SeqIO.write(gb, name + '.fa', 'fasta')
        # Index with samtools
        cmd = ['samtools', 'faidx', name + '.fa']
        try:
            subprocess.run(cmd, check=True)
        except (subprocess.CalledProcessError) as exception:
            print("Error: {}".format(exception))
            sys.exit("Error running samtools")
        return name


    def run_hmmsearch(self,name,hmm):
        out = tempfile.NamedTemporaryFile('w')
        cmd = ['hmmsearch', '--noali', '-o', out.name, hmm, name + '.fa']
        try:
            subprocess.run(cmd, check=True)
        except (subprocess.CalledProcessError) as exception:
            print("Error: {}".format(exception))
            sys.exit("Error running hmmsearch")
        bestscore = 0
        # Get HSP with highest score
        for qresult in SearchIO.parse(out.name, 'hmmer3-text'):
            for hit in qresult:
                for hsp in hit:
                    if hsp.bitscore > bestscore:
                        besthit = hsp
                        bestscore = hsp.bitscore
        return besthit


    def write(self):
        for seq in self.bed:
            with open(seq + '.bed', 'w') as out:
                for line in self.bed[seq]:
                    out.write('\t'.join(line) + '\n')


if __name__ == "__main__":
    main()
