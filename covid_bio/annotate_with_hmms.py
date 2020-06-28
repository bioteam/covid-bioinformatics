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
Coronavirus reference genome NC_045512 proteins:

NC_045512.2	266	13483	+	43740578	ORF1ab	GU280_gp01	YP_009725295.1	4405	orf1a polyprotein
NC_045512.2	21563	25384	+	43740568	S	GU280_gp02	YP_009724390.1	1273	surface glycoprotein
NC_045512.2	25393	26220	+	43740569	ORF3a	GU280_gp03	YP_009724391.1	275	ORF3a protein
NC_045512.2	26245	26472	+	43740570	E	GU280_gp04	YP_009724392.1	75	envelope protein
NC_045512.2	26523	27191	+	43740571	M	GU280_gp05	YP_009724393.1	222	membrane glycoprotein
NC_045512.2	27202	27387	+	43740572	ORF6	GU280_gp06	YP_009724394.1	61	ORF6 protein
NC_045512.2	27394	27759	+	43740573	ORF7a	GU280_gp07	YP_009724395.1	121	ORF7a protein
NC_045512.2	27756	27887	+	43740574	ORF7b	GU280_gp08	YP_009725318.1	43	ORF7b
NC_045512.2	27894	28259	+	43740577	ORF8	GU280_gp09	YP_009724396.1	121	ORF8 protein
NC_045512.2	28274	29533	+	43740575	N	GU280_gp10	YP_009724397.2	419	nucleocapsid phosphoprotein
NC_045512.2	29558	29674	+	43740576	ORF10



BED Format (https://m.ensembl.org/info/website/upload/bed.html)
BED example:

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
    annotator.annotate()
    annotator.write()


class Annotate_With_Hmms:
    def __init__(self, verbose, hmmdir, files):
        self.verbose = verbose
        self.hmmdir = hmmdir
        self.files = files
        # COV2 HMMs
        self.hmms = ['ORF1a', 'ORF1ab', 'S', 'E', 'M', 'N',
                     'NS1', 'NS2', 'NS3', 'NS4', 'NS5', 'NS6', 'NS7', 'NS8', 'NS9',
                     'NS10', 'NS11', 'NS12', 'NS13', 'NS14', 'NS15', 'NS16',
                     'ORF3a', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'ORF9b', 'ORF10']
        self.bed = dict()


    def annotate(self):
        for file in self.files:
            name = self.write_fasta(file)
            self.bed[name] = []
            trackline = "track name='{0}' description='HMM-based annotation of COV sequence {1}' itemRgb='on'".format(name,name)
            self.bed[name].append(trackline)
            for hmm in self.hmms:
                feat = os.path.basename(hmm).split('-')[0]
                hit = self.run_hmmsearch(name, hmm)
                # Create feature line with thicker line for any ATG
                thickStart = str(hit.hit_start - 1) if 'ORF' in feat or feat in 'EMNS' else ''
                thickEnd = str(hit.hit_start + 2) if 'ORF' in feat or feat in 'EMNS' else ''
                # Color structural proteins, ORFs, and NSPs
                if feat in 'EMNS':
                    color = '66,245,173'
                elif 'ORF' in feat:
                    color = '66,123,245'
                else:
                    color = '245,66,197'

                featureline = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}".format(
                    name,
                    hit.hit_start,
                    hit.hit_end,
                    feat,
                    hit.bitscore,
                    '+',
                    thickStart,
                    thickEnd,
                    color,
                    '',
                    '',
                    '')

                self.bed[name].append(featureline)


    def get_hmms(self):
        # Get all nucleotide HMMs
        # self.hmms = [f for f in glob.glob(self.hmmdir + '/*nt.hmm', recursive=False)]
        if self.verbose:
            print("HMMs: {}".format(self.hmms))


    def write_fasta(self, file):
        gb = SeqIO.read(file, 'gb')
        # Create fasta "reference" sequence
        SeqIO.write(gb, gb.id + '.fa', 'fasta')
        # Index with samtools
        cmd = ['samtools', 'faidx', gb.id + '.fa']
        try:
            subprocess.run(cmd, check=True)
        except (subprocess.CalledProcessError) as exception:
            print("Error: {}".format(exception))
            sys.exit("Error running samtools")
        return gb.id


    def run_hmmsearch(self, name, hmm):
        out = tempfile.NamedTemporaryFile('w')
        hmmpath = self.hmmdir + '/' + hmm + '-nt.hmm'
        cmd = ['hmmsearch', '--noali', '-o', out.name, hmmpath, name + '.fa']
        if self.verbose:
            print("Command: {0}".format(cmd))
        try:
            subprocess.run(cmd, check=True)
        except (subprocess.CalledProcessError) as exception:
            print("Error: {}".format(exception))
            sys.exit("Error running hmmsearch using {}".format(hmm))
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
                    out.write(line + '\n')


if __name__ == "__main__":
    main()
