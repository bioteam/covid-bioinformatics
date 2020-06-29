#!/usr/bin/env python3

import argparse
import sys
import os
import glob
import tempfile
import subprocess
from Bio import SearchIO
from Bio import SeqIO
import tmhmm

'''
Annotate COV GenBank files using a collection of HMMs.

Example BED files showing genes:

track name='NC_045512.2' description='HMM-based annotation of COV sequence NC_045512.2' itemRgb='on'
NC_045512.2	265	13483	ORF1a	17560.1	+	264	267	66,123,245			
NC_045512.2	265	21555	ORF1ab	28346.8	+	264	267	66,123,245			
NC_045512.2	21562	25384	S	4981.5	+	21561	21564	66,245,173			
NC_045512.2	26244	26472	E	271.3	+	26243	26246	66,245,173			
NC_045512.2	26522	27191	M	886.4	+	26521	26524	66,245,173			
NC_045512.2	28273	29533	N	1656.6	+	28272	28275	66,245,173			
NC_045512.2	265	805	NS1	715.7	+			245,66,197			
NC_045512.2	805	2719	NS2	2524.5	+			245,66,197			
...
NC_045512.2	20658	21552	NS16	1104.6	+			245,66,197			
NC_045512.2	25392	26220	ORF3a	1070.6	+	25391	25394	66,123,245			
NC_045512.2	25813	26281	ORF3b	279.1	+	25812	25815	66,123,245			
NC_045512.2	27201	27387	ORF6	210.4	+	27200	27203	66,123,245			
NC_045512.2	27393	27759	ORF7a	463.2	+	27392	27395	66,123,245			
NC_045512.2	27755	27887	ORF7b	140.6	+	27754	27757	66,123,245			
NC_045512.2	27893	28259	ORF8	416.6	+	27892	27895	66,123,245			
NC_045512.2	28283	28577	ORF9b	288.1	+	28282	28285	66,123,245			
NC_045512.2	29557	29674	ORF10	123.3	+	29556	29559	66,123,245	

BED Format (https://m.ensembl.org/info/website/upload/bed.html)
'''

parser = argparse.ArgumentParser()
parser.add_argument('-verbose', default=False, action='store_true', help="Verbose")
parser.add_argument('-hmmdir', default=os.path.dirname(os.path.abspath(__file__)), help="HMM directory")
parser.add_argument('files', nargs='+', help='File names')
args = parser.parse_args()


def main():
    annotator = Annotate_With_Hmms(args.verbose, args.hmmdir, args.files)
    annotator.find_genes()
    annotator.find_tms()
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


    def find_genes(self):
        self.bed['genes'] = dict()
        for file in self.files:
            name = self.write_fasta(file)
            self.bed['genes'][name] = []
            trackline = "track name='{0} genes' \
                        description='HMM-based gene detection of COV sequence {1}' \
                        itemRgb='on'".format(name,name)
            self.bed['genes'][name].append(trackline)
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

                self.bed['genes'][name].append(featureline)


    def find_tms(self):
        '''
        Translate predicted gene sequences and predict their TM regions
        '''
        self.bed['tms'] = dict()
        for file in self.files:
            gb = SeqIO.read(file,'gb')
            self.bed['tms'][gb.id] = []
            trackline = "track name='{0} Transmembrame regions' \
                        description='tmhmm-based TM detection of COV sequence {1}' \
                        itemRgb='on'".format(name,name)
            self.bed['tms'][gb.id].append(trackline)
            for track in self.bed['genes']:



    def run_tmhmm(self):=
        annotation, posterior = tmhmm.predict(sequence, 'TMHMM2.0.model')


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
            sys.exit("Error running samtools faidx")
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
        for tracktype in self.bed:
            for seq in self.bed[tracktype]:
                with open(seq + '.bed', 'w') as out:
                    for line in self.bed[tracktype][seq]:
                        out.write(line + '\n')


if __name__ == "__main__":
    main()
