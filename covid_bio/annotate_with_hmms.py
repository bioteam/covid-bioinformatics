#!/usr/bin/env python3

import argparse
import sys
import os
import glob
import tempfile
import subprocess
import re
from Bio import SearchIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import tmhmm

'''
Annotate COV GenBank files using a collection of HMMs. Example BED files showing genes:

track name='NC_045512.2' description='HMM-based annotation of COV sequence NC_045512.2' itemRgb='on'
NC_045512.2	265	13483	ORF1a	17560.1	+	264	267	66,123,245			
...
NC_045512.2	29557	29674	ORF10	123.3	+	29556	29559	66,123,245	

BED Format (https://m.ensembl.org/info/website/upload/bed.html)
'''

parser = argparse.ArgumentParser()
parser.add_argument('-verbose', action='store_true', help="Verbose")
parser.add_argument('-hmmdir', default=os.path.dirname(os.path.abspath(__file__)), help="HMM directory")
parser.add_argument('-rfamfile', default='cov_allvirus.cm', help="Rfam covariance model file")
parser.add_argument('files', nargs='+', help='File names')
args = parser.parse_args()


def main():
    annotator = Annotate_With_Hmms(args.verbose, args.hmmdir, args.rfamfile, args.files)
    for file in annotator.files:
        gb = annotator.write_fasta(file)
        annotator.find_genes(gb)
        annotator.find_rfam(gb)
        annotator.find_proteins(gb)
        # annotator.find_tms(gb)
        annotator.write_bed(gb)


class Annotate_With_Hmms:
    '''
    Create tracks for genes using HMMs ('genes'), Rfam hits ('rfam'), tmhmm predictions ('tms')
    '''
    def __init__(self, verbose, hmmdir, rfamfile, files):
        self.verbose = verbose
        self.hmmdir = hmmdir
        self.rfamfile = rfamfile
        self.files = files
        # COV2 HMMs
        self.covproteins = ['ORF1a', 'ORF1ab', 'S', 'E', 'M', 'N',
                            'NS1', 'NS2', 'NS3', 'NS4', 'NS5', 'NS6', 'NS7', 'NS8', 'NS9',
                            'NS10', 'NS11', 'NS12', 'NS13', 'NS14', 'NS15', 'NS16',
                            'ORF3a', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'ORF9b', 'ORF10']
        # Text for BED files
        self.beds = dict()
        # Gene positions on given genome
        self.genes = dict()
        # Corresponding protein sequences in given genome
        self.proteins = dict()

    def write_fasta(self, file):
        '''
        Create genome fasta file using GenBank id
        '''
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
        return gb


    def find_genes(self, gb):
        self.beds[gb.id] = dict()
        self.beds[gb.id]['genes'] = []
        self.genes[gb.id] = dict()
        trackline = "track name='{0} genes' \
                    description='HMM-based gene detection of COV sequence {1}' \
                    itemRgb='on'".format(gb.id, gb.id)
        self.beds[gb.id]['genes'].append(trackline)
        for protein in self.covproteins:
            hit = self.run_hmmsearch(gb.id, protein)
            # Create feature line with thicker line for any ATG
            thickStart = str(hit.hit_start - 1) if 'ORF' in protein or protein in 'EMNS' else ''
            thickEnd = str(hit.hit_start + 2) if 'ORF' in protein or protein in 'EMNS' else ''
            # Color structural proteins, ORFs, and NSPs
            if protein in 'EMNS':
                color = '66,245,173'
            elif 'ORF' in protein:
                color = '66,123,245'
            else:
                color = '245,66,197'

            featureline = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}".format(
                    gb.id,
                    hit.hit_start,
                    hit.hit_end,
                    protein,
                    hit.bitscore,
                    '+',
                    thickStart,
                    thickEnd,
                    color,
                    '','','')

            self.beds[gb.id]['genes'].append(featureline)
            self.genes[gb.id][protein] = (hit.hit_start, hit.hit_end)


    def find_rfam(self, gb):
        if gb.id not in self.beds.keys():
            self.beds[gb.id] = dict()
        self.beds[gb.id]['rfam'] = []
        trackline = "track name='{0} Rfam hits' \
                    description='COV Rfam hits in sequence {1}' \
                    useScore=1".format(gb.id, gb.id)
        self.beds[gb.id]['rfam'].append(trackline)
        hits = self.run_cmscan(gb.id)
        for hit in hits:
            # ['Sarbecovirus-3UTR', 'RF03125', 'NC_045512.2', '-', 'cm', '1', '335', '29536', '29870',
            #  '+', 'no', '1', '0.40', '0.0', '415.9', '7.4e-127', '!', 'Sarbecovirus', "3'UTR"]
            featureline = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}".format(
                    gb.id,
                    hit[7],
                    hit[8],
                    hit[0],
                    hit[14],
                    hit[9],
                    '','','','','','')
            self.beds[gb.id]['rfam'].append(featureline)


    def find_proteins(self, gb):
        '''
        Translate predicted gene sequences and account for frameshift if necessary
        '''
        # Get gene nucleotide and protein sequences
        for protein in self.covproteins:
            ntstr = str(gb.seq)[self.genes[gb.id][protein][0]:self.genes[gb.id][protein][1]]
            aaseq = self.translate_orf(ntstr, protein, gb)
            # Get gene coordinates for a given genome
            # for num, hmm in enumerate(self.genes[gb.id]):
            #     # Skip ORF1a and ORF1ab
            #     if 'ORF1a' in hmm:
            #         continue
            #     # Skip track line
            #     if num == 0:
            #         continue
            #     # Get gene nucleotide and protein sequences
            #     ntstr = str(gb.seq)[self.genes[gb.id][hmm][0]:self.genes[gb.id][hmm][1]]
            #     aaseq = self.translate_orf(ntstr)
            if self.verbose:
                print("{0} translation: {1}".format(hmm, str(aaseq)))
       

    def find_tms(self, gb):
        '''
        Predict TM regions
        '''
        self.beds[gbid]['tms'] = dict()
        for file in self.files:
            gb = SeqIO.read(file,'gb')
            self.beds[gbid]['tms']= []
            trackline = "track name='{0} Transmembrame regions' \
                        description='tmhmm-based TM detection of COV sequence {1}' \
                        itemRgb='on'".format(gb.id,gb.id)
            self.beds[gbid]['tms'].append(trackline)
            # Get gene coordinates for a given genome
            for num, hmm in enumerate(self.genes[gb.id]):
                # Skip ORF1a and ORF1ab
                if 'ORF1a' in hmm:
                    continue
                # Skip track line
                if num == 0:
                    continue
                if self.verbose:
                    print("{0} translation: {1}".format(hmm, str(aaseq)))
                tms = self.run_tmhmm(str(aaseq))


    def run_tmhmm(self, aaseq):
        annotation, posterior = tmhmm.predict(str(aaseq))
        if 'M' not in annotation:
            return None
        return


    def translate_orf(self, ntstr, protein, gb):
        """
        COV2 ORF1ab contains a -1 frameshift so any translation has to handle this.
        To do: confirm that this approach works for all COV. Possible this method
        has to use an Rfam match rather than simple search and replace.
        """
        noframeshift = 'TTAAACGGG'
        frameshift =   'TTAAACCGGG'
        # Remove trailing *
        aastr = str(aaseq)[:-1]
        if '*' in aaseq:
            if self.verbose:
                print("Stop codon found in {0} {1}".format(gb.id, protein))
            if len(re.findall(noframeshift, ntstr)) == 1:
                if self.verbose:
                    print("Found 1 'slip sequence'")
                nstr.replace(noframeshift, frameshift)
                aaseq = Seq(ntstr, IUPAC.unambiguous_dna).translate(to_stop=False)
            else:
                sys.exit("More than 1 'slip sequence' found in {0} {1}".format(gb.id, protein))
        if '*' in aaseq:
            sys.exit("Problem translating {0} {1}".format(gb.id, protein))
        return aaseq
        

    def run_hmmsearch(self, name, hmm):
        '''
        Rum hmmsearch and return the highest scoring hit
        '''
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


    def run_cmscan(self, name):
        out = tempfile.NamedTemporaryFile('w')
        rfampath = os.path.join(self.hmmdir, self.rfamfile)
        cmd = ['cmscan', '--tblout', out.name, rfampath, name + '.fa']
        if self.verbose:
            print("Command: {0}".format(cmd))
        try:
            subprocess.run(cmd, check=True)
        except (subprocess.CalledProcessError) as exception:
            print("Error: {}".format(exception))
            sys.exit("Error running cmscan using {}".format(name))
        hits = self.parse_cmscan(out.name)
        return hits


    def parse_cmscan(self, file):
        with open(file, 'r') as fin:
            return [line.strip().split() for line in fin if not line.startswith('#')]
        #     lines = [line.strip().split() for line in fin if not line.startswith('#')]
        # return lines


    def write_bed(self, gb):
        for bed in self.beds[gb.id]:
            with open(gb.id + '-' + bed + '.bed', 'w') as out:
                for line in self.beds[gb.id][bed]:
                    out.write(line + '\n')


if __name__ == "__main__":
    main()
