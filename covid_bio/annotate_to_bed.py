#!/usr/bin/env python3

import argparse
import sys
import os
import tempfile
import subprocess
import re
import tmhmm
from Bio import SearchIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from vars import COV_DIR

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
parser.add_argument('-rfam_file', default='cov_allvirus.cm', help="Rfam covariance model file")
parser.add_argument('-cov_dir', default=COV_DIR, help="Destination directory")
parser.add_argument('files', nargs='+', help='File names')
args = parser.parse_args()


def main():
    annotator = Annotate_With_Hmms(args.verbose, args.rfam_file, args.cov_dir, args.files)
    for file in annotator.files:
        gb = annotator.write_fasta(file)
        annotator.find_genes(gb)
        annotator.find_rfam(gb)
        annotator.find_proteins(gb)
        annotator.find_tms(gb)
        annotator.write_bed(gb)

class Annotate_With_Hmms:
    '''
    Create tracks for genes using HMMs ('genes'), Rfam hits ('rfam'), tmhmm predictions ('tms')
    '''
    def __init__(self, verbose, rfam_file, cov_dir, files):
        self.verbose = verbose
        self.rfam_file = rfam_file
        self.cov_dir = cov_dir
        self.files = files
        # COV2 HMMs
        self.cov_proteins = ['ORF1a', 'ORF1ab', 'S', 'E', 'M', 'N',
                            'NS1', 'NS2', 'NS3', 'NS4', 'NS5', 'NS6', 'NS7', 'NS8', 'NS9',
                            'NS10', 'NS11', 'NS12', 'NS13', 'NS14', 'NS15', 'NS16',
                            'ORF3a', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'ORF9b', 'ORF10']
        # Text for BED files
        self.beds = dict()
        # Gene positions on given genome
        self.gene_positions = dict()
        # Protein,gene, and nucleotide TM sequences in given genome
        self.protein_strs = dict()
        self.gene_strs = dict()
        self.tm_positions = dict()

    def write_fasta(self, file):
        '''
        Create genome fasta file using GenBank id
        '''
        gb = SeqIO.read(file, 'gb')
        # Create fasta "reference" sequence
        SeqIO.write(gb, os.path.join(self.cov_dir, gb.id + '.fa'), 'fasta')
        # Index with samtools
        cmd = ['samtools', 'faidx', os.path.join(self.cov_dir, gb.id + '.fa')]
        try:
            subprocess.run(cmd, check=True)
        except (subprocess.CalledProcessError) as exception:
            print("Error: {}".format(exception))
            sys.exit("Error running samtools faidx on {}.fa".format(gb.id))
        return gb

    def find_genes(self, gb):
        self.beds[gb.id] = dict()
        self.beds[gb.id]['genes'] = []
        self.gene_positions[gb.id] = dict()
        trackline = "track name='{0} genes' \
                    description='HMM-based gene detection of COV sequence {1}' \
                    itemRgb='on'".format(gb.id, gb.id)
        self.beds[gb.id]['genes'].append(trackline)
        for protein in self.cov_proteins:
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
            self.gene_positions[gb.id][protein] = (hit.hit_start, hit.hit_end)

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
        self.protein_strs[gb.id] = dict()
        self.gene_strs[gb.id] = dict()
        # Get gene nucleotide and protein sequences
        for protein in self.cov_proteins:
            ntstr = str(gb.seq)[self.gene_positions[gb.id][protein][0]:self.gene_positions[gb.id][protein][1]]
            aastr = self.translate_orf(ntstr, protein, gb)
            if self.verbose:
                print("{0} translation: {1}".format(protein, aastr))
            self.protein_strs[gb.id][protein] = aastr
            self.gene_strs[gb.id][protein] = ntstr

    def find_tms(self, gb):
        '''
        Predict TM regions using tmhmm.py
        '''
        self.beds[gb.id]['tms'] = []
        track_line = "track name='{0} Transmembrame regions' \
                        description='tmhmm-based TM detection of COV sequence {1}' \
                        itemRgb='on'".format(gb.id,gb.id)
        self.beds[gb.id]['tms'].append(track_line)
        for protein in self.protein_strs[gb.id]:
            # Do not analyze polyprotein
            if 'ORF1a' in protein:
                continue
            self.run_tmhmm(self.protein_strs[gb.id][protein], protein, gb)
            if protein in self.tm_positions[gb.id].keys():
                for tm in self.tm_positions[gb.id][protein]:
                    # Get gene coordinates for a given genome
                    feature_line = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}".format(
                        gb.id,
                        self.get_nt_position(tm[0], protein, gb),
                        self.get_nt_position(tm[1], protein, gb),
                        '','','','','','','','','')
                    self.beds[gb.id]['tms'].append(feature_line)

    def get_nt_position(self, num, protein, gb):
        '''
        Return nt position in genome given aa position in a specified protein
        '''
        return (3 * num) + self.gene_positions[gb.id][protein][0]

    def run_tmhmm(self, aastr, protein, gb):
        tm_annotation, _ = tmhmm.predict(aastr)
        if 'M' not in tm_annotation:
            return None
        self.parse_annotation(tm_annotation, protein, gb)

    def parse_annotation(self, annotation, protein, gb):
        # 'm' in 'oooommmmmiiiiimmmmiii': [4, 5, 6, 7, 8, 14, 15, 16, 17]
        # matchs = [i for i, char in enumerate(annotation) if query == char]
        p = re.compile("M+")
        # Matches are zero-based, so the start() of the 'bc' search() in 'abcdef' = 1
        for m in p.finditer(annotation):
            if gb.id not in self.tm_positions.keys():
                self.tm_positions[gb.id] = dict()
            if protein not in self.tm_positions[gb.id].keys():
                self.tm_positions[gb.id][protein] = []
            if self.verbose:
                print("TM: {0} {1} {2}-{3}".format(gb.id, protein, m.start(), m.end()))
            self.tm_positions[gb.id][protein].append((m.start(), m.end()))

    def translate_orf(self, ntstr, protein, gb):
        """
        COV2 ORF1ab contains a -1 frameshift so any translation has to handle this.
        To do: confirm that this approach works for all COV. Possible this method
        has to use an Rfam match rather than simple string search and replace.
        """
        noframeshift = 'TTAAACGGG'
        frameshift =   'TTAAACCGGG'
        aastr = self.translate(ntstr)
        if '*' in aastr:
            if self.verbose:
                print("Stop codon found in {0} {1}".format(gb.id, protein))
            if len(re.findall(noframeshift, ntstr)) == 1:
                if self.verbose:
                    print("Found 1 'slip sequence'")
                ntstr = ntstr.replace(noframeshift, frameshift)
                aastr = self.translate(ntstr)
            else:
                sys.exit("More than 1 'slip sequence' found in {0} {1}".format(gb.id, protein))
        if '*' in aastr:
            sys.exit("Stop codon found in {0} {1}: {2}".format(gb.id, protein, aastr))
        return aastr

    def translate(self, ntstr):
        aaseq = Seq(ntstr, IUPAC.unambiguous_dna).translate(to_stop=False)
        # Remove trailing '*'
        return str(aaseq)[:-1]

    def run_hmmsearch(self, name, hmm):
        '''
        Rum hmmsearch and return the highest scoring hit
        '''
        out = tempfile.NamedTemporaryFile('w')
        hmmpath = os.path.join(self.cov_dir, hmm + '-nt.hmm')
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
        '''
        Run cmscan and return list of all hits
        '''
        out = tempfile.NamedTemporaryFile('w')
        rfampath = os.path.join(self.cov_dir, self.rfam_file)
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
        '''
        Parse cmscan table output and return list of lists
        '''
        with open(file, 'r') as fin:
            return [line.strip().split() for line in fin if not line.startswith('#')]

    def write_bed(self, gb):
        for bed in self.beds[gb.id]:
            with open(os.path.join(self.cov_dir, gb.id + '-' + bed + '.bed'), 'w') as out:
                for line in self.beds[gb.id][bed]:
                    out.write(line + '\n')


if __name__ == "__main__":
    main()
