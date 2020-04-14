#!/usr/local/bin/python3

import argparse
import os.path
import yaml
import tempfile
from Bio import SeqIO
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from Bio.Alphabet import IUPAC

parser = argparse.ArgumentParser()
# parser.add_argument('-f', default='fasta', dest='format', help="Output format")
# parser.add_argument('-split', default=False, dest='split', help="Split into separate files")
# parser.add_argument('-no-split', action='store_false', dest='split', help="Create one file")
parser.add_argument('-analyze', default=False, action='store_true', help="Parse files without file creation")
parser.add_argument('-verbose', default=False, action='store_true', help="Verbose")
parser.add_argument('-aligner', default='muscle', help="Alignment application")
parser.add_argument('files', nargs='+', help='File names')
args = parser.parse_args()


def main():
    builder = Seqs_To_Aligns_And_Hmms(args.analyze, args.verbose, args.aligner, args.files)
    builder.read()
    builder.align()
    # builder.create_objects()
    # builder.write()


class Seqs_To_Aligns_And_Hmms:

    def __init__(self, analyze, verbose, aligner, files):
        # self.seq_format = seq_format
        self.analyze = analyze
        self.verbose = verbose
        self.aligner = aligner
        self.files = files
        self.seqs = dict()

    def read(self):
        full_paths = [os.path.join(os.getcwd(), path) for path in self.files]
        for path in full_paths:
            seqs = []
            try:
                # Each one of these is a multiple fasta file
                for index, fa in enumerate(SeqIO.parse(path,"fasta")):
                    seqs.append(fa)
                seqs = self.remove_dups(seqs)
                self.seqs[os.path.basename(path)] = seqs
            except (RuntimeError) as exception:
                print("Error parsing sequences in '" +
                      str(path) + "':" + str(exception))

    def remove_dups(self, seqs):
        d = dict()
        for seq in seqs:
            d[str(seq.seq)] = seq
        return list(d.values())
        
    def align(self):
        for filename in self.seqs:
            # print("{}".format(self.seqs[filename]))
            fp = tempfile.TemporaryFile()
            SeqIO.write(self.seqs[filename], fp, 'fasta')
   
        # self.get_host()
        # self.get_date()
        # self.get_organism()
        # self.standardize_cds()
        # self.standardize_mat_peptide()

    def create_objects(self):
        if self.analyze:
            return
        '''
        '''


    def write(self):
        if self.analyze:
            return
        if self.split:
            return
            # for acc in self.accs:
            #     for feat in self.accs[acc]:
            #         seqfile = acc + '.' + self.seq_format
            # SeqIO.write(record, seqfile, self.seq_format)
            # seqfile = record + '-CDS' + '.' + self.seq_format
            # SeqIO.write(seqs, seqfile, self.seq_format)
        else:
            for name in self.aa.keys():
                seqfile = name + '-aa.' + self.seq_format
                SeqIO.write(self.aa[name], seqfile, self.seq_format)
            for name in self.nt.keys():
                seqfile = name + '-nt.' + self.seq_format
                SeqIO.write(self.nt[name], seqfile, self.seq_format)

if __name__ == "__main__":
    main()
