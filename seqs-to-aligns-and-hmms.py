#!/usr/local/bin/python3

import argparse
import os.path
import tempfile
import subprocess
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('-verbose', default=False, action='store_true', help="Verbose")
parser.add_argument('-aligner', default='muscle', help="Alignment application")
parser.add_argument('-hmmbuilder', default='hmmbuild', help="HMM build application")
parser.add_argument('files', nargs='+', help='File names')
args = parser.parse_args()


def main():
    builder = Seqs_To_Aligns_And_Hmms(args.verbose, args.aligner, args.hmmbuilder, args.files)
    builder.read()
    builder.make_align()
    builder.make_hmm()

class Seqs_To_Aligns_And_Hmms:

    def __init__(self, verbose, aligner, hmmbuild, files):
        self.verbose = verbose
        self.aligner = aligner
        self.hmmbuild = hmmbuild
        self.files = files
        self.seqs, self.alns, self.hmms = dict(), dict(), dict()

    def read(self):
        full_paths = [os.path.join(os.getcwd(), path) for path in self.files]
        for path in full_paths:
            seqs = []
            # For example "S-aa"
            name = os.path.basename(path).split('.')[0]
            try:
                # Each one of these is a multiple fasta file
                for index, fa in enumerate(SeqIO.parse(path,"fasta")):
                    seqs.append(fa)
                self.seqs[name] = self.remove_dups(seqs)
            except (RuntimeError) as exception:
                print("Error parsing sequences in '" +
                      str(path) + "':" + str(exception))

    def remove_dups(self, seqs):
        d = dict()
        for seq in seqs:
            d[str(seq.seq)] = seq
        return list(d.values())
        
    def make_align(self):
        for name in self.seqs:
            seqfile = tempfile.NamedTemporaryFile('w',delete=False)
            print(seqfile.name)
            SeqIO.write(self.seqs[name], seqfile, 'fasta')
            cmd = self.make_align_cmd(seqfile.name, name)
            subprocess.run(cmd)
            self.alns[name] = name + '.aln'

    def make_align_cmd(self, infile, name):
        if self.aligner == 'muscle':
            cmd = [self.aligner, '-in', infile, '-out', name + '.aln']
        if self.aligner == 'mafft':
            cmd = [self.aligner, infile, '>', name + '.aln']
        return cmd

    def make_hmm(self):
        for name in self.alns:
            # Either --amino or --dna
            opt = '--amino' if '-aa' in name else '--dna'
            subprocess.run([self.hmmbuild, opt, name + '.hmm', self.alns[name]])
            self.hmms[name] = name + '.hmm'

if __name__ == "__main__":
    main()
