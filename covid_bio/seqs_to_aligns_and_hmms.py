#!/usr/bin/env python3

import argparse
import sys
import os
import re
import tempfile
import subprocess
from Bio import AlignIO
from Bio import SeqIO
from myconfig import COV_DIR

parser = argparse.ArgumentParser()
parser.add_argument('-verbose', default=False, action='store_true', help="Verbose")
parser.add_argument('-aligner', default='mafft', help="Alignment application")
parser.add_argument('-hmmbuilder', default='hmmbuild', help="HMM build application")
parser.add_argument('-skip', default='ORF1a-aa,ORF1a-nt,ORF1ab-aa,ORF1ab-nt', help="Do not align")
parser.add_argument('-json', action='store_true', help="Create JSON for Gen3")
parser.add_argument('-maf', action='store_true', help="Create additional MAF format alignments")
parser.add_argument('-cov_dir', default=COV_DIR, help="Destination directory")
parser.add_argument('files', nargs='+', help='File names')
args = parser.parse_args()


def main():
    builder = Seqs_To_Aligns_And_Hmms(args.verbose, args.aligner, args.hmmbuilder, args.skip, 
        args.json, args.maf, args.files)
    builder.read()
    builder.make_align()
    builder.make_hmm()
    builder.write_json()


class Seqs_To_Aligns_And_Hmms:
    def __init__(self, verbose, aligner, hmmbuild, skip, json, maf, files):
        self.verbose = verbose
        self.aligner = aligner
        self.hmmbuild = hmmbuild
        self.skip = skip
        self.make_json = json
        self.maf = maf
        self.files = files
        self.seqs, self.alns, self.hmms = dict(), dict(), dict()

    def read(self):
        for path in [f for f in self.files if os.path.isfile(f)]:
            seqs = []
            # Get basename without suffix, for example "S-aa"
            name = os.path.basename(path).split('.')[0]
            if name in self.skip:
                continue
            if 'invalid' in name:
                continue
            try:
                if self.verbose:
                    print("Reading Fasta file: {}".format(path))
                # Each one of these is a multiple fasta file
                for fa in enumerate(SeqIO.parse(path, "fasta")):
                    seqs.append(fa)
                self.seqs[name] = self.remove_dups(seqs)
            except (RuntimeError) as exception:
                print("Error parsing sequences in '" +
                      str(path) + "':" + str(exception))

    def remove_dups(self, seqs):
        d = dict()
        for seq in seqs:
            # Do not put sequences containing "X" in alignments
            if 'X' in str(seq.seq):
                continue
            d[str(seq.seq)] = seq
        return list(d.values())

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
                    print("Error running '{}': ".format(self.aligner) + str(exception))

            # Create additional Maf format alignment
            if self.maf:
                AlignIO.convert(name + '.fasta', 'fasta', name + '.maf', 'maf', alphabet=None)

            self.alns[name] = name + '.fasta'

    def make_align_cmd(self, infile, name):
        '''
        time mafft --auto M-aa.fasta > M-aa.aln
            real	0m2.254s
            user	0m2.036s
            sys	0m0.170s
        time muscle -in M-aa.fasta -out M-aa.aln
            real	0m19.963s
            user	0m19.594s
            sys	0m0.238s
        time clustalo -i M-aa.fasta -o x.aln --outfmt=fasta
            real	0m23.045s
            user	0m18.665s
            sys	0m4.255s
        '''
        if self.aligner == 'muscle':
            return [self.aligner, '-quiet','-in', infile, '-out', name + '.fasta'], None
        elif self.aligner == 'mafft':
            return [self.aligner, '--auto', infile], name + '.fasta'
        elif self.aligner == 'clustalo':
            return [self.aligner, '-i', infile, '-o', name + '.fasta', '--outfmt=fasta'], None
        else:
            sys.exit("No command for aligner {}".format(self.aligner))

    def make_hmm(self):
        for name in self.alns:
            if self.verbose:
                print("{0} input file is '{1}'".format(self.hmmbuild, self.alns[name]))
            # Either --amino or --dna
            opt = '--amino' if '-aa' in name else '--dna'
            try:
                subprocess.run([self.hmmbuild, opt, name + '.hmm', self.alns[name]])
            except (subprocess.CalledProcessError) as exception:
                print("Error running {}: ".format(self.hmmbuild) + str(exception))

            self.hmms[name] = name + '.hmm'

    def write_json(self):
        if not self.make_json:
            return
        from make_json import make_hmm_json
        for name in self.hmms:
            json = make_hmm_json(self.hmms[name])
            with open(os.path.join(COV_DIR, name + '-hmm.json'), 'w') as out:
                out.write(json)
        from make_json import make_alignment_json
        for name in self.alns:
            json = make_alignment_json(self.alns[name], self.aligner)
            with open(os.path.join(COV_DIR, name + '-aln.json'), 'w') as out:
                out.write(json)

if __name__ == "__main__":
    main()
