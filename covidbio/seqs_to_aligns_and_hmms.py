#!/usr/bin/env python3

import argparse
import os
import re
import subprocess
import sys
import tempfile
from Bio import AlignIO
from Bio import SeqIO
from covidbio.utilities import read_config

config = read_config()

parser = argparse.ArgumentParser()
parser.add_argument('-verbose', default=False, action='store_true', help="Verbose")
parser.add_argument('-aligner', default='clustalo', help="Alignment application")
parser.add_argument('-hmmbuilder', default='hmmbuild', help="HMM build application")
parser.add_argument('-skip', default=['invalid'], nargs='+', help="Do not align if filename contains string")
parser.add_argument('-json', action='store_true', help="Create JSON for Gen3")
parser.add_argument('-maf', action='store_true', help="Create additional MAF format alignments")
parser.add_argument('-strain', help="Strain name", default=config['STRAIN'])
parser.add_argument('-data_dir', help="Location for strain-specific directories", default=config['DATA_DIR'])
parser.add_argument('files', nargs='+', help='File names')
args = parser.parse_args()

'''

'''

def main():
    builder = Seqs_To_Aligns_And_Hmms(args.verbose, args.aligner, args.hmmbuilder, args.skip, 
        args.json, args.maf, args.strain, args.data_dir, args.files)
    for f in builder.files:
        if not self.is_file(f):
            continue
        name = os.path.basename(f).split('.')[0]
        if self.should_skip(name):
            continue
        builder.make_align(f, name)
        builder.make_maf(name)
        builder.make_hmm(name)
        builder.write_json(name)


class Seqs_To_Aligns_And_Hmms:
    def __init__(self, verbose, aligner, hmmbuild, skip, json, maf, strain, data_dir, files):
        self.strain = strain
        self.data_dir = data_dir
        self.verbose = verbose
        self.aligner = aligner
        self.hmmbuild = hmmbuild
        self.skip = skip
        self.make_json = json
        self.maf = maf
        self.cov_dir = os.path.join(self.data_dir, self.strain)
        self.files = files
        if not os.path.isdir(self.cov_dir):
            sys.exit("Directory {} does not exist".format(self.cov_dir))

    '''
    Make a list of unique sequences, create temporary input file and
    make the alignment (*.fasta)
    '''
    def make_align(self, f, name):
        align_name = os.path.join(self.cov_dir, name + '.fasta')
        if self.file_exists(align_name):
            return
        seqs = self.read(f)
        # Most aligners will reject a file with a single sequence so
        # duplicate the *fa file to make the alignment *fasta file
        if len(seqs) == 1:
            cmd = ['cp', 
                    os.path.join(self.cov_dir, name + '.fa'), 
                    align_name]
            subprocess.run(cmd, check=True)
        else:
            seqs = self.remove_dups(seqs, name)
            tmpfasta = tempfile.NamedTemporaryFile()
            SeqIO.write(seqs, tmpfasta.name, 'fasta')
            # out_filename is used to redirect the STDOUT to file
            # when self.aligner is "mafft" and it requires redirect to file
            cmd, out_filename = self.make_align_cmd(tmpfasta.name, align_name)
            if self.verbose:
                print("{0} command: '{1}'".format(self.aligner, cmd))
            try:
                if out_filename:
                    with open(out_filename, "w") as f:
                        subprocess.run(cmd, check=True, stdout=f)
                else:
                    subprocess.run(cmd, check=True)
            except (subprocess.CalledProcessError) as exception:
                print("Error running '{}': ".format(cmd) + str(exception))

    '''
    Read one *fa file, return array of unique sequences.
    '''
    def read(self, path):
        try:
            if self.verbose:
                print("Reading Fasta file: {}".format(path))
            seqs = list(SeqIO.parse(path, "fasta")
        except (RuntimeError) as exception: 
            print("Error parsing sequences in '" +
                str(path) + "':" + str(exception))
        return seqs

    def remove_dups(self, records, filename):
        d = dict()
        # Do not put sequences containing 'X', 'N', or '*' in alignments
        if '-aa' in filename:
            invalid = re.compile(r'[xX*]')
        elif '-nt' in filename:
            invalid = re.compile(r'[xXnN*]')
        else:
            sys.exit("Cannot determine sequence type in {}".format(filename))
        for record in records:
            seq_string = str(record.seq)
            if re.search(invalid, seq_string):
                continue
            d[seq_string] = record
        return list(d.values())

    def file_exists(self, f):
        if os.path.exists(f) and os.stat(f).st_size > 0:
            if self.verbose:
                print("File {} already exists".format(f))
            return True

    def should_skip(self, f):
        for str in builder.skip:
            if str in f:
                if self.verbose:
                    print('Skipping {0}'.format(f))
                return True

    def is_file(self, f):
        if not os.path.isfile(f):
            if self.verbose:
                print('{0} is not a file'.format(f))
            return False

    def make_maf(self, name):
        '''
        Create additional Maf format alignment if needed
        '''
        if self.maf:
            AlignIO.convert(os.path.join(self.cov_dir, name + '.fasta'),
                            'fasta',
                            os.path.join(self.cov_dir, name + '.maf'),
                            'maf')

    def make_align_cmd(self, infile, align_name):
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
            return [self.aligner, '-quiet','-in', infile, '-out', align_name], None
        elif self.aligner == 'mafft':
            return [self.aligner, '--quiet', '--auto', infile], align_name
        elif self.aligner == 'clustalo':
            return [self.aligner, '-i', infile, '-o', align_name, '--outfmt=fasta', '--auto'], None
        else:
            sys.exit("Cannot make command for unknown aligner {}".format(self.aligner))

    def make_hmm(self, name):
        hmm_name = os.path.join(self.cov_dir, name + '.hmm')
        if self.file_exists(hmm_name):
            return
        # Either --amino or --dna
        opt = '--amino' if '-aa' in name else '--dna'
        cmd = [self.hmmbuild, 
                opt, 
                hmm_name, 
                os.path.join(self.cov_dir, name + '.fasta')]
        if self.verbose:
            print("{0} command: {1}".format(self.hmmbuild, cmd))
        try:
            subprocess.run(cmd, check=True)
        except (subprocess.CalledProcessError) as exception:
            print("Error running {}: ".format(self.hmmbuild) + str(exception))

    def write_json(self, name):
        if not self.make_json:
            return
        from utilities import make_hmm_json
        json = make_hmm_json(name)
        with open(os.path.join(self.cov_dir, name + '-hmm.json'), 'w') as out:
            out.write(json)
        from utilities import make_alignment_json
        json = make_alignment_json(name, self.aligner)
        with open(os.path.join(self.cov_dir, name + '-aln.json'), 'w') as out:
            out.write(json)

if __name__ == "__main__":
    main()