#!/usr/bin/env python3

import argparse
import sys
import os
import re
import glob
import subprocess
import numpy
from Bio import SearchIO
from Bio import Entrez
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('-verbose', default=False, action='store_true', help="Verbose")
parser.add_argument('-download', default=False, action='store_true', help="Download hit sequences")
parser.add_argument('-align', default=False, action='store_true', help="Align hits to HMM")
parser.add_argument('-taxfilter', default=None, help="Exclude clade")
parser.add_argument('-chunk', default=10, help="Number of ids to send to Elink")
parser.add_argument('files', nargs='+', help='File names')
args = parser.parse_args()


def main():
    query = Parse_Hmmsearch(args.verbose, args.download, args.align, args.taxfilter, args.chunk, args.files)
    query.parse()
    query.filter()
    query.download_hits()
    query.align_hits_to_hmm()

class Parse_Hmmsearch:

    def __init__(self, verbose, download, align, taxfilter, chunk, files):
        self.verbose = verbose
        self.download = download
        self.align = align
        self.taxfilter = taxfilter
        self.chunk = chunk
        self.files = files
        # The primary keys for hits{} and fasta{} are file names, the secondary keys are the hits
        self.hits = dict()
        self.fasta = dict()
        self.email = 'briano@bioteam.net'


    def parse(self):
        for file in self.files:
            base = os.path.basename(file).split('.')[0]
            # self.hits[base] = dict()
            matches = re.match(r'(\w+-\w+)_(\w+-\w+)', base)
            for qresult in SearchIO.parse(file, 'hmmer3-tab'):
                if not qresult:
                    if self.verbose:
                        print("No match:\t{0}\t{1}".format(matches[1], matches[2]))
                    continue
                pids = [hit.id for hit in qresult]
                taxarray = self.get_taxid(pids)
                self.hits[base] = self.get_lineage(taxarray)


    def get_taxid(self, pids):
        Entrez.email = self.email
        # Split the list of ids into "batches" of ids for Entrez
        num_chunks = int(len(pids)/self.chunk) + 1
        taxarray = []
        errorarray =[]
        for id_chunk in numpy.array_split(numpy.array(pids), num_chunks):
            # Some protein ids do not have taxonomy ids according to Elink
            try:
                if self.verbose:
                    print("Protein ids: {}".format(id_chunk))
                handle = Entrez.elink(dbfrom="protein", db="taxonomy", id=id_chunk)
                results = Entrez.read(handle)
                handle.close()
                # Protein id, Taxonomy id
                for num, result in enumerate(results):
                    taxarray.append([id_chunk[num], result["LinkSetDb"][0]["Link"][0]["Id"]])
            except:
                if self.verbose:
                    print("Problem getting taxids for: {}".format(id_chunk))
                # Collect protein ids without taxonomy ids
                errorarray = errorarray + list(id_chunk)                
        return taxarray


    def download_hits(self):
        if not self.download:
            return
        for file in self.hits:
            pids = ','.join(self.hits[file].keys())
            try:
                if self.verbose:
                    print("Downloading records: {}".format(pids))
                handle = Entrez.efetch(
                    db="protein",
                    rettype='fasta',
                    retmode="text",
                    id=pids )
                self.fasta[file] = list(SeqIO.parse(handle, 'fasta'))
                handle.close()
            except (RuntimeError) as exception:
                print("Error retrieving sequences using id '" +
                    str(pids) + "':" + str(exception))
        self.write()


    def write(self):
        for file in self.fasta:
            seqfile = file + '-hits.fa'
            if self.verbose:
                print("Writing {}".format(seqfile))
            SeqIO.write(self.fasta[file], seqfile, 'fasta')


    def align_hits_to_hmm(self):
        if not self.align:
            return
        for file in self.hits:
            fastafile = file + '-hits.fa'
            if os.path.exists(fastafile) and os.stat(fastafile).st_size != 0:
                gene = re.match(r'(\w+-\w+)_\w+', file)[1]
                # No guarantee that the HMM is in the current dir so look for it
                hmm = [f for f in glob.glob('**/' + gene + '.hmm', recursive=True)][0]
                if self.verbose:
                    print("Creating alignment with hmmalign: {}".format(file + '-hits.sto'))
                subprocess.run(['hmmalign','--amino', '-o', file + '-hits.sto',
                    hmm, fastafile, ], check=True)


    def filter(self):
        if self.taxfilter:
            for file in self.hits:
                filtered = dict()
                for hit in self.hits[file]:
                    if self.taxfilter not in self.hits[file][hit]:
                        filtered[hit] = self.hits[file][hit]
                self.hits[file] = filtered
        if self.verbose:
            for file in self.hits:
                for hit in self.hits[file]:
                    print("Lineage is: {}".format(self.hits[file][hit]))


    def get_lineage(self, taxarray):
        Entrez.email = self.email
        taxids = [ elem[1] for elem in taxarray ]
        handle = Entrez.efetch(db="taxonomy", id=','.join(taxids))
        results = Entrez.read(handle)
        handle.close()
        # With "LineageEx" you get the NCBI taxon identifiers of the clades
        taxdict = {}
        for num, result in enumerate(results):
            taxdict[taxarray[num][0]] = result["Lineage"]
        return taxdict

if __name__ == "__main__":
    main()

