#!/usr/bin/env python3

import argparse
from Bio import SearchIO
from Bio import Entrez

parser = argparse.ArgumentParser()
parser.add_argument('-verbose', default=False, action='store_true', help="Verbose")
parser.add_argument('-download', default=False, action='store_true', help="Download hit sequences")
parser.add_argument('-align', default=False, action='store_true', help="Align hits to HMM")
parser.add_argument('-taxfilter', default=None, help="Exclude clade")
parser.add_argument('files', nargs='+', help='File names')
args = parser.parse_args()


def main():
    query = Parse_Hmmsearch(args.verbose, args.download, args.align, args.taxfilter, args.files)
    query.parse()
    query.filter()
    query.download_hits()
    query.align_hits_to_hmm()

class Parse_Hmmsearch:

    def __init__(self, verbose, download, align, taxfilter, files):
        self.verbose = verbose
        self.download = download
        self.align = align
        self.taxfilter = taxfilter
        self.files = files
        self.hits = dict()
        self.email = 'briano@bioteam.net'
 
    def parse(self):
        for file in self.files:
            matches = re.match(r'(\w+-\w+)_(\w+-\w+)', os.path.basename(file))
            for qresult in SearchIO.parse(file, 'hmmer3-tab'):
                if not qresult:
                    if self.verbose:
                        print("No match:\t{0}\t{1}".format(matches[1], matches[2]))
                    continue
                for hit in qresult:
                    taxid = self.get_taxid(hit.id)
                    # if self.verbose:
                    #     print("Hit taxonomy id: {}".format(taxid))
                    self.hits[hit.id] = self.get_lineage(taxid)


    def get_taxid(self, pid):
        Entrez.email = self.email
        handle = Entrez.elink(dbfrom="protein", db="taxonomy", id=pid)
        result = Entrez.read(handle)
        handle.close()
        taxid = result[0]["LinkSetDb"][0]["Link"][0]["Id"]
        return taxid

    def download_hits(self):
        if not self.download:
            return
        try:
            if self.verbose:
                print("Downloading records: {}".format(id_chunk))
            handle = Entrez.efetch(
                    db="protein",
                    rettype='fasta',
                    retmode="text",
                    id=','.join(self.hits.keys())
            )
            self.fasta = handle.read()
            handle.close()
        except (RuntimeError) as exception:
            print("Error retrieving sequences using id '" +
                  str(pids) + "':" + str(exception))

    def align_hits_to_hmm(self):
        if not self.align:
            return

    def filter(self):
        if self.taxfilter:
            filtered = dict()
            for hit in self.hits:
                if self.taxfilter not in self.hits[hit]:
                    filtered[hit] = self.hits[hit]
            self.hits = filtered
        if self.verbose:
            for hit in self.hits:
                print("Lineage is: {}".format(self.hits[hit]))

    def get_lineage(self, taxid):
        Entrez.email = self.email
        handle = Entrez.efetch(db="taxonomy", id=taxid)
        records = Entrez.read(handle)
        # With "LineageEx" you get the NCBI taxon identifiers of the clades
        return records[0]["Lineage"]
 

if __name__ == "__main__":
    main()

