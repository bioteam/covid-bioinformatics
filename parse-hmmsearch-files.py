#!/usr/bin/env python3

import argparse
from pathlib import Path
import os
import re
from Bio import SearchIO
from Bio import Entrez

parser = argparse.ArgumentParser()
parser.add_argument('-verbose', default=False, action='store_true', help="Verbose")
parser.add_argument('-taxfilter', help="Exclude")
parser.add_argument('files', nargs='+', help='File names')
args = parser.parse_args()


def main():
    query = Parse_Hmmsearch(args.verbose, args.taxfilter, args.files)
    query.read()

class Parse_Hmmsearch:

    def __init__(self, verbose, taxfilter, files):
        self.verbose = verbose
        self.taxfilter = taxfilter
        self.files = files
        self.email = 'briano@bioteam.net'
 
    def read(self):
        for file in self.files:
            matches = re.match(r'(\w+-\w+)_(\w+-\w+)', os.path.basename(file))
            for qresult in SearchIO.parse(file, 'hmmer3-tab'):
                if not qresult:
                    if self.verbose:
                        print("No match:\t{0}\t{1}".format(matches[1], matches[2]))
                    continue
                for hit in qresult:
                    taxid = self.get_taxid(hit.id)
                    if self.verbose:
                        print("Hit taxonomy id: {}".format(taxid))
                    lineage = self.get_lineage(taxid)
                    if self.verbose and lineage:
                        print("Hit lineage: {}".format(lineage))


    def get_taxid(self, pid):
        Entrez.email = self.email
        handle = Entrez.elink(dbfrom="protein", db="taxonomy", id=pid)
        result = Entrez.read(handle)
        handle.close()
        taxid = result[0]["LinkSetDb"][0]["Link"][0]["Id"]
        return taxid


    def get_lineage(self, taxid):
        Entrez.email = self.email
        handle = Entrez.efetch(db="taxonomy", id=taxid)
        records = Entrez.read(handle)
        # With "LineageEx" you’ll get the NCBI taxon identifiers of the lineage entries
        if self.taxfilter:
            if self.taxfilter in records[0]["Lineage"]:
                return None
        return records[0]["Lineage"]
 

if __name__ == "__main__":
    main()

