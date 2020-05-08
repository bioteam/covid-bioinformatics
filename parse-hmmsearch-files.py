#!/usr/bin/env python3

import argparse
from pathlib import Path
import os
import re
from taxadb.taxid import TaxID
from taxadb.accessionid import AccessionID
from Bio import SearchIO

parser = argparse.ArgumentParser()
parser.add_argument('-verbose', default=False, action='store_true', help="Verbose")
parser.add_argument('-db', default=os.path.join(Path.home(),'taxonomy','taxadb.sqlite'), 
    help="Path to database")
parser.add_argument('-filter', default="Viridae", help="Exclude")
parser.add_argument('files', nargs='+', help='File names')
args = parser.parse_args()


def main():
    query = Query_Taxadb(args.verbose, args.db, args.filter, args.files)
    query.read()

class Query_Taxadb:

    def __init__(self, verbose, db, filter, files):
        self.verbose = verbose
        self.db = db
        self.filter = filter
        self.files = files
        self.taxid = TaxID(dbtype='sqlite', dbname=self.db)
        self.accession = AccessionID(dbtype='sqlite', dbname=self.db)
 
 
    def read(self):
        for file in self.files:
            matches = re.match(r'(\w+-\w+)_(\w+-\w+)', os.path.basename(file))
            for qresult in SearchIO.parse(file, 'hmmer3-tab'):
                if not qresult:
                    if self.verbose:
                        print("No match:\t{0}\t{1}".format(matches[1], matches[2]))
                    continue
                hit = qresult[0].hits[0]
                acc = hit.id.split('.')[0]
                lineages = self.query(acc)
                if self.verbose:
                    print("Lineages: {}".format(lineages))

    def query(self, acc):
        lineages = []
        for taxid in self.accession.taxid([acc]):
            lineage = self.taxid.lineage_name(taxid, reverse=True)
            lineages.append(lineage)
        return lineages


if __name__ == "__main__":
    main()

