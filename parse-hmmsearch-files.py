#!/usr/bin/env python3

import argparse
from pathlib import Path
import os
from taxadb.taxid import TaxID
from taxadb.accessionid import AccessionID
from Bio import SearchIO

parser = argparse.ArgumentParser()
parser.add_argument('-verbose', default=False, action='store_true', help="Verbose")
parser.add_argument('-db', default='    ', help="Database")
parser.add_argument('files', nargs='+', help='File names')
args = parser.parse_args()


def main():
    query = Query_Taxadb(args.verbose, args.files)
    query.read()

class Query_Taxadb:

    def __init__(self, verbose, files):
        self.verbose = verbose
        self.files = files
        self.dbpath = os.path.join(Path.home(),'taxonomy','taxadb.sqlite')
        self.taxdb = TaxID(dbtype='sqlite', dbname=self.dbpath)
        self.accession = AccessionID(dbtype='sqlite', dbname=self.dbpath)
 
 
    def read:
        for file in self.files:
            matches = re.match(r'(\w+-\w+)_(\w+-\w+)', os.path.basename(file))
            # Turn a Generator into a list
            qresult = list(SearchIO.parse(file, 'hmmer3-tab'))
            if not qresult:
                print("No match:\t{0}\t{1}".format(matches[1], matches[2]))
                continue
            hit = qresult[0].hits[0]
            taxids = accession.taxid(hit.id)
            lineage = taxid.lineage_name(taxid, reverse=True)

            #     if qresult[0].id.split('-')[0] != hit.id.split('-')[0]:
            #         print("{0}\t{1}\t{2}".format(qresult[0].id, hit.id, hit.evalue))


if __name__ == "__main__":
    main()

