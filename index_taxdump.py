#!/usr/bin/env python3

import argparse
import os
import sys
import codecs
from whoosh import index
from whoosh.fields import Schema, TEXT, NUMERIC, STORED
from whoosh.qparser import QueryParser
from whoosh.index import open_dir
from whoosh.query import *
from whoosh.filedb.filestore import FileStorage


'''
>head prot.accession2taxid
accession	accession.version	taxid	gi
A0A009IHW8	A0A009IHW8.1	470	1835922267
A0A023GPI8	A0A023GPI8.1	232300	1027923628
A0A023GPJ0	A0A023GPJ0.2	716541	765680613
A0A023GS28	A0A023GS28.1	37565	1679377489
A0A023GS29	A0A023GS29.1	37565	1679377490
A0A023IWD9	A0A023IWD9.2	262245	1384595658
A0A023IWE0	A0A023IWE0.1	580329	1373488000
A0A023IWE1	A0A023IWE1.1	67723	1375505955
A0A023IWE2	A0A023IWE2.1	1324310	1372175666
'''

parser = argparse.ArgumentParser()
parser.add_argument('-verbose', default=False, action='store_true', help="Verbose")
parser.add_argument('-create', default=False, action='store_true', help="Create index")
parser.add_argument('-query', help="Query term")
parser.add_argument('-index', default="accession_idx", help="Path to index directory")
parser.add_argument('-file', help='Path to taxdump file')
args = parser.parse_args()

def main():
    builder = Index_Taxdump(args.verbose, args.create, args.query, args.index, args.file)
    if builder.create:
        builder.create_schema()
    if builder.query:
        builder.do_query()

class Index_Taxdump:
    def __init__(self, verbose, create, query, index, file):
        self.verbose = verbose
        self.create = create
        self.query = query
        self.index = index
        self.file = file

    def create_schema(self):
        # Set 'stored' to True if you want that value returned in search results
        schema = Schema(accession=TEXT, accession_version=TEXT, 
            taxid=NUMERIC(stored=True), gi=NUMERIC)
        
        if not os.path.exists(self.index):
            os.mkdir(self.index)

        if not self.file:
            sys.exit("prot.accession2taxid file required")

        ix = index.create_in(self.index, schema)
        
        with codecs.open(self.file, 'r', encoding='utf-8') as f:
            # Skip header line 1
            first_line = f.readline()
            with ix.writer() as writer:
                for index,line in enumerate(f):
                    row = line.strip().split('\t')
                    writer.add_document(accession=row[0],
                        accession_version=row[1],
                        taxid=row[2],
                        gi=row[3])
                    if self.verbose and index % 1000 == 0:
                        print("Line {}".format(index))

    def do_query(self):
        storage = FileStorage(self.index)
        ix = storage.open_index()
        qry = QueryParser('accession', schema=ix.schema).parse(self.query)

        with ix.searcher() as searcher:
            results = searcher.search(qry,limit=1)
            if self.verbose:
                print("{}".format(results[0]))

if __name__ == "__main__":
    main()
