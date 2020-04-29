#!/usr/local/bin/python3

import argparse
import os.path
import subprocess
import sys
from Bio import SearchIO


parser = argparse.ArgumentParser()
parser.add_argument('-verbose', default=False, action='store_true', help="Verbose")
# parser.add_argument('-aligner', default='clustalo', help="Alignment application")
parser.add_argument('files', nargs='+', help='File names')
args = parser.parse_args()


def main():
    builder = Parse_Hmms_Results(args.verbose, args.files)
    builder.read()

class Parse_Hmms_Results:

    def __init__(self, verbose, files):
        self.verbose = verbose
        self.files = files
        # self.seqs, self.alns, self.hmms = dict(), dict(), dict()

    def read(self):
        full_paths = [os.path.join(os.getcwd(), path) for path in self.files]
        for path in full_paths:
            for qresult in SearchIO.parse(path, 'hmmer3-tab'):
                print("{0}\t{1}\t{2}".format(qresult.id, qresult.accession, len(qresult.hits)))
                hit = qresult.hits[0]
                print("{0}\t{1}\t{2}".format(hit.id, hit.evalue, hit.description))


if __name__ == "__main__":
    main()


'''
QueryResult	
accession	query accession (if present)
description	query sequence description
id	query name

Hit
accession	hit accession
bias	hit-level bias
bitscore	hit-level score
description	hit sequence description
cluster_num	clu column
domain_exp_num	exp column
domain_included_num	inc column
domain_obs_num	dom column
domain_reported_num	rep column
env_num	env column
evalue	hit-level evalue
id	target name
overlap_num	ov column
region_num	reg column

HSP	
bias	bias of the best domain
bitscore	bitscore of the best domain
evalue	evalue of the best domain
'''