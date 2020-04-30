#!/usr/local/bin/python3

import argparse
import os.path
import sys
import re
from Bio import SearchIO


parser = argparse.ArgumentParser()
parser.add_argument('-verbose', default=False,
                    action='store_true', help="Verbose")
parser.add_argument('files', nargs='+', help='File names')
args = parser.parse_args()


def main():
    builder = Parse_Hmms_Results(args.verbose, args.files)
    builder.read()


class Parse_Hmms_Results:

    def __init__(self, verbose, files):
        self.verbose = verbose
        self.files = files

    def read(self):
        # full_paths = [os.path.join(os.getcwd(), path) for path in self.files]
        for file in self.files:
            matches = re.match(r'(\w+-\w+)_(\w+-\w+)', file)
            # Turn a Generator into a list
            qresult = list(SearchIO.parse(file, 'hmmer3-tab'))
            if not qresult:
                print("No match:\t{0}\t{1}".format(matches[1], matches[2]))
            else:
                print("{0}\t{1}\t{2}".format(qresult[0].id, qresult[0].accession, len(qresult[0].hits)))
                hit = qresult[0].hits[0]
                print("{0}\t{1}".format(hit.id, hit.evalue))


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
