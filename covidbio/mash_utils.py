#!/usr/bin/env python3

import sys
import subprocess
import os
from Bio import SeqIO

'''
Use the mash application to do pairwise sequence comparision using 
kmers, and use similarity data to identify closely related sequences 
that are removed to create an "nr" file.
'''
def make_nr(fname, informat='fasta', threshold=0.1):
    results = run_mash(fname)
    # Arbitrarily pick first of sorted pair of ids to skip and ignore sequence compared to itself
    skip = set([sorted([e[0], e[1]])[0] for e in results if float(e[2]) < threshold and e[0] != e[1]])
    # For example, input is viruses.dat, output is viruses-nr.fa
    with open(os.path.basename(fname).split('.')[0] + '-nr.fa', 'w') as out:
        for seq in SeqIO.parse(fname, informat):
            if seq.id not in skip:
                SeqIO.write(seq, out, 'fasta')

'''
mash results example:

[['CYS1_DICDI', 'CYS1_DICDI', '0', '0', '335/335'], 
['CATH_HUMAN', 'CYS1_DICDI', '0.567663', '2.38275e-14', '2/660'], 
['CATL_HUMAN', 'CATH_HUMAN', '1', '1', '0/652'],
['CATH_HUMAN', 'CATH_RAT', '0.120713', '0', '110/542']]
['CATH_RAT', 'CATH_HUMAN', '0.120713', '0', '110/542']]
'''
def run_mash(fname):
    # The same file acts as both query and reference
    cmd = ['mash', 'dist', '-i', '-a', fname, fname]
    try:
        proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except (subprocess.CalledProcessError) as exception:
        print("Error: {}".format(exception))
        sys.exit("Error running 'mash dist' on {}".format(fname))
    # Return mash results as list of lists
    return [e.split('\t') for e in proc.stdout.split('\n') if e]

