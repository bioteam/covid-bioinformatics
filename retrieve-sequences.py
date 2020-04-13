#!/usr/bin/env python

"""Usage: python3 retrieve-sequences.py sequences.csv

Downloads all sequences listed in the provided CSV to individual
GenBank files in the current working directory.

See the file `genbank-sequences.csv` in the `seq` directory for an
example format.

"""

import sys
from Bio import Entrez

sequence_file = sys.argv[1]
Entrez.email = 'karl@bioteam.net'

try:
    ids = []
    fp = open(sequence_file)
    line = fp.readline()
    assert line.split(',')[0] == 'Accession'
    for line in fp:
        ids.append(line.split(',')[0])

    # Change db to 'protein' if retrieving protein sequences.
    handle = Entrez.efetch(
        db='nuccore', id=','.join(ids),
        rettype='gb', retmode='text')

    fp = None
    for line in handle.read().splitlines():
        if line.startswith('LOCUS'):
            acc = line.split()[1]
            print(f'Writing {acc}.gb...')
            fp = open(acc + '.gb', 'w')

        if fp:
            fp.write(line + '\n')
            
        if line.startswith('//'):
            fp.close()
            fp = None
            
except IOError as ex:
    print(str(exception))
