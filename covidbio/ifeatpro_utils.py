#!/usr/bin/env python3

import os
import tempfile
from ifeatpro.features import get_feature

'''
Returns a single array or "feature vector" using ifeatpro given:

- Feature name, e.g. 'ctriad'
- Protein sequence as string, e.g. 'MAVKLSSTTRD'
- Protein id, e.g. 'ARAB_ABC5'
'''
def run_ifeatpro(feature, seqstr, pid):
    tmpfasta = tempfile.NamedTemporaryFile()
    # Write fasta file
    with open(tmpfasta.name, 'w') as f:
        f.write(">" + pid + "\n" + seqstr)
    tmpdir = tempfile.TemporaryDirectory()
    get_feature(tmpfasta.name, feature, tmpdir.name)
    with open(os.path.join(tmpdir.name, feature + ".csv")) as f:
        l = f.readline()
    return l.split(',')[1:]

