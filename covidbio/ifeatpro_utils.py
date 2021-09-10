#!/usr/bin/env python3

import os
import tempfile
from ifeatpro.features import get_feature

def run_ifeatpro(feature, seqstr, pid):
    tmpfasta = tempfile.NamedTemporaryFile()
    with open(tmpfasta.name, 'w') as f:
        f.write(">" + pid + "\n" + seqstr)
    tmpdir = tempfile.TemporaryDirectory()
    get_feature(tmpfasta.name, feature, tmpdir.name)
    with open(os.path.join(tmpdir.name, feature + ".csv")) as f:
        l = f.readline()
    return l.split(',')[1:]

