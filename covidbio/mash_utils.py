#!/usr/bin/env python3

import sys
import subprocess

def run_mash(fname):
    cmd = ['mash', 'dist', '-i', '-a', fname, fname]
    try:
        proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except (subprocess.CalledProcessError) as exception:
        print("Error: {}".format(exception))
        sys.exit("Error running 'mash dist' on {}".format(fname))
    # Return list of lists
    return [e.split('\t') for e in proc.stdout.split('\n') if e]

