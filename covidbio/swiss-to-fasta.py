#!/usr/bin/env python3

import sys
import os

from Bio import SearchIO
from Bio import SeqIO
from Bio.Seq import Seq

#from covidbio.utilities import read_strains, read_config

#short script to convert swissprot to fasta

records = SeqIO.parse("THIS_IS_YOUR_INPUT_FILE.swiss", "swiss")
count = SeqIO.write(records, "THIS_IS_YOUR_OUTPUT_FILE.fasta", "fasta")
print("Converted %i records" % count)
