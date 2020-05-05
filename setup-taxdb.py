#!/usr/bin/env python3

import argparse
import sys
from ete2 import NCBITaxa

parser = argparse.ArgumentParser()
parser.add_argument('-verbose', default=False, action='store_true', help="Verbose")


def main():
    builder = Setup_Taxdb(args.verbose)


class Sequp_Taxdb:

    def __init__(self, verbose):
        self.verbose = verbose
        ncbi = NCBITaxa()


if __name__ == "__main__":
    main()
