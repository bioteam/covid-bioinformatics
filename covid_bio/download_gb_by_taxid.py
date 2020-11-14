#!/usr/bin/env python3

import argparse
import sys
import os
import re
import numpy
import yaml
from Bio import Entrez
from Bio import SeqIO
from vars import COV_DIR, EMAIL

# NCBI Taxonomy ids:
# 333387: single record for testing
# 2697049: Severe acute respiratory syndrome coronavirus 2
# 694009 (parent of 2697049): Severe acute respiratory syndrome-related coronavirus
# https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?lvl=0&id=694009
parser = argparse.ArgumentParser()
parser.add_argument('-t, -strain', default='COV2', dest='strain', help="Taxonomy id")
parser.add_argument('-r', default=True, dest='recurse', help="Recursive retrieval of child tax IDs", type=bool)
parser.add_argument('-f', default='gb', dest='format', help="Input and output format")
parser.add_argument('-e', default=EMAIL, dest='email', help="Email for Entrez")
parser.add_argument('-min_len', default=25000, type=int, help="Minimum length")
parser.add_argument('-max_len', type=int, help="Maximum length")
parser.add_argument('-split', default=True, help="Split into separate files")
parser.add_argument('-no-split', dest='split', action='store_false', help="Create one 'taxid' file")
parser.add_argument('-verbose', action='store_true', help="Verbose")
parser.add_argument('-retmax', default=10000, type=int, help="Entrez retmax")
parser.add_argument('-chunk', default=500, type=int, help="eFetch batch size")
parser.add_argument('-api_key', help="Entrez API key")
parser.add_argument('-json', action='store_true', help="Create JSON for Gen3")
parser.add_argument('-no-fetch', action='store_false', dest='fetch', help="Do not download")
parser.add_argument('-cov_dir', default=COV_DIR, help="Destination directory")
args = parser.parse_args()

def main():
    entrez = DownloadGbByTaxid(args.email, args.strain, args.format, args.min_len, args.max_len,
                                args.split, args.recurse, args.verbose, args.retmax,
                                args.chunk, args.api_key, args.json, args.fetch, args.cov_dir)
    entrez.search()
    entrez.efetch()
    entrez.filter()
    entrez.write()


class DownloadGbByTaxid:

    def __init__(self, email, strain, format, min_len, max_len, split, recurse, verbose, retmax,
                 chunk, api_key, json, fetch, cov_dir):
        self.email = email
        self.strain = strain
        self.format = format
        self.min_len = min_len
        self.max_len = max_len
        self.split = split
        self.recurse = recurse
        self.verbose = verbose
        self.retmax = retmax
        self.chunk = chunk
        self.api_key = api_key
        self.fetch = fetch
        self.json = json
        self.cov_dir = cov_dir
        self.nt_ids = []
        self.records = []
        if not self.api_key and 'NCBI_API_KEY' in os.environ.keys():
            self.api_key = os.environ['NCBI_API_KEY']
        self.read_strains()

    def read_strains(self):
        y = os.path.dirname(os.path.abspath(__file__)) + '/cov_strains.yaml'
        with open(y) as file:
            synonyms = yaml.load(file, Loader=yaml.FullLoader)
        self.taxid = synonyms['taxid']

    def search(self):
        nummatch = re.match(r'^\d+$', str(self.taxid))
        if not nummatch and self.verbose:
            print("String '" + self.taxid + "' is not an NCBI taxon id")
            return
        Entrez.email = self.email
        Entrez.api_key = self.api_key
        if self.recurse:
            try:
                handle = Entrez.esearch(db="nuccore", 
                                        idtype="acc",
                                        retmax=self.retmax, 
                                        term="txid{}[Organism:exp]".format(self.taxid))
                records = Entrez.read(handle)
                handle.close()
            except (RuntimeError) as exception:
                print("Error retrieving sequence ids using Taxonomy id '" +
                      str(self.taxid) + "'" + str(exception))

            for link in records['IdList']:
                self.nt_ids.append(link)
        else:
            try:
                links = Entrez.read(
                    Entrez.elink(dbfrom="taxonomy",
                                 db="nucleotide",
                                 idtype="acc",
                                 id=self.taxid))
            except (RuntimeError) as exception:
                print("Error retrieving sequence ids using Taxonomy id '" +
                      str(self.taxid) + "'" + str(exception))

            if len(links[0]["LinkSetDb"]) == 0:
                if self.verbose:
                    print("No sequences found with id " + self.taxid)
                return

            for link in links[0]["LinkSetDb"][0]["Link"]:
                self.nt_ids.append(link["Id"])

        if self.verbose:
            print("Esearch id count for Taxonomy id {0}: {1}".format(
                self.taxid, len(self.nt_ids)))

    def efetch(self):
        if not self.fetch:
            sys.exit(0)
        Entrez.email = self.email
        Entrez.api_key = self.api_key
        # Split the list of ids into "batches" of ids for Entrez
        num_chunks = int(len(self.nt_ids)/self.chunk) + 1

        for id_chunk in numpy.array_split(numpy.array(self.nt_ids), num_chunks):
            try:
                if self.verbose:
                    print("Downloading records: {}".format(id_chunk))
                handle = Entrez.efetch(
                    db="nucleotide",
                    rettype=self.format,
                    retmode="text",
                    id=','.join(id_chunk)
                )
                records = SeqIO.parse(handle, self.format)
                self.records = self.records + list(records)
                handle.close()
            except (RuntimeError) as exception:
                print("Error retrieving sequences using id '" +
                  str(self.taxid) + "':" + str(exception))

    def filter(self):
        if self.min_len:
            filtered = []
            for record in self.records:
                if self.verbose:
                    print("Filtering {}".format(record.id))
                if len(record) >= self.min_len:
                    filtered.append(record)
            self.records = filtered

    def write(self):
        if self.split:
            for record in self.records:
                seqfile = os.path.join(self.cov_dir, record.name + '.' + self.format)
                if self.verbose:
                    print("Writing {}".format(seqfile))
                SeqIO.write(record, seqfile, self.format)
                if self.json:
                    from make_json import make_genome_json
                    genome_json = make_genome_json(record.name)
                    with open(record.name + '.json', 'w') as out:
                        out.write(genome_json)
        else:
            seqfile = os.path.join(self.cov_dir, 'taxid-' + str(self.taxid) + '.' + self.format)
            if self.verbose:
                print("Writing {}".format(seqfile))
            SeqIO.write(self.records, seqfile, self.format)


if __name__ == "__main__":
    main()
