#!/usr/bin/env python3

import argparse
import sys
import os
import re
import numpy
import yaml
from Bio import Entrez
from Bio import SeqIO
from covidbio.utilities import read_strains, read_config

'''
NCBI Taxonomy ids:
333387: single record for testing
2697049: Severe acute respiratory syndrome coronavirus 2
694009 (parent of 2697049): Severe acute respiratory syndrome-related coronavirus
https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?lvl=0&id=694009
'''

config = read_config()

parser = argparse.ArgumentParser()
parser.add_argument('-r', default=True, dest='recurse', help="Recursive retrieval of child tax IDs", type=bool)
parser.add_argument('-f', default='gb', dest='format', help="Input and output format")
parser.add_argument('-email', default=config['EMAIL'], help="Email for Entrez")
parser.add_argument('-min_len', default=25000, type=int, help="Minimum length")
parser.add_argument('-max_len', type=int, help="Maximum length")
parser.add_argument('-verbose', action='store_true', help="Verbose")
parser.add_argument('-retmax', default=50000, type=int, help="Entrez retmax")
parser.add_argument('-chunk', default=500, type=int, help="eFetch batch size")
parser.add_argument('-api_key', help="Entrez API key")
parser.add_argument('-json', action='store_true', help="Create JSON for Gen3")
parser.add_argument('-no-fetch', action='store_false', dest='fetch', help="Do not download")
parser.add_argument('-strain', default=config['STRAIN'], help="Strain name")
parser.add_argument('-date_filter', action='store_true', help="Filter by years in cov_strains.yaml")
parser.add_argument('-data_dir', default=config['DATA_DIR'], help="Location for all strain-specific directories")
args = parser.parse_args()

def main():
    entrez = DownloadGbByTaxid(args.email, args.format, args.min_len, args.max_len,
                                args.recurse, args.verbose, args.retmax,
                                args.chunk, args.api_key, args.json, args.fetch, args.strain,
                                args.date_filter, args.data_dir)
    entrez.search()
    entrez.efetch()
    entrez.filter()
    entrez.write()


class DownloadGbByTaxid:

    def __init__(self, email, format, min_len, max_len, recurse, verbose, retmax,
                 chunk, api_key, json, fetch, strain, date_filter, data_dir):
        self.email = email
        self.strain = strain
        self.data_dir = data_dir
        self.format = format
        self.min_len = min_len
        self.max_len = max_len
        self.recurse = recurse
        self.verbose = verbose
        self.retmax = retmax
        self.chunk = chunk
        self.api_key = api_key
        self.fetch = fetch
        self.json = json
        self.date_filter = date_filter
        self.nt_ids = list()
        self.records = list()
        if not self.api_key and 'NCBI_API_KEY' in os.environ.keys():
            self.api_key = os.environ['NCBI_API_KEY']
        self.cov_dir = os.path.join(self.data_dir, self.strain)
        if not os.path.isdir(self.cov_dir):
            sys.exit("Directory {} does not exist".format(self.cov_dir))
        # Get details about the specific strain (e.g. MERS, COV2)
        strains = read_strains()
        self.taxid = strains[self.strain]['taxid']
        if self.date_filter:
            self.date_filter = strains[self.strain]['years']

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
        filtered = list()
        for record in self.records:
            # Minimum length requirement
            if self.min_len:
                if len(record) < self.min_len:
                    continue
            # Filter by year, from cov_strains.yaml
            if self.date_filter:
                pub_year = int(record.annotations['date'].split('-')[2])
                if pub_year not in self.date_filter:
                    continue
            filtered.append(record)
        self.records = filtered

    def write(self):
        for record in self.records:
            seqfile = os.path.join(
                self.cov_dir, record.name + '.' + self.format)
            if self.verbose:
                print("Writing {}".format(seqfile))
            SeqIO.write(record, seqfile, self.format)
            if self.json:
                from utilities import make_genome_json
                genome_json = make_genome_json(record.name)
                with open(record.name + '.json', 'w') as out:
                    out.write(genome_json)


if __name__ == "__main__":
    main()
