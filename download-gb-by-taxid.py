#!/usr/local/bin/python3

import re
import os
import argparse
import numpy
import yaml
from Bio import Entrez
from Bio import SeqIO

# NCBI Taxonomy ids:
# 1042633: single record for testing
# 2697049: Severe acute respiratory syndrome coronavirus 2
# 694009 (parent of 2697049): Severe acute respiratory syndrome-related coronavirus
# https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?lvl=0&id=694009
parser = argparse.ArgumentParser()
parser.add_argument('-t', default=694009, dest='taxid', help="Taxonomy id")
parser.add_argument('-r', default=True, dest='recurse', help="Recursive retrieval of child tax IDs", type=bool)
parser.add_argument('-f', default='gb', dest='format', help="Input and output format")
parser.add_argument('-e', default="briano@bioteam.net", dest='email', help="Email")
parser.add_argument('-min', default=25000, type=int, help="Minimum length")
parser.add_argument('-max', type=int, help="Maximum length")
parser.add_argument('-split', dest='split', default=True, help="Split into separate files")
parser.add_argument('-no-split', dest='split', action='store_false', help="Create one 'taxid' file")
parser.add_argument('-verbose', action='store_true', dest='verbose', help="Verbose")
parser.add_argument('-quiet', default=False, dest='verbose', help="Quiet")
parser.add_argument('-cloud', default=False, type=bool, help="Cloud mode")
args = parser.parse_args()

def main():
    entrez = DownloadGbByTaxid(args.email, args.taxid, args.format, args.min, args.max,
                                     args.split, args.recurse, args.verbose, args.cloud)
    entrez.search()
    entrez.efetch()
    entrez.filter()
    entrez.write()


class DownloadGbByTaxid:

    def __init__(self, email, taxid, format, min_len, max_len, split, recurse, verbose, cloud):
        self.email = email
        self.taxid = taxid
        self.format = format
        self.min_len = min_len
        self.max_len = max_len
        self.split = split
        self.recurse = recurse
        self.verbose = verbose
        self.nt_ids = []
        self.records = []
        self.retmax = 100
        self.cloud = cloud
        if self.cloud:
            with open('config.yaml') as file:
                config = yaml.load(file, Loader=yaml.FullLoader)
            scriptname = os.path.basename(__file__)
            self.email = config[scriptname]['email']
            self.min_len = config[scriptname]['min_len']
            self.max_len = config[scriptname]['max_len']
            self.split = config[scriptname]['split']
            self.recurse = config[scriptname]['recurse']
            self.verbose = config[scriptname]['verbose']
            self.taxid = config[scriptname]['taxid']
            self.format = config[scriptname]['format']
            self.retmax = config[scriptname]['retmax']


    def search(self):
        nummatch = re.match(r'^\d+$', str(self.taxid))
        if not nummatch and self.verbose:
            print("String '" + self.taxid + "' is not an NCBI taxon id")
            return

        Entrez.email = self.email

        if self.recurse == True:
            try:
                handle = Entrez.esearch(db="nuccore", 
                                        idtype="acc",
                                        retmax=5000, 
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
                print("No sequences found with id " + self.taxid)
                return

            for link in links[0]["LinkSetDb"][0]["Link"]:
                self.nt_ids.append(link["Id"])

        if self.verbose:
            print("Esearch id count for Taxonomy id {0}: {1}".format(
                self.taxid, len(self.nt_ids)))

    def efetch(self):
        Entrez.email = self.email
        # Split the list of ids into "batches" of ids for Entrez
        num_chunks = int(len(self.nt_ids)/self.retmax) + 1

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
                seqfile = record.name + '.' + self.format
                if self.verbose:
                    print("Writing {}".format(seqfile))
                SeqIO.write(record, seqfile, self.format)
        else:
            seqfile = 'taxid-' + str(self.taxid) + '.' + self.format
            if self.verbose:
                print("Writing {}".format(seqfile))
            SeqIO.write(self.records, seqfile, self.format)


if __name__ == "__main__":
    main()
