#!/usr/bin/env python3

import argparse
import sys
import os
import re
import subprocess
import numpy
from Bio import SearchIO
from Bio import Entrez
from Bio import SeqIO
from covidbio.utilities import read_config

config = read_config()

parser = argparse.ArgumentParser()
parser.add_argument('-verbose', default=False, action='store_true', help="Verbose")
parser.add_argument('-download', default=False, action='store_true', help="Download hit sequences")
parser.add_argument('-align', default=False, action='store_true', help="Align hits to HMM")
parser.add_argument('-taxfilter', default=None, help="Exclude clade using NCBI clade name")
parser.add_argument('-lexfilter', default=None, help="Exclude clade using search string")
parser.add_argument('-chunk', default=10, help="Number of ids to send to Elink")
parser.add_argument('-data_dir', default=config['DATA_DIR'], help="Location for all strain-specific directories")
parser.add_argument('-email', default=config['EMAIL'], help="Email for Entrez")
parser.add_argument('-api_key', help="Entrez API key")
parser.add_argument('-strain', default=config['STRAIN'], help="Strain name")
parser.add_argument('files', nargs='+', help='File names')
args = parser.parse_args()

'''

'''

def main():
    query = Parse_Hmmsearch(args.verbose, args.download, args.align, args.taxfilter, 
        args.lexfilter, args.chunk, args.data_dir, args.email, args.api_key, args.strain, args.files)
    for f in query.files:
        pids, fname = query.parse(f)
        if pids == None or pids == []:
            continue
        taxarray = query.get_taxid(pids)
        lineages = query.get_lineage(taxarray)
        filtered_pids = query.filter(lineages)
        query.download_hits(filtered_pids, fname)
        query.align_hits_to_hmm(fname)

class Parse_Hmmsearch:

    def __init__(self, verbose, download, align, taxfilter, lexfilter, chunk, data_dir, email, api_key, strain, files):
        self.email = email
        self.strain = strain
        self.data_dir = data_dir
        self.verbose = verbose
        self.download = download
        self.align = align
        self.taxfilter = taxfilter
        self.lexfilter = lexfilter
        self.chunk = chunk
        self.api_key = api_key
        self.files = files
        self.cov_dir = os.path.join(self.data_dir, self.strain)
        if not os.path.isdir(self.cov_dir):
            sys.exit("Directory {} does not exist".format(self.cov_dir))
        if not self.api_key and 'NCBI_API_KEY' in os.environ.keys():
            self.api_key = os.environ['NCBI_API_KEY']

    def parse(self, file):
        '''
        Parse hmmsearch output and filter by arbitrary search string 
        if one is specified by -lexfilter.
        '''
        fname = os.path.basename(file).split('.')[0]
        if os.stat(file).st_size == 0:
            return None, fname
        matches = re.match(r'(\w+-\w+)_([^.]+)', fname)
        try:
            qresult = SearchIO.read(file, 'hmmer3-tab')
        except:
            if self.verbose:
                print("No match:\t{0}\t{1}".format(matches[1], matches[2]))
            return None, fname
        if self.lexfilter:
            patt = re.compile(self.lexfilter, re.IGNORECASE)
            return [ hit.id for hit in qresult if not re.search(patt, hit.description)], fname
        else:
            return [ hit.id for hit in qresult], fname

    def get_taxid(self, pids):
        '''
        Returns an array of tuples where tup[0] is a protein id 
        and tup[1] is the Taxonomy id, for example:
        [('P0DTC2.1', '2697049'), ('P59594.1', '694009')]
        '''
        Entrez.email = self.email
        Entrez.api_key = self.api_key
        # Split the list of ids into "batches" of ids for Entrez
        num_chunks = len(pids)/int(self.chunk) + 1
        taxarray = list()
        errorarray = list()
        for id_chunk in numpy.array_split(numpy.array(pids), num_chunks):
            # Some protein ids do not have taxonomy ids according to Elink
            try:
                if self.verbose:
                    print("Protein ids: {}".format(id_chunk))
                handle = Entrez.elink(dbfrom="protein", db="taxonomy", id=id_chunk)
                results = Entrez.read(handle)
                handle.close()
                # Protein id, Taxonomy id
                for num, result in enumerate(results):
                    taxarray.append((id_chunk[num], result["LinkSetDb"][0]["Link"][0]["Id"]))
            except:
                if self.verbose:
                    print("Problem getting taxids for: {}".format(id_chunk))
                # Collect protein ids without taxonomy ids
                errorarray = errorarray + list(id_chunk)                
        return taxarray

    def download_hits(self, pids, fname):
        if not self.download:
            return
        pids = ','.join(pids)
        Entrez.email = self.email
        Entrez.api_key = self.api_key
        try:
            if self.verbose:
                print("Downloading records: {}".format(pids))
            handle = Entrez.efetch(
                    db="protein",
                    rettype='fasta',
                    retmode="text",
                    id=pids )
            records = list(SeqIO.parse(handle, 'fasta'))
            handle.close()
        except (RuntimeError) as exception:
            print("Error retrieving sequences using id '" +
                str(pids) + "':" + str(exception))
        self.write(records, fname)

    def align_hits_to_hmm(self, fname):
        if not self.align:
            return
        fastafile = fname + '-hits-no-' + self.taxfilter + '.fa' if self.taxfilter else fname + '-hits.fa'
        fastafile = os.path.join(self.cov_dir, fastafile)
        if os.path.exists(fastafile) and os.stat(fastafile).st_size != 0:
            gene = re.match(r'(\w+-\w+)_\w+', fname)[1]
            hmm = os.path.join(self.cov_dir, gene + '.hmm')
            outfile = fname + '-hits-no-' + self.taxfilter + '.sto' if self.taxfilter else fname + '-hits.sto'
            outfile = os.path.join(self.cov_dir, outfile)
            if self.verbose:
                print("Creating alignment with hmmalign: {}".format(outfile))
            subprocess.run(['hmmalign','--amino', '-o', outfile,
                            hmm, fastafile, ], check=True)

    def filter(self, lineages):
        if self.taxfilter:
            filtered = [ pid for pid in lineages.keys() if self.taxfilter not in lineages[pid] ]
            if self.verbose:
                print("Filtered lineages: {}".format(filtered))
            return filtered
        else:
            return lineages

    def get_lineage(self, taxarray):
        '''
        Returns a dict where the key is a Protein id and the value is a lineage,
        for example:
        {'P0DTC2.1': 'Viruses; Riboviria; Orthornavirae; Pisuviricota; Pisoniviricetes; \
        Nidovirales; Cornidovirineae; Coronaviridae; Orthocoronavirinae; Betacoronavirus; \
        Sarbecovirus; Severe acute respiratory syndrome-related coronavirus'}
        '''
        Entrez.email = self.email
        Entrez.api_key = self.api_key
        taxids = [ elem[1] for elem in taxarray ]
        if self.verbose:
            print("Fetching from 'taxonomy': {}".format(taxids))
        handle = Entrez.efetch(db="taxonomy", id=','.join(taxids))
        results = Entrez.read(handle)
        handle.close()
        # With "LineageEx" you get the NCBI taxon identifiers of the clades
        taxdict = {}
        for num, result in enumerate(results):
            taxdict[taxarray[num][0]] = result["Lineage"]
        if self.verbose:
            print("Retrieved from 'taxonomy': {}".format(taxdict))
        return taxdict

    def write(self, records, fname):
        seqfile = fname + '-hits-no-' + self.taxfilter + '.fa' if self.taxfilter else fname + '-hits.fa'
        seqfile = os.path.join(self.cov_dir, seqfile)
        if self.verbose:
            print("Writing {}".format(seqfile))
        SeqIO.write(records, seqfile, 'fasta')

if __name__ == "__main__":
    main()

