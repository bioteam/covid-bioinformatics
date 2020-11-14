#!/usr/bin/env python3

import argparse
import sys
import os
import re
import yaml
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from vars import COV_DIR
from utils import read_synonyms
from utils import read_variants


parser = argparse.ArgumentParser()
parser.add_argument('-f', default='fasta', dest='format', help="Output format")
parser.add_argument('-split', default=False, dest='split', help="Split into separate files")
parser.add_argument('-no-split', action='store_false', dest='split', help="Create one file")
parser.add_argument('-analyze', default=False, action='store_true', help="Parse files without file creation")
parser.add_argument('-verbose', default=False, action='store_true', help="Verbose")
parser.add_argument('-json', default=False, action='store_true', help="Create JSON for Gen3")
parser.add_argument('-cov_dir', default=COV_DIR, help="Destination directory")
parser.add_argument('-host_filter', help="Host name to filter")
parser.add_argument('files', nargs='+', help='File names')
args = parser.parse_args()

'''
Example SeqRecord:

ID: MT123293.2
Name: MT123293
Description: Severe acute respiratory syndrome coronavirus 2 isolate SARS-CoV-2/IQTC03/human/2020/CHN, complete genome
Number of features: 23
/molecule_type=RNA
/topology=linear
/data_file_division=VRL
/date=13-MAR-2020
/accessions=['MT123293']
/sequence_version=2
/keywords=['']
/source=Severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2)
/host=Homo sapiens
/organism=Severe acute respiratory syndrome coronavirus 2
/taxonomy=['Viruses', 'Riboviria', 'Nidovirales', 'Cornidovirineae', 'Coronaviridae', 'Orthocoronavirinae', 
    'Betacoronavirus', 'Sarbecovirus']
/references=[Reference(title='Direct Submission', ...)]
/comment=On Mar 13, 2020 this sequence version replaced MT123293.1.
/structured_comment=OrderedDict([('Assembly-Data', OrderedDict([('Assembly Method', 'SPAdes v. v3.13.0'), 
    ('Sequencing Technology', 'Illumina')]))])
Seq('GGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT...AAA', IUPACAmbiguousDNA())

Example SeqRecord CDS feature:

type: CDS
location: [29551:29668](+)
qualifiers:
    Key: codon_start, Value: ['1']
    Key: gene, Value: ['ORF10']
    Key: product, Value: ['ORF10 protein']
    Key: protein_id, Value: ['QIE07489.1']
    Key: translation, Value: ['MGYINVFAFPFTIYSLLLCRMNSRNYIAQVDVVNFNLT']
'''


def main():
    extractor = Feature_To_Gene_And_Protein(args.format, args.split, args.analyze, 
        args.verbose, args.json, args.cov_dir, args.host_filter, args.files)
    extractor.read()
    extractor.standardize()
    extractor.create_objects()
    extractor.remove_invalid()
    extractor.write()


class Feature_To_Gene_And_Protein:
    from utils import make_sequence_json

    def __init__(self, seq_format, split, analyze, verbose, json, cov_dir, host_filter, files):
        self.seq_format = seq_format
        self.split = split
        # Do not write to any files
        self.analyze = analyze
        self.verbose = verbose
        self.make_json = json
        self.cov_dir = cov_dir
        self.host_filter = host_filter
        self.files = files
        # Initial collection of all features keyed by accession
        self.accs = dict()
        # Features sorted by standard gene/protein names
        self.sorted_cds = dict()
        self.sorted_mats = dict()
        self.json = dict()
        # Nucleotide and protein features for writing to fasta
        self.feats = dict()
        # Synonyms that are matched
        self.matches = dict()
        # Dictionaries
        self.synonyms = read_synonyms()
        self.variants = read_variants()


    def remove_invalid(self):
        '''
        Features may be incorrectly named, or named used older names, we'll call these
        "invalid". The "invalid" sequences are usually different lengths from the "valid"
        sequences and not alignable with the "valid" sequences, though they share the same name. 
        All "valid" sequences are written to files using a standard gene name (e.g. "NS8") and all 
        "invalid" sequences are written to "invalid" files (e.g. "NS8-nt-invalid.fasta").
        '''
        # Copy from self.feats to an empty dict so we don't need to handle changes to self.feats
        feats = dict()

        for gene in self.feats:
            for feat in self.feats[gene].keys():
                feat_len = len(self.feats[gene][feat]['aa'])
                # If a given feature length is not in the right length range
                if feat_len < self.variants[gene][0] or feat_len > self.variants[gene][1]:
                    if self.verbose:
                        print("'{0}' is invalid: {1} not in {2}".format(feat,
                                    feat_len, str(self.variants[gene])))
                    invalid = gene + '-invalid'
                    if invalid not in feats.keys():
                        feats[invalid] = dict()
                    if feat not in feats[invalid].keys():
                        feats[invalid][feat] = dict()
                    # Copy feature
                    feats[invalid][feat]['aa'] = self.feats[gene][feat]['aa']
                    feats[invalid][feat]['nt'] = self.feats[gene][feat]['nt']
                # Check that ORF1a, ORF1ab, NS1 begin with M
                elif gene in ['ORF1a','ORF1ab','NS1'] and str(self.feats[gene][feat]['aa'].seq)[0] != 'M':
                        if self.verbose:
                            print("'{0}' is invalid: no M".format(feat))
                        invalid = gene + '-invalid'
                        if invalid not in feats.keys():
                            feats[invalid] = dict()
                        if feat not in feats[invalid].keys():
                            feats[invalid][feat] = dict()
                        # Copy feature
                        feats[invalid][feat]['aa'] = self.feats[gene][feat]['aa']
                        feats[invalid][feat]['nt'] = self.feats[gene][feat]['nt']
                else:
                    if gene not in feats.keys():
                        feats[gene] = dict()
                    feats[gene][feat] = self.feats[gene][feat]

        self.feats = feats

    def read(self):
        paths = [f for f in self.files if os.path.isfile(f)]
        for path in paths:
            try:
                gbs = [rec for rec in SeqIO.parse(path, "gb")]
            except (RuntimeError) as exception:
                print("Error parsing sequences in '" +
                      str(path) + "':" + str(exception))
            for gb in gbs:
                self.accs[gb.id] = dict()
                self.accs[gb.id]['cds'] = []
                self.accs[gb.id]['mat_peptide'] = []
                # Collect specific features without standard names
                for feat in [feat for feat in gb.features]: 
                    if feat.type == 'CDS':
                        self.accs[gb.id]['cds'].append(feat)
                    if feat.type == 'mat_peptide':
                        self.accs[gb.id]['mat_peptide'].append(feat)
                    # Get 'host' from 'source'
                    if feat.type == 'source':
                        self.accs[gb.id]['source'] = feat
                # Collect the annotations
                self.accs[gb.id]['annotations'] = gb.annotations
                # Collect nucleotide sequence
                self.accs[gb.id]['seq'] = str(gb.seq)

    def standardize(self):
        '''
        1. Look for gene and protein names in 'product'
        2. Look for matches to 'product' in the cov_dictionary
        '''
        self.get_host()
        self.get_date()
        self.get_organism()
        self.standardize_cds()
        self.standardize_mat_peptide()

    def standardize_cds(self):
        '''
        If we get a standard name for a feature from the dictionary we
        rename it and sort it according to the standard name.
        '''
        for acc in self.accs:
            for cds in self.accs[acc]['cds']:
                # Skip a CDS feature without a 'product' tag
                if 'product' not in cds.qualifiers.keys():
                    continue
                # Skip a CDS feature without a 'translation' tag
                if 'translation' not in cds.qualifiers.keys():
                    continue                
                # Skip 'hypothetical protein', 'putative protein', etc.
                if cds.qualifiers["product"][0] in self.synonyms['SKIP']:
                    continue
                id = self.get_standard_name(cds.qualifiers["product"][0], acc)
                if id:
                    if id not in self.sorted_cds.keys():
                        self.sorted_cds[id] = []
                    # For example: "ORF1ab-GU553365.1"
                    cds.id = id + '-' + acc
                    self.sorted_cds[id].append(cds)
    
    def standardize_mat_peptide(self):
        for acc in self.accs:
            for pep in self.accs[acc]['mat_peptide']:
                # Skip a mat_peptide feature without a 'product' tag
                if 'product' not in pep.qualifiers.keys():
                    continue
                # Skip 'hypothetical protein', 'putative protein', etc.
                if pep.qualifiers["product"][0] in self.synonyms['SKIP']:
                    continue
                id = self.get_standard_name(pep.qualifiers["product"][0], acc)
                if id:
                    if id not in self.sorted_mats.keys():
                        self.sorted_mats[id] = []
                    pep.id = id + '-' + acc
                    self.sorted_mats[id].append(pep)

    def get_standard_name(self, product, acc):
        for name in self.synonyms:
            if product in self.synonyms[name]:
                if self.verbose:
                    print("Gene {0} found: {1}, '{2}'".format(name, acc, product))
                return name
        if self.verbose:
            print("No gene for product {0}, '{1}'".format(acc, product))
        if self.analyze:
            match = self.synonyms[name][self.synonyms[name].index(product)]
            self.matches[match] = 1
        return None

    def get_host(self):
        # Iterate over list since we may be deleting keys
        for acc in list(self.accs):
            if 'host' in self.accs[acc]['source'].qualifiers.keys():
                self.accs[acc]['host'] = self.accs[acc]['source'].qualifiers['host'][0]
            elif 'isolation_source' in self.accs[acc]['source'].qualifiers.keys():
                self.accs[acc]['host'] = self.accs[acc]['source'].qualifiers['isolation_source'][0]
            else:
                self.accs[acc]['host'] = 'unknown'
            # Filter by host if specified
            if self.host_filter:
                if self.host_filter not in self.accs[acc]['host']:
                    self.accs.pop(acc)

    def get_date(self):
        for acc in self.accs.keys():
            if 'date' in self.accs[acc]['annotations'].keys():
                self.accs[acc]['date'] = self.accs[acc]['annotations']['date']
            else:
                self.accs[acc]['date'] = 'unknown'

    def get_organism(self):
        for acc in self.accs.keys():
            if 'organism' in self.accs[acc]['annotations'].keys():
                self.accs[acc]['organism'] = self.accs[acc]['annotations']['organism']
            else:
                self.accs[acc]['organism'] = 'unknown'

    def get_location(self, feat):
        for acc in self.accs.keys():
            if 'organism' in self.accs[acc]['annotations'].keys():
                self.accs[acc]['organism'] = self.accs[acc]['annotations']['organism']
            else:
                self.accs[acc]['organism'] = 'unknown'

    def create_objects(self):
        '''
            Create SeqRecords for aa and nt that will written out as fasta.
            Fasta format metadata is made up of "id" and "description".
            All other SeqRecord fields are ignored when Biopython makes fasta.
        '''
        for name in self.sorted_cds.keys():
            self.feats[name] = dict()
            if self.make_json:
                self.json[name] = dict()
            for feat in self.sorted_cds[name]:
                self.feats[name][feat.id] = dict()
                acc = feat.id.split('-')[1]
                # nt
                ntseq = feat.extract(self.accs[acc]['seq'])
                nt = SeqRecord(Seq(ntseq, IUPAC.ambiguous_dna),
                                  id=feat.id,
                                  description=self.make_desc(feat, ntseq, 'nt'))
                self.feats[name][feat.id]['nt'] = nt
                # aa
                aaseq = feat.qualifiers["translation"][0]
                aa = SeqRecord(Seq(aaseq, IUPAC.extended_protein),
                                  id=feat.id,
                                  description=self.make_desc(feat, aaseq, 'aa'))
                self.feats[name][feat.id]['aa'] = aa

                if self.make_json:
                    from utils import make_sequence_json
                    self.json[name][feat.id] = dict()
                    self.json[name][feat.id]['aa'] = make_sequence_json(aa, 'aa')
                    self.json[name][feat.id]['nt'] = make_sequence_json(nt, 'nt')

        for name in self.sorted_mats.keys():
            self.feats[name] = dict()
            if self.make_json:
                self.json[name] = dict()
            for feat in self.sorted_mats[name]:
                self.feats[name][feat.id] = dict()
                acc = feat.id.split('-')[1]
                # nt
                ntseq = feat.extract(self.accs[acc]['seq'])
                nt = SeqRecord(Seq(ntseq, IUPAC.ambiguous_dna),
                                  id=feat.id,
                                  description=self.make_desc(feat, ntseq, 'nt'))
                self.feats[name][feat.id]['nt'] = nt
                # aa
                aaseq = str(nt.translate().seq)
                aa = SeqRecord(Seq(aaseq, IUPAC.extended_protein),
                                  id=feat.id,
                                  description=self.make_desc(feat, aaseq, 'aa'))
                self.feats[name][feat.id]['aa'] = aa

                if self.make_json:
                    from utils import make_sequence_json
                    self.json[name][feat.id] = dict()
                    self.json[name][feat.id]['aa'] = make_sequence_json(aa, 'aa')
                    self.json[name][feat.id]['nt'] = make_sequence_json(nt, 'nt')

    def make_desc(self, feat, seq, seqtype):
        '''
        The description will have:
        - length ('aa' or 'nt')
        - nucleotide coordinates (all features)
        - date (annotation)
        - host (source feature)
        - organism (annotation)
        - protein id (CDS feature only)
        '''
        slen = str(len(seq)) + seqtype
        acc = feat.id.split('-')[1]
        loc = str(feat.location).replace(' ','').replace('join{','').replace('}','').replace('(+)','')
        desc = " ".join([slen,
                        loc,
                        self.accs[acc]['date'],
                        "host='{}'".format(self.accs[acc]['host']),
                        "org='{}'".format(self.accs[acc]['organism']) ])
        if 'protein_id' in feat.qualifiers.keys():
            desc = desc + ' pid=' + feat.qualifiers["protein_id"][0]
        return desc

    def write(self):
        if self.analyze:
            return
        if self.split:
            return
            # for acc in self.accs:
            #     for feat in self.accs[acc]:
            #         seqfile = acc + '.' + self.seq_format
            # SeqIO.write(record, seqfile, self.seq_format)
            # seqfile = record + '-CDS' + '.' + self.seq_format
            # SeqIO.write(seqs, seqfile, self.seq_format)
        else:
            for name in self.feats.keys():
                aaseqfile = os.path.join(self.cov_dir, name + '-aa.fa')
                ntseqfile = os.path.join(self.cov_dir, name + '-nt.fa')
                aahandle = open(aaseqfile, "w")
                nthandle = open(ntseqfile, "w")
                for feat in self.feats[name].keys():
                    SeqIO.write(self.feats[name][feat]['aa'], aahandle, self.seq_format)
                    SeqIO.write(self.feats[name][feat]['nt'], nthandle, self.seq_format)
                aahandle.close()
                nthandle.close()
                if self.make_json and 'invalid' not in name:
                    # No JSON for invalid sequences
                    aafile = os.path.join(self.cov_dir, name + '-aa-fasta.json')
                    ntfile = os.path.join(self.cov_dir, name + '-nt-fasta.json')
                    aahandle = open(aafile, "w")
                    nthandle = open(ntfile, "w")
                    for feat in self.feats[name].keys():
                        aahandle.write(self.json[name][feat]['aa'])
                        nthandle.write(self.json[name][feat]['nt'])
                    aahandle.close()
                    nthandle.close()

if __name__ == "__main__":
    main()
