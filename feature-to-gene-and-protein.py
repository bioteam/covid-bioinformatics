#!/usr/local/bin/python3

import argparse
import os
import yaml
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

parser = argparse.ArgumentParser()
parser.add_argument('-f', default='fasta', dest='format', help="Output format")
parser.add_argument('-split', default=False, dest='split', help="Split into separate files")
parser.add_argument('-no-split', action='store_false', dest='split', help="Create one file")
parser.add_argument('-analyze', default=False, action='store_true', help="Parse files without file creation")
parser.add_argument('-verbose', default=False, action='store_true', help="Verbose")
parser.add_argument('files', nargs='+', help='File names')
args = parser.parse_args()

'''
SeqRecord:

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
/organism=Severe acute respiratory syndrome coronavirus 2
/taxonomy=['Viruses', 'Riboviria', 'Nidovirales', 'Cornidovirineae', 'Coronaviridae', 'Orthocoronavirinae', 'Betacoronavirus', 'Sarbecovirus']
/references=[Reference(title='Direct Submission', ...)]
/comment=On Mar 13, 2020 this sequence version replaced MT123293.1.
/structured_comment=OrderedDict([('Assembly-Data', OrderedDict([('Assembly Method', 'SPAdes v. v3.13.0'), ('Sequencing Technology', 'Illumina')]))])
Seq('GGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT...AAA', IUPACAmbiguousDNA())

SeqRecord Features:

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
    extractor = Feature_To_Gene_And_Protein(args.format, args.split, args.analyze, args.verbose, args.files)
    extractor.read()
    extractor.standardize()
    extractor.create_objects()
    extractor.write()


class Feature_To_Gene_And_Protein:

    def __init__(self, seq_format, split, analyze, verbose, files):
        self.seq_format = seq_format
        self.split = split
        self.analyze = analyze
        self.verbose = verbose
        self.files = files
        # Initial collection of features
        self.accs = dict()
        # Features sorted by standard gene/protein names
        self.sorted_cds = dict()
        self.sorted_mats = dict()
        # Nucleotide and protein features for writing to fasta
        self.nt = dict()
        self.aa = dict()
        # Features that cannot be identified by name
        self.unidentified = []
        self.synonyms = self.read_synonyms()

    def read_synonyms(self):
        y = os.path.dirname(os.path.abspath(__file__)) + '/cov_dictionary.yaml'
        with open(y) as file:
            synonyms = yaml.load(file, Loader=yaml.FullLoader)
        return synonyms

    def read(self):
        full_paths = [os.path.join(os.getcwd(), path) for path in self.files]
        for path in full_paths:
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
                    print("Product {0} from {1} found".format(product, acc))
                return name
        if self.verbose:
            print("Product {0} from {1} not found".format(product, acc))
        self.unidentified.append(product + '-' + acc)
        return None

    def get_host(self):
        for acc in self.accs.keys():
            if 'host' in self.accs[acc]['source'].qualifiers.keys():
                self.accs[acc]['host'] = self.accs[acc]['source'].qualifiers['host'][0]
            elif 'isolation_source' in self.accs[acc]['source'].qualifiers.keys():
                self.accs[acc]['host'] = self.accs[acc]['source'].qualifiers['isolation_source'][0]
            else:
                self.accs[acc]['host'] = 'unknown'

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
        if self.analyze:
            return
        '''
            Create SeqRecords for aa and nt that will written out as fasta.
            Fasta format metadata is made up of "id" and "description".
            All other SeqRecord fields are ignored when Biopython makes fasta.
        '''
        for name in self.sorted_cds.keys():
            if name not in self.aa.keys():
                self.aa[name] = []
            if name not in self.nt.keys():
                self.nt[name] = []
            for feat in self.sorted_cds[name]:
                desc = self.make_desc(feat)
                acc = feat.id.split('-')[1]
                # nt
                ntseq = SeqRecord(Seq(feat.extract(self.accs[acc]['seq']), IUPAC.ambiguous_dna),
                                  id=feat.id,
                                  description=desc)
                self.nt[name].append(ntseq)
                # aa
                aaseq = SeqRecord(Seq(feat.qualifiers["translation"][0], IUPAC.extended_protein),
                                  id=feat.id,
                                  description=desc)
                self.aa[name].append(aaseq)

        for name in self.sorted_mats.keys():
            if name not in self.aa.keys():
                self.aa[name] = []
            if name not in self.nt.keys():
                self.nt[name] = []
            for feat in self.sorted_mats[name]:
                desc = self.make_desc(feat)
                acc = feat.id.split('-')[1]
                # nt
                ntseq = SeqRecord(Seq(feat.extract(self.accs[acc]['seq']), IUPAC.ambiguous_dna),
                                  id=feat.id,
                                  description=desc)
                self.nt[name].append(ntseq)
                # aa
                aaseq = SeqRecord(Seq(str(ntseq.translate().seq), IUPAC.extended_protein),
                                  id=feat.id,
                                  description=desc)
                self.aa[name].append(aaseq)

    def make_desc(self, feat):
        '''
        The description will have:
        - nucleotide coordinates (all features)
        - date (annotation)
        - host (source feature)
        - organism (annotation)
        - protein id (CDS feature only)
        '''
        acc = feat.id.split('-')[1]
        loc = str(feat.location).replace(' ','').replace('join{','').replace('}','').replace('(+)','')
        desc = " ".join([loc,
                        self.accs[acc]['date'],
                        "host='{}'".format(self.accs[acc]['host']),
                        "org='{}'".format(self.accs[acc]['organism']) ])
        if 'protein_id' in feat.qualifiers.keys():
            desc = desc + ' ' + feat.qualifiers["protein_id"][0]
        if self.verbose:
            print("Description: {}".format(desc))
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
            for name in self.aa.keys():
                seqfile = name + '-aa.' + self.seq_format
                SeqIO.write(self.aa[name], seqfile, self.seq_format)
            for name in self.nt.keys():
                seqfile = name + '-nt.' + self.seq_format
                SeqIO.write(self.nt[name], seqfile, self.seq_format)

if __name__ == "__main__":
    main()
