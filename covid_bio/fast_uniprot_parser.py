#!/usr/bin/env python3

import argparse
import sys
import os
import re
import itertools

parser = argparse.ArgumentParser()
parser.add_argument('-verbose', default=False, action='store_true', help="Verbose")
parser.add_argument('-matrix', default=False, action='store_true', help="Make matrix for ChemProp")
parser.add_argument('files', nargs='+', help='File names')
args = parser.parse_args()

'''
Example Uniprot file:

ID   022L_IIV3               Reviewed;         225 AA.
AC   Q197D8;
DT   16-JUN-2009, integrated into UniProtKB/Swiss-Prot.
DT   11-JUL-2006, sequence version 1.
DT   22-APR-2020, entry version 33.
DE   RecName: Full=Transmembrane protein 022L;
GN   ORFNames=IIV3-022L;
OS   Invertebrate iridescent virus 3 (IIV-3) (Mosquito iridescent virus).
OC   Viruses; Iridoviridae; Betairidovirinae; Chloriridovirus.
OX   NCBI_TaxID=345201;
OH   NCBI_TaxID=7163; Aedes vexans (Inland floodwater mosquito) (Culex vexans).
OH   NCBI_TaxID=42431; Culex territans.
OH   NCBI_TaxID=332058; Culiseta annulata.
OH   NCBI_TaxID=310513; Ochlerotatus sollicitans (eastern saltmarsh mosquito).
OH   NCBI_TaxID=329105; Ochlerotatus taeniorhynchus (Black salt marsh mosquito) (Aedes taeniorhynchus).
OH   NCBI_TaxID=7183; Psorophora ferox.
RN   [1]
RP   NUCLEOTIDE SEQUENCE [LARGE SCALE GENOMIC DNA].
RX   PubMed=16912294; DOI=10.1128/jvi.00464-06;
RA   Delhon G., Tulman E.R., Afonso C.L., Lu Z., Becnel J.J., Moser B.A.,
RA   Kutish G.F., Rock D.L.;
RT   "Genome of invertebrate iridescent virus type 3 (mosquito iridescent
RT   virus).";
RL   J. Virol. 80:8439-8449(2006).
CC   -!- SUBCELLULAR LOCATION: Host membrane {ECO:0000305}; Multi-pass membrane
CC       protein {ECO:0000305}.
DR   EMBL; DQ643392; ABF82052.1; -; Genomic_DNA.
DR   RefSeq; YP_654594.1; NC_008187.1.
DR   GeneID; 4156271; -.
DR   KEGG; vg:4156271; -.
DR   Proteomes; UP000001358; Genome.
DR   GO; GO:0033644; C:host cell membrane; IEA:UniProtKB-SubCell.
DR   GO; GO:0016021; C:integral component of membrane; IEA:UniProtKB-KW.
PE   4: Predicted;
KW   Host membrane; Membrane; Reference proteome; Transmembrane;
KW   Transmembrane helix.
FT   CHAIN           1..225
FT                   /note="Transmembrane protein 022L"
FT                   /id="PRO_0000377944"
FT   TRANSMEM        2..22
FT                   /note="Helical"
FT                   /evidence="ECO:0000255"
FT   TRANSMEM        53..73
FT                   /note="Helical"
FT                   /evidence="ECO:0000255"
FT   TRANSMEM        80..100
FT                   /note="Helical"
FT                   /evidence="ECO:0000255"
FT   TRANSMEM        111..131
FT                   /note="Helical"
FT                   /evidence="ECO:0000255"
FT   TRANSMEM        136..156
FT                   /note="Helical"
FT                   /evidence="ECO:0000255"
FT   TRANSMEM        162..182
FT                   /note="Helical"
FT                   /evidence="ECO:0000255"
SQ   SEQUENCE   225 AA;  25107 MW;  3BD60B1CA8C7D7F5 CRC64;
     MSFVHKLPTF YTAGVGAIIG GLSLRFNGAK FLSDWYINKY NDSVPAWSLQ TCHWAGIALY
     CVGWVTLASV IYLKHRDNSI LKGSILSCIV ISAVWSILEY NQDMFVSNPK LPLISCAMLV
     SSLAALVALK YHIKDIFTIL GAAIIIILAE YVVLPYQRQY NIVDGIGLPL LLLGFFILYQ
     VFSVPNPSTP TGVMVPKPED EWDIEMAPLN HRDRQVPESE LENVK
//

'''


def main():
    extractor = Fast_Uniprot_Parser(args.verbose, args.matrix, args.files)
    extractor.read()
    if extractor.matrix:
        extractor.write_matrix()


class Fast_Uniprot_Parser:

    def __init__(self, verbose, matrix, files):
        self.verbose = verbose
        self.matrix = matrix
        self.files = files
        self.data = dict()

    def read(self):
        paths = [f for f in self.files if os.path.isfile(f)]
        for path in paths:
            with open(path, 'r') as f:
                for linenum,line in enumerate(f):
                    # Capture Uniprot id
                    matches = re.match(r'^ID\s+(\S+)', line)
                    if matches:
                        pid = matches[1]
                        data = dict()
                        continue
                    # Look for GO terms
                    matches = re.match(r'^DR   GO; (GO:\d+)', line)
                    if matches:
                        if 'GO' not in data.keys():
                            data['GO'] = []
                            data['SQ'] = []
                        data['GO'].append(matches[1])
                        continue
                    # Look for KEGG terms
                    matches = re.match(r'^DR   KEGG; ([a-z]+:[^;]+)', line)
                    if matches:
                        if 'KEGG' not in data.keys():
                            data['KEGG'] = []
                            data['SQ'] = []
                        data['KEGG'].append(matches[1])
                        continue
                    # Capture sequence if there's a GO or KEGG term
                    matches = re.match(r'^     ([\sA-Z]+)', line)
                    if matches:
                        if 'SQ' in data.keys():
                            data['SQ'].append(matches[1].replace(' ','').strip())
                            continue
                    # End of entry
                    matches = re.match(r'^//', line)
                    if matches:
                        if 'SQ' in data.keys():
                            data['SQ'] = ''.join(data['SQ'])
                            self.data[pid] = data
                    if self.verbose and linenum % 10000 == 0:
                        print("Line {}".format(linenum))

    def write_matrix(self):
        goterms, uniquegoterms = self.parse_terms('GO')
        keggterms, uniquekeggterms = self.parse_terms('KEGG')
        for pid in self.data.keys():
            ntfile = ''
            nthandle = open(ntfile, "w")
            for feat in self.feats[name].keys():
                nthandle.write(self.json[name][feat]['nt'])
            nthandle.close()
                    
    def parse_terms(self, ontology):
        '''
        >>> data
        {3: {'GO': [5, 9]}, 5: {'KEGG': [7, 4]}}
        >>> [ ids.get('GO') for ids in [trms for trms in data.values() ]]
        [[5, 9], None]  
        '''
        # Array of arrays of all terms plus None
        terms = [ gids.get(ontology) for gids in [pids for pids in self.data.values() ]]
        # Remove None
        terms = [i for i in terms if i] 
        # Merge array of arrays
        terms = list(itertools.chain.from_iterable(terms))
        uniqueterms = set(terms)
        return terms, uniqueterms

if __name__ == "__main__":
    main()
