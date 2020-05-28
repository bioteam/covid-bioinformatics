#!/usr/bin/env python3

import argparse
import sys
import os
import re

parser = argparse.ArgumentParser()
parser.add_argument('-verbose', default=False, action='store_true', help="Verbose")
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
    extractor = Fast_Uniprot_Parser(args.verbose, args.files)
    extractor.read()
    # extractor.write()


class Fast_Uniprot_Parser:

    def __init__(self, verbose, files):
        self.verbose = verbose
        self.files = files
        self.data = dict()

    def read(self):
        paths = [f for f in self.files if os.path.isfile(f)]
        for path in paths:
            with open(path, 'r') as f:
                for linenum,line in enumerate(f):
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
                    # Only collect sequence if there's a GO term
                    if 'GO' in data.keys():
                        matches = re.match(r'^     ([\sA-Z]+)', line)
                        if matches:
                            data['SQ'].append(matches[1].replace(' ','').strip())
                            continue
                    matches = re.match(r'^//', line)
                    if matches:
                        if 'GO' in data.keys():
                            data['SQ'] = ''.join(data['SQ'])
                            self.data[pid] = data
                    if self.verbose and linenum % 10000 == 0:
                        print("Line {}".format(linenum))

                    
            print(self.data)


    def write(self):
        for name in self.feats.keys():
                aaseqfile = name + '-aa.fa'
                ntseqfile = name + '-nt.fa'
                aahandle = open(aaseqfile, "w")
                nthandle = open(ntseqfile, "w")
                for feat in self.feats[name].keys():
                    SeqIO.write(self.feats[name][feat]['aa'], aahandle, self.seq_format)
                    SeqIO.write(self.feats[name][feat]['nt'], nthandle, self.seq_format)
                aahandle.close()
                nthandle.close()
                if self.make_json and 'invalid' not in name:
                    # No JSON for invalid sequences
                    aafile = name + '-aa-fasta.json'
                    ntfile = name + '-nt-fasta.json'
                    aahandle = open(aafile, "w")
                    nthandle = open(ntfile, "w")
                    for feat in self.feats[name].keys():
                        aahandle.write(self.json[name][feat]['aa'])
                        nthandle.write(self.json[name][feat]['nt'])
                    aahandle.close()
                    nthandle.close()
                    
    
if __name__ == "__main__":
    main()
