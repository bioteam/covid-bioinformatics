#!/bin/bash

wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
python3 generate_dataset.py uniprot_sprot.fasta ../../DEEPred/TrainTestDatasets/ ctriad ./