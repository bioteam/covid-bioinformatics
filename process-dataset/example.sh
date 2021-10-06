#!/bin/bash

wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
python3 generate_dataset.py -fasta_file_path uniprot_sprot.fasta -go_annots_path ../../DEEPred/TrainTestDatasets/ -output_dir .