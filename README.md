# covid-bioinformatics
Software tools to collect and analyze Coronavirus sequences, including code to create gene and protein entries 
starting with COV (coronavirus) genome files downloaded from NCBI, then create gene and protein-specific 
collections, alignments, and HMMs.


# typical usage
* ./download_gb_by_taxid.py
* ./feature_to_gene_and_protein.py *.gb
* ./seqs_to_aligns_and_hmms.py *.fasta


The COV genes and proteins are parsed from the GenBank files as features and assigned standard names based on 
their *product* tags. The synonyms for these standard names are listed in *cov_dictionary.yaml*.


# requirements
* Python3 and packages, including Biopython
* Aligner (muscle, or clustalo, or mafft)
* HMMER 3.3


# to_do
* HMM-based annotation and QC
* Visualization (e.g. for [IGV](https://igvteam.github.io/igv-webapp/fileFormats.html))
