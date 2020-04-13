# covid-bioinformatics
Software tools to collect and analyze Coronavirus sequences, including Python code to create gene and protein entries 
starting with COV (coronavirus) genome files downloaded from NCBI.

Typical usage:

*./download-gb-by-taxid.py -no-split*

*./feature-to-gene-and-protein.py taxid-694009.gb*

The COV genes and proteins are parsed from the GenBank files and sorted by their *product* tags in the source NCBI files.
Synonyms for different *product* values are listed in *cov_dictionary.yaml*.


# To-do
* Remove duplicate sequences before creating alignments
* Analyze alignment and remove sequences without similarity
* Trim alignments before creating HMMs
* Code to create alignments and HMMs
