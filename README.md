# covid-bioinformatics
Software tools to collect and analyze Coronavirus sequences, including code that extracts gene and protein
sequences from COV (coronavirus) genome files downloaded from NCBI, and creates gene and protein-specific 
sequence collections, alignments, and Hidden Markov Models (HMMs). The HMMs can be used to annotate COV
sequences and create BED files for genome visualization. A *cov_strains* file describes the expected genes
for SARS, MERS, and COV2 so that collections are strain-specific.

# example usage
* `./download_gb_by_taxid.py`
* `./feature_to_gene_and_protein.py COV2/*.gb`
* `./seqs_to_aligns_and_hmms.py COV2/*.fa`
* `./annotate_to_bed.py COV2/*.gb`

The COV genes and proteins are parsed from the GenBank files as features and assigned standard names based on 
their *product* or *mat_peptide* tags. Possible synonyms for these standard names are listed in *cov_dictionary.yaml*. A
QC step compares all the COV protein sequences to expected lengths listed in the *cov_length_variants.yaml* file.

# requirements
* Python3 and packages, including [Biopython](https://biopython.org/)
* Sequence aligner ([muscle](https://drive5.com/muscle/), or [clustalo](http://www.clustal.org/omega/), or [mafft](https://mafft.cbrc.jp/alignment/software/))
* [HMMER 3.3](http://hmmer.org)
* To make BED files:
    * [Infernal](http://eddylab.org/infernal/) and [RFAM COV](https://xfam.wordpress.com/2020/04/27/rfam-coronavirus-release/)
    * [tmhmm.py](https://github.com/dansondergaard/tmhmm.py)

# recommended for best performance
* Get an [NCBI API key](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/)
* Pass the key as a command-line argument or configure as an `env` variable (e.g. `export NCBI_API_KEY=8cc3fffffff2b4444492e68a8167aaaa08`)

# output directory
The code is capable of creating large numbers of files. In order to keep these files organized you can specify 
a directory with `-data_dir` that should contain strain-specific directories (i.e. *COV2, SARS, MERS*). The 
data dir defaults to the current directory, and the default strain is *COV2*.

# file formats created
* GenBank sequence: `.gb`
* Fasta sequence: `.fa`
* Fasta index file: `.fai`
* Fasta alignment: `.fasta`
* Stockholm alignment: `.sto`
* MAF alignment: `.maf`
* HMMER profile HMM: `.hmm`
* hmmsearch table: `.tblout`
* BED track file: `.bed`


# to do
* HMM-based annotation and QC
* Visualization (e.g. for [IGV](https://igvteam.github.io/igv-webapp/fileFormats.html))
