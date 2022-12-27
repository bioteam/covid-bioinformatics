# covid-bioinformatics
Software tools to collect and analyze Coronavirus (COV) sequences, including code that extracts gene and protein
sequences from COV genomes from GenBank, and creates gene and protein-specific 
sequence collections, alignments, and Hidden Markov Models (HMMs) for protein and nucleotide COV sequences.
The code uses the [cov_strains.yaml](https://github.com/bioteam/covid-bioinformatics/blob/master/covid_bio/cov_strains.yaml) file that lists the expected genes for *SARS, MERS*, and *COV2* so that collections are strain-specific.

Code is included that annotates COV genomes using the nucleotide HMMs, [pyTMHMM](https://github.com/bosborne/pyTMHMM),
and [Infernal](http://eddylab.org/infernal/) for RNA motif prediction and creates BED files for genome visualization. 

# example usage
* `mkdir COV2`
* `./download_gb_by_taxid.py`
* `./feature_to_gene_and_protein.py COV2/*.gb`
* `./seqs_to_aligns_and_hmms.py COV2/*.fa`
* `./annotate_to_bed.py -rfam_file ../COV2/cov_allvirus.cm COV2/NC_045512.2.gb` 

The COV genes and proteins are parsed from the GenBank files as features and assigned standard names based on 
their *product* or *mat_peptide* tags. Possible synonyms for these standard names are listed in [cov_features.yaml](https://github.com/bioteam/covid-bioinformatics/blob/master/covid_bio/cov_features.yaml). A QC step compares the COV sequences to expected lengths listed in the [cov_length_variants.yaml](https://github.com/bioteam/covid-bioinformatics/blob/master/covid_bio/cov_length_variants.yaml) file.

# requirements
* Python3 and packages, including [Biopython](https://biopython.org/)
* A sequence aligner, [muscle](https://drive5.com/muscle/), or [clustalo](http://www.clustal.org/omega/), or [mafft](https://mafft.cbrc.jp/alignment/software/)) (default: `clustalo`)
* [HMMER 3.3](http://hmmer.org)

# optional
* To make BED files:
    * [Infernal](http://eddylab.org/infernal/) and [RFAM COV](https://xfam.wordpress.com/2020/04/27/rfam-coronavirus-release/) for RNA motif prediction
    * [pyTMHMM](https://github.com/bosborne/pyTMHMM) for trans-membrane region prediction

# recommended for best performance
* Get an [NCBI API key](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/)
* Pass the key as a command-line argument or configure as an `env` variable (e.g. `export NCBI_API_KEY=8cc3fffffff2b4444492e68a8167aaaa08`)

# output directory
The code is capable of creating large numbers of files. In order to keep these files organized you should specify 
a directory with `-data_dir` that is strain-specific (i.e. *COV2, SARS, MERS*), and the
code will write all files to the strain-specific directory as needed. The data dir location is the current directory
by default, and the default strain is *COV2*.

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
