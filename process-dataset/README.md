# Introduction

This directory contains scripts to generate the datasets that DEEPred needs.

The `generate_dataset.py` is used to generate the input datasets. Usage:

`python3 generate_dataset.py <path/to/fasta/file> <path/to/annotations/dir> <feature type> <path/to/output/dir>`

Example:

`python3 generate_dataset.py uniprot_sprot.fasta ../../DEEPred/TrainTestDatasets/ ctriad ./`

Note: Currently only `ctriad` is supported.

## Output

The script will output a file named `dataset.txt`.

In order for DEEPRed to use it as an inpute feature vector it should be copied to `<DEEPRed repo>/FeatureVectors/Parsed_CTriadFeatures_CAFA2.txt` and `<DEEPRed repo>/FeatureVectors/Parsed_CTriadFeatures_uniprot_training_test_set.txt`

# Example

The `example.sh` contains a bash script that downloads the current release of the UniprotKB database in FASTA format and runs the `generate_dataset.py` script.

Note: Need to adjust path to annotations directory accordingly.

# Virtual Environment Requirements

The virtual environment needs to have Python3.9+. It requires the following packages:
- biopython
- ifeatpro

- The virtual environment can be created using the following command:

```
virtualenv <virtenv_name> -p 3.9
. <virtenv_name>/bin/activate
pip3 install requirements.txt
```
