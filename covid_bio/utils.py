import json
import hashlib
import yaml
import os


def make_genome_json(name):
    # Genomes
    return json.dumps({
        "data_category": 'Annotated Sequence',
        "data_type": 'Complete Genomic Sequence',
        "data_format": 'gb',
        "source": 'NCBI',
        "submitter_id": name + '_gb',
        "file_name": name + '.gb'
    })


def make_sequence_json(obj, seqtype):
    # Sequences
    # Data Category: Protein or Nucleotide
    category = 'Protein' if seqtype == 'aa' else 'Nucleotide'
    return json.dumps({
        "data_category": category,
        "data_type": 'Sequence',
        "data_format": 'fasta',
        "submitter_id": obj.id + '_fasta',
        "file_name": obj.id + '.fasta'
    })


def make_alignment_json(file, aligner):
    # Alignments
    virus_sequence_alignment_id = file.replace('.', '_')
    # Data Category: Protein or Nucleotide
    category = 'Protein' if '-aa' in file else 'Nucleotide'
    return json.dumps({
        "data_category": category,
        "data_type": 'Sequence Alignment',
        "data_format": 'fasta',
        "submitter_id": virus_sequence_alignment_id,
        "file_name": file,
        "alignment_tool": aligner
    })


def make_hmm_json(file):
    # HMMs
    virus_sequence_hmm_id = file.replace('.', '_')
    # Data Category: Protein or Nucleotide
    category = 'Protein' if '-aa' in file else 'Nucleotide'
    return json.dumps({
        "data_category": category,
        "data_type": 'Sequence HMM',
        "data_format": 'hmmer',
        "submitter_id": virus_sequence_hmm_id,
        "file_name": file
    })


def checksum(self, filename):
    with open(filename, "rb") as f:
        bytes = f.read()
    return hashlib.md5(bytes).hexdigest()


def read_strains():
    y = os.path.dirname(os.path.abspath(__file__)) + '/cov_strains.yaml'
    with open(y) as file:
        strains = yaml.load(file, Loader=yaml.FullLoader)
    return strains


def read_synonyms(self):
    y = os.path.dirname(os.path.abspath(__file__)) + '/cov_features.yaml'
    with open(y) as file:
        synonyms = yaml.load(file, Loader=yaml.FullLoader)
    return synonyms


def read_variants(self):
    y = os.path.dirname(os.path.abspath(__file__)) + \
        '/cov_length_variants.yaml'
    with open(y) as file:
        variants = yaml.load(file, Loader=yaml.FullLoader)
    return variants
