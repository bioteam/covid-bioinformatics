import os
import sys
parentdir = os.path.dirname(__file__)
sys.path.insert(0,parentdir)

from Bio import SeqIO
from ifeatpro_utils import run_ifeatpro
from pathlib import Path


def build_dataset(fasta_file_path, go_annots_path, feature_type, output_dir):

    # Create a list of all the sequences
    lst_all_test_prot_ids = set()

    for subdir, _, test_fls in os.walk(go_annots_path):
        for test_fl in test_fls:
            if test_fl.startswith("test"):
                test_fl_path = Path(os.path.join(subdir, test_fl))
                all_test_ids = open(test_fl_path, "r")
                lst_all_test_ids = all_test_ids.read().split("\n")
                all_test_ids.close()
                while "" in lst_all_test_ids:
                    lst_all_test_ids.remove("")
                for ln in lst_all_test_ids:
                    lst_all_test_prot_ids.add(ln)

    print("INFO: Found a total of {} proteins".format(len(lst_all_test_prot_ids)))

    seq_records = {}
    # We need to parse the id of the protein since it is e.g. id=sp|Q9R157|ADA18_MOUSE' and we want id=Q9R157
    for seq in SeqIO.parse(fasta_file_path, "fasta"):
        seq_records[seq.id.split("|")[1]] = seq

    # Create a dict dict[prot_id] = str(sequence)
    protein_seqs = {}
    for protein in lst_all_test_prot_ids:
        record = seq_records.get(protein)
        if record:
            protein_seqs[protein] = str(record.seq)
        else:
            print(f"No record found in FASTA for protein: {protein}")

    print(f"INFO: Protein to seq dictionary has {len(protein_seqs)} proteins")

    convert_to_deepred_input_format(output_dir, protein_seqs, feature_type)


def convert_to_deepred_input_format(output_dir, protein_seqs, feature_type):
    with open(output_dir + "/" + "dataset.txt", "w") as output_file:
        for protein_id, seq in protein_seqs.items():
            feature_vec = run_ifeatpro(feature_type, seq, protein_id)
            # Write single line for each protein


if __name__ == "__main__":

    if len(sys.argv) != 5:
        print("Usage: python3 generate_dataset.py <path/to/fasta/file> <path/to/annotations/dir> <feature type> <path/to/output/dir>")
        exit(1)

    fasta_file_path = sys.argv[1]
    go_annots_path = sys.argv[2]
    feature_type = sys.argv[3]
    output_dir = sys.argv[4]

    build_dataset(fasta_file_path, go_annots_path, feature_type, output_dir)
