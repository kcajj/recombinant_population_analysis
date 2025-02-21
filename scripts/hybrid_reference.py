from Bio import AlignIO
import numpy as np
import random


def create_hybrid_ref(msa_path):
    """
    Generate a hybrid reference sequence from an MSA.

    extracts the two sequences and iterates over each position. If one of the
    sequences has a gap at a given position, the nucleotide from the other sequence is used.
    In cases where both sequences have valid nucleotides, one is chosen at random with 50%
    probability.

    Args:
        msa_path (str): The filepath to the two-sequence MSA in FASTA format.

    Returns:
        str: The generated hybrid reference sequence.
    """

    alignment = AlignIO.read(open(msa_path), "fasta")
    l = alignment.get_alignment_length()
    msa_matrix = np.zeros([2, l], dtype=str)
    for i, record in enumerate(alignment):
        for pos, nuc in enumerate(record.seq):
            msa_matrix[i][pos] = nuc

    hybrid_ref_seq = ""
    for pos in range(l):
        if msa_matrix[0][pos] == "-":
            hybrid_ref_seq += msa_matrix[1][pos]
        elif msa_matrix[1][pos] == "-":
            hybrid_ref_seq += msa_matrix[0][pos]
        else:
            if random.random() < 0.5:
                hybrid_ref_seq += msa_matrix[0][pos]
            else:
                hybrid_ref_seq += msa_matrix[1][pos]

    return hybrid_ref_seq


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="Create a hybrid reference sequence from a pairwise alignment",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--msa", help="path of the paiwise alignment")
    parser.add_argument("--out", help="output hybrid reference")

    args = parser.parse_args()
    msa_path = args.msa
    hybrid_ref_path = args.out

    hybrid_ref_seq = create_hybrid_ref(msa_path)

    with open(hybrid_ref_path, "w") as output_handle:
        output_handle.write(">hybrid_ref\n")
        output_handle.write(hybrid_ref_seq)
