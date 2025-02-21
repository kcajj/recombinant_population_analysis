from Bio import SeqIO, AlignIO
import numpy as np
import matplotlib.pyplot as plt


def read_msa(path):
    alignment = AlignIO.read(open(path), "fasta")
    print(alignment)
    l = alignment.get_alignment_length()
    msa_matrix = np.zeros([3, l], dtype=str)
    for i, record in enumerate(alignment):
        for pos, nuc in enumerate(record.seq):
            msa_matrix[i][pos] = nuc
    return msa_matrix


def get_evidences_distributions(msa_matrix, i_ref1=0, i_ref2=1, i_extra=2):
    """
    Generate a distribution of evidence types at each position in a multiple sequence alignment (MSA).

    This function compares a read sequence (i_extra) with two reference sequences (i_ref1, i_ref2).
    For each position:
        - If all three nucleotides (read, first ref, second ref) match (and are not "-"), the evidence is 0 (no difference).
        - If the read nucleotide differs from both references, the evidence is 1.
        - If the read nucleotide matches the first reference but differs from the second, the evidence is 2.
        - If the read nucleotide matches the second reference but differs from the first, the evidence is 3.
    
    ! Any position where one or more sequences have "-" (gap) remains 0 by default.

    Parameters:
        msa_matrix (numpy.ndarray): A 2D array representing the MSA. Each row is a sequence, and each
            column is a position.
        i_ref1 (int, optional): Integer indicating the position of the first reference sequence in the msa. Defaults to 0.
        i_ref2 (int, optional): Integer indicating the position of the second reference sequence in the msa. Defaults to 1.
        i_extra (int, optional): Integer indicating the position of the read in the msa. Defaults to 2.

    Returns:
        numpy.ndarray: A 1D array of integers, where each position indicates the evidence type (0, 1, 2, or 3).
    """
    l = len(msa_matrix[0])
    e_distribution = np.zeros(l, dtype=int)

    for pos, array in enumerate(msa_matrix[0]):
        nuc_extra = msa_matrix[i_extra, pos]
        nuc_first_ref = msa_matrix[i_ref1, pos]
        nuc_second_ref = msa_matrix[i_ref2, pos]
        # if nuc_first_ref!='-' and nuc_second_ref!='-':
        if nuc_extra != "-" and nuc_first_ref != "-" and nuc_second_ref != "-":
            if nuc_extra == nuc_first_ref and nuc_extra == nuc_second_ref:
                continue
            elif nuc_extra != nuc_first_ref and nuc_extra != nuc_second_ref:
                e_distribution[pos] = 1
            elif nuc_extra == nuc_first_ref and nuc_extra != nuc_second_ref:
                e_distribution[pos] = 2
            elif nuc_extra != nuc_first_ref and nuc_extra == nuc_second_ref:
                e_distribution[pos] = 3
        # elif nuc_first_ref=="-" and nuc_second_ref=="-":
        #    if nuc_extra==nuc_first_ref:
        #        continue
        #    elif nuc_extra!=nuc_first_ref:
        #        e_distribution[pos]=1

    return e_distribution


def add_to_msa(msa_path, seq, mapping_start, mapping_end):
    """
    Add a given sequence to a specified region within an existing multiple sequence alignment (MSA).

    Parameters:
        msa_path (str): 
            Path to the MSA file in FASTA format.
        seq (str): 
            The sequence to be inserted or replaced in the MSA region.
        mapping_start (int): 
            Start index (inclusive) of the region within the MSA to be replaced.
        mapping_end (int): 
            End index (exclusive) of the region within the MSA to be replaced.

    Returns:
        numpy.ndarray: 
            A 2D array (rows representing sequences, columns representing positions) 
            containing the updated MSA slice with the specified region replaced by the
            provided sequence.
    """
    alignment = AlignIO.read(open(msa_path), "fasta")
    l = alignment.get_alignment_length()

    msa_matrix = np.zeros([3, l], dtype=str)
    for i, record in enumerate(alignment):
        for pos, nuc in enumerate(record.seq):
            msa_matrix[i][pos] = nuc

    cut_msa_matrix = msa_matrix[:, mapping_start:mapping_end]
    for pos in range(len(cut_msa_matrix[2])):
        cut_msa_matrix[2][pos] = seq[pos]

    return cut_msa_matrix


def extract_references_names(msa_path):
    alignment = AlignIO.read(open(msa_path), "fasta")
    names = []
    for record in alignment:
        names.append(record.id)
    return names


def length_seq(seq_path):
    seq = SeqIO.read(seq_path, "fasta")
    return len(seq.seq)
