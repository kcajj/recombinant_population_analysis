from Bio import SeqIO,AlignIO
import numpy as np
import random

if __name__ == "__main__":

    msa_path = "results/msa/refs_msa.fasta"
    hybrid_ref_path = "results/msa/hybrid_ref.fasta"

    alignment = AlignIO.read(open(msa_path), "fasta")
    print(alignment)
    l=alignment.get_alignment_length()
    msa_matrix=np.zeros([2,l],dtype=str)
    for i,record in enumerate(alignment):
        for pos,nuc in enumerate(record.seq):
            msa_matrix[i][pos]=nuc
    
    hybrid_ref_seq=''
    for pos in range(l):
        if msa_matrix[0][pos]=='-':
            hybrid_ref_seq+=msa_matrix[1][pos]
        elif msa_matrix[1][pos]=='-':
            hybrid_ref_seq+=msa_matrix[0][pos]
        else:
            if random.random()<0.5:
                hybrid_ref_seq+=msa_matrix[0][pos]
            else:
                hybrid_ref_seq+=msa_matrix[1][pos]

    with open(hybrid_ref_path, "w") as output_handle:
        output_handle.write(">hybrid_ref\n")
        output_handle.write(hybrid_ref_seq)