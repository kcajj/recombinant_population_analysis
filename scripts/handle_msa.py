from Bio import SeqIO,AlignIO
import numpy as np
import matplotlib.pyplot as plt

def read_msa(path):
    alignment = AlignIO.read(open(path), "fasta")
    print(alignment)
    l=alignment.get_alignment_length()
    msa_matrix=np.zeros([3,l],dtype=str)
    for i,record in enumerate(alignment):
        for pos,nuc in enumerate(record.seq):
            msa_matrix[i][pos]=nuc
    return msa_matrix

def get_evidences_distributions(msa_matrix,i_ref1,i_ref2,i_extra):
    l=len(msa_matrix[0])
    e_distribution = np.zeros(l, dtype=int)

    for pos,array in enumerate(msa_matrix[0]):
        nuc_extra=msa_matrix[i_extra,pos]
        nuc_first_ref=msa_matrix[i_ref1,pos]
        nuc_second_ref=msa_matrix[i_ref2,pos]
        if nuc_first_ref!='-' and nuc_second_ref!='-':
        #if nuc_extra!='-' and nuc_first_ref!='-' and nuc_second_ref!='-':
            if  nuc_extra==nuc_first_ref and nuc_extra==nuc_second_ref:
                continue
            elif nuc_extra!=nuc_first_ref and nuc_extra!=nuc_second_ref:
                e_distribution[pos]=1
            elif nuc_extra==nuc_first_ref and nuc_extra!=nuc_second_ref:
                e_distribution[pos]=2
            elif nuc_extra!=nuc_first_ref and nuc_extra==nuc_second_ref:
                e_distribution[pos]=3

    return e_distribution

def map_refcoord_msacoord(ref_path,refs_msa_path,i_ref_in_msa):
    map={}
    i=0
    j=0
    ref_seq=SeqIO.read(ref_path, "fasta").seq
    alignment=AlignIO.read(open(refs_msa_path), "fasta")
    while i<len(ref_seq):
        if alignment[i_ref_in_msa,j]!='-':
            map[i]=j
            i+=1
        j+=1
    return map