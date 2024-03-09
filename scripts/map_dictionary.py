from Bio import SeqIO,AlignIO


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