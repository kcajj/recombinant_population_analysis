from Bio import AlignIO
import numpy as np
import subprocess
from handle_msa import get_evidences_distributions, add_to_msa
import time
from collections import defaultdict
import pysam

phages={'EM11':0,'EM60':1}

refs_msa_path="results/msa/refs_msa.fasta"

null_prob=defaultdict(list)
a_prob=defaultdict(list)
b_prob=defaultdict(list)
time_spent=defaultdict(list)

for phage in phages:

    bam_file=f"data/test/hybrid_test_{phage}_new_chemistry.bam" #for test dataset
    
    c = 0
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if not(read.is_secondary):

                read_seq = read.query_sequence
                read_msa_seq = ''

                mapping_start=read.reference_start
                mapping_end=read.reference_end

                start_time=time.time()

                alignment_array = read.get_aligned_pairs()

                for (read_pos,ref_pos) in alignment_array:
                    if ref_pos!=None:
                        if read_pos==None:
                            read_msa_seq+='-'
                        else:
                            read_msa_seq+=read_seq[read_pos].lower()

                msa_matrix = add_to_msa(refs_msa_path, read_msa_seq, mapping_start, mapping_end)

                e_distribution_to_plot = get_evidences_distributions(msa_matrix)

                e_distribution = np.where(e_distribution_to_plot > 0, e_distribution_to_plot-1, e_distribution_to_plot)

                l = len(e_distribution)
                null_prob[phage].append(np.count_nonzero(e_distribution == 0)/l)
                a_prob[phage].append(np.count_nonzero(e_distribution == 1)/l)
                b_prob[phage].append(np.count_nonzero(e_distribution == 2)/l)
                
                end_time=time.time()
                time_spent[phage].append((end_time-start_time)/l)
                    
                c+=1

                print(c)

print("mean null probability")
for k,v in null_prob.items():
    print(k," ",np.mean(v))
print("")

print("mean a probability")
for k,v in a_prob.items():
    print(k," ",np.mean(v))
print("")

print("mean b probability")
for k,v in b_prob.items():
    print(k," ",np.mean(v))
print("")

print("mean time spent per base")
for k,v in time_spent.items():
    print(k," ",np.mean(v))
print("")
