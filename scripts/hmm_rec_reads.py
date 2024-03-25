import numpy as np
from handle_msa import get_evidences_distributions, add_to_msa, length_msa
from viterbi import viterbi_algorithm
import pysam
from collections import defaultdict
import time
    
initial_p={"A":0.5,"B":0.5}

transition_p_froma={"A":0.999,"B":0.001}
transition_p_fromb={"A":0.001,"B":0.999}

emission_p_froma={".":0.969,"a":0.03,"b":0.001}
emission_p_fromb={".":0.969,"a":0.001,"b":0.03}

ip_np=np.array(list(initial_p.values()))
tp_np=np.array([list(transition_p_froma.values()),
                list(transition_p_fromb.values())])
ep_np=np.array([list(emission_p_froma.values()),
                list(emission_p_fromb.values())])

'''
takes a MSA, extracts the evidences (handle msa) and gives them in input to viterbi algorithm
'''

if __name__ == "__main__":

    population='P2'
    timestep='7'

    phages={'EM11':0,'EM60':1}

    refs_msa_path="results/msa/refs_msa.fasta"

    bam_file=f"results/alignments/{population}_{timestep}.bam"

    output_path=f"results/genomewide_recombination_arrays/{population}_{timestep}.npz"

    length_treshold=5000

    l_msa=length_msa(refs_msa_path)
    recombination_distribution=np.zeros(l_msa,dtype=int)

    time_spent_per_read=defaultdict(list)
    time_spent_per_base=defaultdict(list)
    tot_t_start=time.time()
    c=0
    c1=0
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if not(read.is_secondary):
                read_seq = read.query_sequence
                read_msa_seq = ''

                l_read=read.infer_query_length()
                if l_read>length_treshold:

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

                    hmm_prediction = viterbi_algorithm(e_distribution, tp_np, ep_np, ip_np)

                    pre_status = hmm_prediction[0]
                    for i in range(len(hmm_prediction)):
                        post_status = hmm_prediction[i]
                        if pre_status != post_status:
                            recombination_distribution[i] += 1
                        pre_status = post_status

                    end_time=time.time()
                    time_spent_per_read[population].append(end_time-start_time)
                    time_spent_per_base[population].append((end_time-start_time)/l_read)

                    print(c)
                    c1+=1
            c+=1
            if c%1000==0:
                print("mean time spent (per read and per base)")
                for k,v in time_spent_per_read.items():
                    print(k," ",np.mean(v))
                for k,v in time_spent_per_base.items():
                    print(k," ",np.mean(v))
                print("")
                print("total reads", c)
                print("reads used", c1)

                np.savez(output_path,recombination_distribution)

print("mean time spent (per read and per base)")
for k,v in time_spent_per_read.items():
    print(k," ",np.mean(v))
for k,v in time_spent_per_base.items():
    print(k," ",np.mean(v))
print("total reads", c)
print("reads used", c1)

np.savez(output_path,recombination_distribution)

tot_t_finish=time.time()
tot_t=tot_t_finish-tot_t_start
print("total time", tot_t)