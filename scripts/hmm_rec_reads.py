from Bio import AlignIO
import numpy as np
import matplotlib.pyplot as plt
from handle_msa import read_msa, get_evidences_distributions, map_refcoord_msacoord, add_to_msa
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

    time_spent=defaultdict(list)

    bam_file=f"data/test/hybrid_test_{population}_{timestep}.bam"

    c=0
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if not(read.is_secondary) and not(read.is_supplementary):

                plot_path = f"results/plots/reads/{population}_{timestep}_{c}.png"

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

                hmm_prediction = viterbi_algorithm(e_distribution, tp_np, ep_np, ip_np)

                end_time=time.time()
                time_spent[population].append(end_time-start_time)

                hmm_plot, (evidences, prediction) = plt.subplots(2, 1, figsize=(10, 5))
                hmm_plot.suptitle(f'HMM read {c}, {population}, t{timestep}')

                colours = np.where(e_distribution_to_plot == 0, "green", np.where(e_distribution_to_plot == 1, "red", np.where(e_distribution_to_plot == 2, "orange", "blue")))
                evidences.scatter(range(mapping_start,len(e_distribution_to_plot)+mapping_start), e_distribution_to_plot, c=colours, marker='|', alpha=0.5)
                evidences.set_title('evidence distribution (0:same, 1:err, 2:a, 3:b)')
                evidences.set_xlabel("basepair")
                evidences.set_ylabel("visible states")

                colours = np.where(hmm_prediction == 0, "orange", "blue")
                prediction.scatter(range(mapping_start,len(e_distribution_to_plot)+mapping_start), hmm_prediction, c=colours, marker='|', alpha=0.5)
                prediction.set_title(f'HMM prediction (0:A, 1:B)')
                prediction.set_xlabel("basepair")
                prediction.set_ylabel("hidden states")

                hmm_plot.tight_layout()
                hmm_plot.savefig(plot_path)
                
                print(c)
            c+=1

print("mean time spent")
for k,v in time_spent.items():
    print(k," ",np.mean(v))
print("")
