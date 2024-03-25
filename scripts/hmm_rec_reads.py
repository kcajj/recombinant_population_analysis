import numpy as np
from handle_msa import get_evidences_distributions, add_to_msa, length_msa
from viterbi import viterbi_algorithm
import pysam
from collections import defaultdict
import time
import matplotlib.pyplot as plt
    
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
    #bam_file=f"data/test/hybrid_test_{population}_{timestep}.bam" #for test dataset

    output_path=f"results/genomewide_recombination_arrays/{population}_{timestep}.npz"
    #output_path=f"results/genomewide_recombination_arrays/test_{population}_{timestep}.npz" #for test dataset

    length_treshold=5000 #just for the total dataset

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

                l_alignment=read.query_alignment_length
                if l_alignment>length_treshold: #just for the total dataset

                    mapping_start=read.reference_start
                    mapping_end=read.reference_end

                    start_time=time.time()

                    alignment_array = read.get_aligned_pairs()

                    for (read_pos,ref_pos) in alignment_array:
                        if ref_pos!=None: #clips and gaps ar
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

                    ''' to create a plot for each read
                    plot_path = f"results/plots/reads/{population}_{timestep}_{c}.png"

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
                    '''

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