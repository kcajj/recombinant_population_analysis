from Bio import AlignIO
import numpy as np
import matplotlib.pyplot as plt
from handle_msa import read_msa, get_evidences_distributions, map_refcoord_msacoord
from viterbi import viterbi_algorithm
import pysam
from collections import defaultdict
import time
import subprocess
    
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
    maps_refs_msa={}
    for phage in phages:
        maps_refs_msa[phage]=map_refcoord_msacoord(f"data/references/{phage}_assembly.fasta",refs_msa_path,i_ref_in_msa=phages[phage])

    time_spent=defaultdict(list)

    bam_file=f"data/test/test_{population}_{timestep}.bam"

    c=0
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if not(read.is_secondary) and not(read.is_supplementary):

                phage=read.reference_name.split('_')[0]

                map_ref_msa=maps_refs_msa[phage]

                temp_fasta_path = f"results/temp/{c}.fasta"
                temp_total_msa_path = f"results/temp/{population}_{timestep}_{c}_msa.fasta"
                temp_refs_msa_path = f"results/temp/refs_{population}_{timestep}_{c}_msa.fasta"

                plot_path = f"results/plots/MAFFT_reads/{population}_{timestep}_{c}.png"

                read_sequence=read.query_sequence

                start_time=time.time()

                #find the start and end of mapping. index the ref seuqneces with start and end indexes.
                mapping_start=read.reference_start
                mapping_end=read.reference_end-1
                alignment=AlignIO.read(open(refs_msa_path), "fasta")
                cut_alignment=alignment[:,map_ref_msa[mapping_start]:map_ref_msa[mapping_end]]
                AlignIO.write(cut_alignment, temp_refs_msa_path, "fasta")

                # Create a temporary fasta file with the id and sequence of the read
                with open(temp_fasta_path, "w") as temp_fasta:
                    temp_fasta.write(f">{read.query_name}\n{read_sequence}\n")
                    temp_fasta.close()

                # create the alignment
                msa_command = f"mafft --retree 1 --maxiterate 0 --add {temp_fasta_path} --keeplength {temp_refs_msa_path} > {temp_total_msa_path}"
                subprocess.run(msa_command, shell=True)

                msa_matrix = read_msa(temp_total_msa_path)

                print()
                print()
                print()
                print(read.query_name)
                print(mapping_start,mapping_end)
                print(map_ref_msa[mapping_start],map_ref_msa[mapping_end])
                print(c)
                print()
                print()
                print()

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

                #remove temporary files
                rm_command = f"rm {temp_fasta_path} {temp_total_msa_path} {temp_refs_msa_path}"
                subprocess.run(rm_command, shell=True)
                    
            c+=1

print("mean time spent")
for k,v in time_spent.items():
    print(k," ",np.mean(v))
print("")