import numpy as np
import matplotlib.pyplot as plt
from handle_msa import read_msa, get_evidences_distributions
from viterbi import viterbi_algorithm
import subprocess
    
initial_p={"A":0.5,"B":0.5}

transition_p_froma={"A":0.99999,"B":0.00001}
transition_p_fromb={"A":0.00001,"B":0.99999}

emission_p_froma={".":0.967,"a":0.03,"b":0.003}
emission_p_fromb={".":0.967,"a":0.003,"b":0.03}

ip_np=np.array(list(initial_p.values()))
tp_np=np.array([list(transition_p_froma.values()),
                list(transition_p_fromb.values())])
ep_np=np.array([list(emission_p_froma.values()),
                list(emission_p_fromb.values())])

'''
takes a MSA, extracts the evidences (handle msa) and gives them in input to viterbi algorithm
'''

if __name__ == "__main__":

    refs_msa_path = f"test/results/msa/msa_refs.fasta"

    populations=['P2','P3']
    clones=['C1','C2','C3','C4']
    references=['EM11','EM60']

    for population in populations:
        for clone in clones:

            temp_total_msa_path = f"test/temp/{population}_{clone}_msa.fasta"

            plot_path = f"test/plots/clones/{population}_{clone}_msa.png"

            assembly_sequence=f"data/clone_assemblies/{population}_{clone}_assembly.fasta"

            # create the alignment
            msa_command = f"mafft --auto --add {assembly_sequence} --keeplength {refs_msa_path} > {temp_total_msa_path}"
            subprocess.run(msa_command, shell=True)

            msa_matrix = read_msa(temp_total_msa_path)

            e_distribution_to_plot = get_evidences_distributions(msa_matrix)

            e_distribution = np.where(e_distribution_to_plot > 0, e_distribution_to_plot-1, e_distribution_to_plot)
            
            hmm_prediction, log_lik= viterbi_algorithm(e_distribution, tp_np, ep_np, ip_np)

            hmm_plot, (evidences, prediction) = plt.subplots(2, 1, figsize=(10, 5))
            hmm_plot.suptitle(f'HMM {population}, {clone}')

            colours = np.where(e_distribution_to_plot == 0, "green", np.where(e_distribution_to_plot == 1, "red", np.where(e_distribution_to_plot == 2, "blue", "orange")))
            evidences.scatter(range(len(e_distribution_to_plot)), e_distribution_to_plot, c=colours, marker='|', alpha=0.5)
            evidences.set_title(f'evidence distribution (0=same, 1=error, 2=evidence for {references[0]}, 3=evidence for {references[1]})')
            evidences.set_xlabel("basepair")
            evidences.set_ylabel("visible states")

            colours = np.where(hmm_prediction == 0, "blue", "orange")
            prediction.scatter(range(len(hmm_prediction)), hmm_prediction, c=colours, marker='|', alpha=0.5)
            prediction.set_title(f'HMM prediction (0={references[0]}, 1={references[1]})')
            prediction.set_xlabel("basepair")
            prediction.set_ylabel("hidden states")

            hmm_plot.tight_layout()
            hmm_plot.savefig(plot_path)

            #remove temporary files
            rm_command = f"rm {temp_total_msa_path}"
            subprocess.run(rm_command, shell=True)