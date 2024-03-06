from Bio import AlignIO
import numpy as np
import matplotlib.pyplot as plt
from handle_msa import read_msa, get_evidences_distributions
from viterbi import viterbi_algorithm
    
initial_p={"A":0.5,"B":0.5}

transition_p_froma={"A":0.999,"B":0.001}
transition_p_fromb={"A":0.001,"B":0.999}

emission_p_froma={".":0.60,"a":0.30,"b":0.10}
emission_p_fromb={".":0.60,"a":0.10,"b":0.30}

ip_np=np.array(list(initial_p.values()))
tp_np=np.array([list(transition_p_froma.values()),
                list(transition_p_fromb.values())])
ep_np=np.array([list(emission_p_froma.values()),
                list(emission_p_fromb.values())])

'''
takes a MSA, extracts the evidences (handle msa) and gives them in input to viterbi algorithm
'''

if __name__ == "__main__":

    populations=['P2','P3']
    timesteps=['1','3','5','7']
    reads=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19']
    
    for population in populations:
        for timestep in timesteps:
            for read in reads:
                file=f'results/msa/P2/7/P2_7_{read}_msa.fasta'
                out_folder=f'results/plots/reads/{population}_{timestep}_{read}_msa.png'

                msa_matrix = read_msa(file)

                e_distribution = get_evidences_distributions(msa_matrix,i_ref1=0,i_ref2=1,i_extra=2)

                hmm_prediction = viterbi_algorithm(e_distribution, tp_np, ep_np, ip_np)

                plt.subplot(2, 1, 1)
                plt.scatter(range(len(e_distribution)), e_distribution, c=e_distribution, marker='|', alpha=0.5)
                plt.title('evidence distribution (purple no evidence, green evidence for A, yellow evidence for B)')

                plt.subplot(2, 1, 2)
                plt.plot(hmm_prediction)
                plt.title(f'hmm prediction read {read}, {population}, t{timestep} (0:A, 1:B)')

                plt.tight_layout()
                plt.savefig(out_folder)
                #plt.show()
                plt.close()
