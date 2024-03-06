from Bio import AlignIO
import numpy as np
import matplotlib.pyplot as plt
from handle_msa import read_msa, get_evidences_distributions
from viterbi import viterbi_algorithm
    
initial_p={"A":0.5,"B":0.5}

transition_p_froma={"A":0.999,"B":0.001}
transition_p_fromb={"A":0.001,"B":0.999}

emission_p_froma={".":0.949,"a":0.05,"b":0.001}
emission_p_fromb={".":0.949,"a":0.001,"b":0.05}

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
    clones=['C1','C2','C3','C4']
    references=['EM11','EM60']

    for population in populations:
        for clone in clones:

            file=f'results/msa/clones/{population}_{clone}_msa.fasta'
            out_folder=f'results/plots/clones/{population}_{clone}_msa.png'

            msa_matrix = read_msa(file)

            e_distribution = get_evidences_distributions(msa_matrix,i_ref1=1,i_ref2=2,i_extra=0)
            
            hmm_prediction = viterbi_algorithm(e_distribution, tp_np, ep_np, ip_np)

            plt.subplot(2, 1, 1)
            plt.scatter(range(len(e_distribution)), e_distribution, c=e_distribution, marker='|', alpha=0.5)
            plt.title('evidence distribution (purple no evidence, green evidence for A, yellow evidence for B)')

            plt.subplot(2, 1, 2)
            plt.plot(hmm_prediction)
            plt.title(f'hmm prediction {population}, {clone} (0:A, 1:B)')

            plt.tight_layout()
            plt.savefig(out_folder)
            #plt.show()
            plt.close()