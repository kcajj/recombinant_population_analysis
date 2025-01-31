import numpy as np
from viterbi import viterbi_algorithm
import time
import csv
import sys
import matplotlib.pyplot as plt

initial_p={"A":0.5,"B":0.5}

emission_p_froma={".":0.967,"a":0.03,"b":0.003}
emission_p_fromb={".":0.967,"a":0.003,"b":0.03}

ip_np=np.array(list(initial_p.values()))
ep_np=np.array([list(emission_p_froma.values()),
                list(emission_p_fromb.values())])

transition_p=[0.000001,0.00001,0.00005,0.0001,0.0002]

if __name__ == "__main__":

    population='P2'
    timestep='7'

    refs_msa_path="results/msa/refs_msa.fasta"

    evidence_file=f"test/results/evidence_arrays/test_{population}_{timestep}.tsv"

    output_path=f"test/plots/parameter_estimation/test_{population}_{timestep}.png"

    csv.field_size_limit(sys.maxsize)

    log_likelihoods=[]
    for tp in transition_p:

        tp_np=np.array([[1-tp,tp],[tp,1-tp]])
        tot_log_lik=0
        c=0
        with open(evidence_file) as file:
            tsv_file = csv.reader(file, delimiter="\t")
            for line in tsv_file:

                mapping_start=int(line[0])
                mapping_end=int(line[1])
                evidence_array_str=line[2][1:-1].replace("\n","").split(" ")
                evidence_array=np.array([int(e) for e in evidence_array_str])

                start_time=time.time()

                hmm_prediction,log_lik = viterbi_algorithm(evidence_array, tp_np, ep_np, ip_np)

                tot_log_lik+=log_lik

                c+=1

        print(tp)
        log_likelihoods.append(tot_log_lik)
        print("total log likelihood",tot_log_lik)

    print(log_likelihoods)
    print(transition_p)
    plt.plot(['0.000001','0.00001','0.00005','0.0001','0.0002'],log_likelihoods)
    plt.savefig(output_path)