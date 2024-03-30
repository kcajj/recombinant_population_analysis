import numpy as np
from handle_msa import length_msa
from viterbi import viterbi_algorithm
from collections import defaultdict
import time
import csv
import sys

initial_p={"A":0.5,"B":0.5}

transition_p_froma={"A":0.9999,"B":0.0001}
transition_p_fromb={"A":0.0001,"B":0.9999}

emission_p_froma={".":0.969,"a":0.03,"b":0.001}
emission_p_fromb={".":0.969,"a":0.001,"b":0.03}

ip_np=np.array(list(initial_p.values()))
tp_np=np.array([list(transition_p_froma.values()),
                list(transition_p_fromb.values())])
ep_np=np.array([list(emission_p_froma.values()),
                list(emission_p_fromb.values())])

if __name__ == "__main__":

    population='P2'
    timestep='7'

    refs_msa_path="results/msa/refs_msa.fasta"

    evidence_file=f"results/evidence_arrays/test_{population}_{timestep}.tsv" #for test dataset

    output_path=f"results/genomewide_recombination_arrays/test_{population}_{timestep}.npz" #for test dataset

    l_msa=length_msa(refs_msa_path)
    recombination_distribution=np.zeros(l_msa,dtype=int)
    
    tot_log_lik=0

    csv.field_size_limit(sys.maxsize)

    time_spent_per_read=defaultdict(list)
    time_spent_per_base=defaultdict(list)
    tot_t_start=time.time()
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

            pre_status = hmm_prediction[0]
            for i in range(len(hmm_prediction)):
                post_status = hmm_prediction[i]
                if pre_status != post_status:
                    recombination_distribution[mapping_start+i] += 1
                pre_status = post_status

            end_time=time.time()
            time_spent_per_read[population].append(end_time-start_time)
            time_spent_per_base[population].append((end_time-start_time)/len(hmm_prediction))

            print(c)
            print("log likelihood",log_lik)
            c+=1

    print("mean time spent (per read and per base)")
    for k,v in time_spent_per_read.items():
        print(k," ",np.mean(v))
    for k,v in time_spent_per_base.items():
        print(k," ",np.mean(v))

    print("total log likelihood",tot_log_lik)

    np.savez(output_path,recombination_distribution)

    tot_t_finish=time.time()
    tot_t=tot_t_finish-tot_t_start
    print("total time", tot_t)