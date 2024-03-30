import numpy as np
from handle_msa import length_msa
from viterbi import viterbi_algorithm
import time
import csv
import sys
from multiprocessing import Pool
from itertools import repeat

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

if __name__ == "__main__":

    population='P2'
    timestep='7'

    cores=4

    refs_msa_path="results/msa/refs_msa.fasta"

    evidence_file=f"results/evidence_arrays/test_{population}_{timestep}.tsv" #for test dataset

    output_path=f"results/genomewide_recombination_arrays/test_{population}_{timestep}.npz" #for test dataset

    l_msa=length_msa(refs_msa_path)
    recombination_distribution=np.zeros(l_msa,dtype=int)
    
    tot_log_lik=0

    csv.field_size_limit(sys.maxsize)

    evidence_arrays=[]
    mapping_starts=[]

    tot_t_start=time.time()
    c=0
    with open(evidence_file) as file:
        tsv_file = csv.reader(file, delimiter="\t")
        for line in tsv_file:

            mapping_start=int(line[0])
            mapping_end=int(line[1])
            evidence_array_str=line[2][1:-1].replace("\n","").split(" ")
            evidence_array=np.array([int(e) for e in evidence_array_str])

            evidence_arrays.append(evidence_array)
            mapping_starts.append(mapping_start)
            c+=1

    with Pool(cores) as p:
        results = p.starmap(viterbi_algorithm, zip(evidence_arrays, repeat(tp_np), repeat(ep_np), repeat(ip_np)))

    for i in range(len(results)):
        hmm_prediction = results[i][0]
        log_lik = results[i][1]
        mapping_start = mapping_starts[i]

        tot_log_lik+=log_lik

        pre_status = hmm_prediction[0]
        for i in range(len(hmm_prediction)):
            post_status = hmm_prediction[i]
            if pre_status != post_status:
                recombination_distribution[mapping_start+i] += 1
            pre_status = post_status

    print("total log likelihood",tot_log_lik)

    np.savez(output_path,recombination_distribution)

    tot_t_finish=time.time()
    tot_t=tot_t_finish-tot_t_start
    print("total time", tot_t)
    print("total reads", c)
    print("time per read", tot_t/c)