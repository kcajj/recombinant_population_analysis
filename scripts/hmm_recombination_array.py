import numpy as np
from handle_msa import length_msa
from viterbi import viterbi_algorithm
import time
import csv
import sys
from multiprocessing import Pool
from itertools import repeat

def build_matrix(input):
    matrix=[]
    rows=input.split("/")
    for i in rows:
        row=i.split(",")
        float_row=[float(r) for r in row]
        matrix.append(float_row)
    return np.array(matrix)

def get_evidence_arrays(evidences_file):
    csv.field_size_limit(sys.maxsize)

    evidence_arrays=[]
    mapping_starts=[]

    c_reads=0
    with open(evidences_file) as file:
        tsv_file = csv.reader(file, delimiter="\t")
        for line in tsv_file:
            #read_name=line[0]
            mapping_start=int(line[1])
            #mapping_end=int(line[2])
            evidence_array_str=line[3][1:-1].split(' ')
            evidence_array=np.array([int(e) for e in evidence_array_str])

            evidence_arrays.append(evidence_array)
            mapping_starts.append(mapping_start)
            c_reads+=1

    return evidence_arrays, mapping_starts, c_reads

def get_recombination_array(results, mapping_starts, refs_msa_path):

    recombination_distribution=np.zeros(length_msa(refs_msa_path),dtype=int)
    tot_log_lik=0

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

    return recombination_distribution, tot_log_lik

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="Makes a prediction on each evidence array and summarises the information in a single recombination array",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--evidences", help="path of the .tsv file containing the evidence arrays")
    parser.add_argument("--msa_refs", help="path of the msa between the references")
    parser.add_argument("--out", help="output path of the .npz file containing the recombination array")
    parser.add_argument("--cores", help="number of cores to use", type=int)
    parser.add_argument("--initial_p", help="initial probabilities of the HMM states")
    parser.add_argument("--transition_p", help="transition probabilities of the HMM states")
    parser.add_argument("--emission_p", help="emission probabilities of the HMM states")

    args = parser.parse_args()
    evidences_file=args.evidences
    refs_msa_path=args.msa_refs
    output_path=args.out
    cores=args.cores
    initial_probability=args.initial_p
    transition_probability=args.transition_p
    emission_probability=args.emission_p

    tot_t_start=time.time()

    initial_probability_matrix=build_matrix(initial_probability)
    transition_probability_matrix=build_matrix(transition_probability)
    emission_probability_matrix=build_matrix(emission_probability)

    evidence_arrays, mapping_starts, c_reads = get_evidence_arrays(evidences_file)

    with Pool(cores) as p:
        results = p.starmap(viterbi_algorithm, zip(evidence_arrays, repeat(transition_probability_matrix), repeat(emission_probability_matrix), repeat(initial_probability_matrix)))

    recombination_distribution, tot_log_lik = get_recombination_array(results, mapping_starts, refs_msa_path)

    np.savez(output_path,recombination_distribution)

    tot_t=time.time()-tot_t_start

    output_stats_path=output_path[:-4]+"_stats.txt"
    with open(output_stats_path, 'w') as f:
        f.write("recombination array prediction run of "+evidences_file+'\n')
        f.write("total log likelihood "+str(tot_log_lik)+'\n')
        f.write("total time "+str(tot_t)+'\n')
        f.write("total reads "+str(c_reads)+'\n')
        f.write("time per read "+str(tot_t/c_reads)+'\n')