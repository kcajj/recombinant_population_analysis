import numpy as np
from viterbi import viterbi_algorithm
import time
import csv
import sys
from multiprocessing import Pool
from itertools import repeat
from array_compression import decompress_array, retrive_compressed_array_from_str, compress_array
from handle_msa import length_seq

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

    read_names=[]
    evidence_arrays=[]
    mapping_starts=[]
    mapping_ends=[]

    c_reads=0
    with open(evidences_file) as file:
        tsv_file = csv.reader(file, delimiter="\t")
        for line in tsv_file:
            read_name=line[0]
            mapping_start=int(line[1])
            mapping_end=int(line[2])
            
            compressed_evidence_array=retrive_compressed_array_from_str(line[3])
            evidence_array=decompress_array(compressed_evidence_array)

            read_names.append(read_name)
            evidence_arrays.append(evidence_array)
            mapping_starts.append(mapping_start)
            mapping_ends.append(mapping_end)
            c_reads+=1

    return read_names, evidence_arrays, mapping_starts, mapping_ends, c_reads

def write_prediction_arrays(output_path, results, read_names, mapping_starts, mapping_ends, refs_msa_path, output_coverage_path):

    coverage = np.zeros(length_seq(refs_msa_path), dtype=int)
    coverage0 = np.zeros(length_seq(refs_msa_path), dtype=int)
    coverage1 = np.zeros(length_seq(refs_msa_path), dtype=int)

    with open(output_path, 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
        for i in range(len(results)):
            hmm_prediction = results[i][0]
            log_lik = results[i][1]
            read_name = read_names[i]
            mapping_start = mapping_starts[i]
            mapping_end = mapping_ends[i]

            #update coverage
            for j in range(len(hmm_prediction)):
                coverage[mapping_start + j] += 1
                if hmm_prediction[j] == 0:
                    coverage0[mapping_start + j] += 1
                else:
                    coverage1[mapping_start + j] += 1

            compressed_hmm_prediction = compress_array(hmm_prediction)

            np.set_printoptions(threshold=np.inf,linewidth=np.inf)
            writer.writerow([read_name, mapping_start, mapping_end, log_lik, compressed_hmm_prediction])
            
    np.savez(output_coverage_path,coverage)
    np.savez(output_coverage_path[:-4]+"_0.npz",coverage0)
    np.savez(output_coverage_path[:-4]+"_1.npz",coverage1)

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="Makes a prediction on each evidence array and summarises the information in a single recombination array",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--evidences", help="path of the .tsv file containing the evidence arrays")
    parser.add_argument("--hybrid_ref", help="path of the hybrid reference")
    parser.add_argument("--predictions_out", help="output path of the .npz file containing the recombination array")
    parser.add_argument("--coverage_out", help="output path of the .npz file containing the coverage array")
    parser.add_argument("--cores", help="number of cores to use", type=int)
    parser.add_argument("--initial_p", help="initial probabilities of the HMM states")
    parser.add_argument("--transition_p", help="transition probabilities of the HMM states")
    parser.add_argument("--emission_p", help="emission probabilities of the HMM states")

    args = parser.parse_args()
    evidences_file=args.evidences
    hybrid_ref_path=args.hybrid_ref
    output_predictions_path=args.predictions_out
    output_coverage_path=args.coverage_out
    cores=args.cores
    initial_probability=args.initial_p
    transition_probability=args.transition_p
    emission_probability=args.emission_p

    tot_t_start=time.time()

    initial_probability_matrix=build_matrix(initial_probability)
    transition_probability_matrix=build_matrix(transition_probability)
    emission_probability_matrix=build_matrix(emission_probability)

    read_names, evidence_arrays, mapping_starts, mapping_ends, c_reads = get_evidence_arrays(evidences_file)

    with Pool(cores) as p:
        results = p.starmap(viterbi_algorithm, zip(evidence_arrays, repeat(transition_probability_matrix), repeat(emission_probability_matrix), repeat(initial_probability_matrix)))

    #write prediction arrays and coverage0 and coverage1 arrays
    write_prediction_arrays(output_predictions_path, results, read_names, mapping_starts, mapping_ends, hybrid_ref_path, output_coverage_path)
    
    tot_t=time.time()-tot_t_start

    output_stats_path=output_predictions_path[:-4]+"_stats.txt"
    with open(output_stats_path, 'w') as f:
        f.write("prediction arrays run of "+evidences_file+'\n')
        f.write("total time "+str(tot_t)+'\n')
        f.write("total reads "+str(c_reads)+'\n')
        if c_reads > 0:
            f.write("time per read "+str(tot_t/c_reads)+'\n')