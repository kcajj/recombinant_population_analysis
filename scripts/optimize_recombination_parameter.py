import numpy as np
from viterbi import viterbi_algorithm
import csv
import sys
from multiprocessing import Pool
from itertools import repeat
from array_compression import compress_array, decompress_array, retrive_compressed_array_from_str
from collections import defaultdict
from matplotlib import pyplot as plt


def build_matrix(input):
    matrix = []
    rows = input.split("/")
    for i in rows:
        row = i.split(",")
        float_row = [float(r) for r in row]
        matrix.append(float_row)
    return np.array(matrix)


def get_evidence_arrays(evidences_file):
    csv.field_size_limit(sys.maxsize)

    ancestral_names = []
    evidence_arrays = []
    mapping_starts = []
    mapping_ends = []

    c_reads = 0
    with open(evidences_file) as file:
        tsv_file = csv.reader(file, delimiter="\t")
        for line in tsv_file:
            ancestral_name = line[0]
            mapping_start = int(line[1])
            mapping_end = int(line[2])

            compressed_evidence_array = retrive_compressed_array_from_str(line[3])
            evidence_array = decompress_array(compressed_evidence_array)

            ancestral_names.append(ancestral_name)
            evidence_arrays.append(evidence_array)
            mapping_starts.append(mapping_start)
            mapping_ends.append(mapping_end)
            c_reads += 1

    return ancestral_names, evidence_arrays, mapping_starts, mapping_ends, c_reads


def write_prediction_arrays(output_path, results, read_names, mapping_starts, mapping_ends):
    with open(output_path, "w", newline="") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t", lineterminator="\n")
        for i in range(len(results)):
            hmm_prediction = results[i][0]
            log_lik = results[i][1]
            read_name = read_names[i]
            mapping_start = mapping_starts[i]
            mapping_end = mapping_ends[i]

            compressed_hmm_prediction = compress_array(hmm_prediction)

            np.set_printoptions(threshold=np.inf, linewidth=np.inf)
            writer.writerow([read_name, mapping_start, mapping_end, log_lik, compressed_hmm_prediction])


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="Makes a prediction on each evidence array and summarises the information in a single recombination array",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--replicates", help="replicates names")
    parser.add_argument("--timesteps", help="list of timesteps")
    parser.add_argument("--evidences", help="path of the folder containing the evidence arrays")
    parser.add_argument("--out", help="output path of the plot with the optimization")
    parser.add_argument("--cores", help="number of cores to use", type=int)
    parser.add_argument("--initial_p", help="initial probabilities of the HMM states")
    parser.add_argument("--transition_p", help="transition probabilities of the HMM states")
    parser.add_argument("--emission_p", help="emission probabilities of the HMM states")
    parser.add_argument("--subsample", help="number of reads to subsample", type=int)

    args = parser.parse_args()
    replicates = args.replicates.split(",")[:-1]
    timesteps = args.timesteps.split(",")[:-1]
    evidences_folder = args.evidences
    output_path = args.out
    cores = args.cores
    initial_probability = args.initial_p
    transition_probabilities = [float(prob) for prob in args.transition_p.split(",")]
    emission_probability = args.emission_p
    subsample = args.subsample

    initial_probability_matrix = build_matrix(initial_probability)
    emission_probability_matrix = build_matrix(emission_probability)
    transition_probabilities.sort()

    log_liks = defaultdict(dict)  # log likelihoods saved for each clone and population
    """
    keys: probability of recombination
    values: keys: population_clone
            values: log likelihood
    """
    for replicate in replicates:
        for timestep in timesteps:

            evidences_file = f"{evidences_folder}/{replicate}/{timestep}.tsv"

            read_names, evidence_arrays, mapping_starts, mapping_ends, c_reads = get_evidence_arrays(evidences_file)

            if subsample<c_reads:
                idx = np.random.choice(c_reads, subsample, replace=False)
                read_names = [read_names[i] for i in idx]
                evidence_arrays = [evidence_arrays[i] for i in idx]
                mapping_starts = [mapping_starts[i] for i in idx]
                mapping_ends = [mapping_ends[i] for i in idx]

            for prob in transition_probabilities:

                transition_probability_matrix = np.array([[1 - prob, prob], [prob, 1 - prob]])

                with Pool(cores) as p:
                    results = p.starmap(
                        viterbi_algorithm,
                        zip(
                            evidence_arrays,
                            repeat(transition_probability_matrix),
                            repeat(emission_probability_matrix),
                            repeat(initial_probability_matrix),
                        ),
                    )

                tot_log_lik = 0
                for i in range(len(results)):
                    hmm_prediction = results[i][0]
                    log_lik = results[i][1]
                    tot_log_lik += log_lik

                log_liks[prob][f"{replicate}_{timestep}"] = tot_log_lik

    mean_log_liks = []
    for prob, samples in log_liks.items():
        prob_mean = []
        for sample, log_lik in samples.items():
            prob_mean.append(log_lik)
        mean_log_liks.append(np.mean(prob_mean))

    plt.plot(transition_probabilities, mean_log_liks)
    plt.xlabel("Recombination probability")
    plt.ylabel("log likelihood")
    plt.title("Optimization of the recombination parameter by log likelihood maximisation")
    plt.savefig(output_path, bbox_inches="tight")
