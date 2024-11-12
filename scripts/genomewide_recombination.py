import numpy as np
from handle_msa import length_seq
import csv
import sys
from array_compression import decompress_array, retrive_compressed_array_from_str


def get_recombination_array(predictions_file, hybrid_ref_path):

    csv.field_size_limit(sys.maxsize)

    recombination_distribution = np.zeros(length_seq(hybrid_ref_path), dtype=int)
    recombination_distribution_01 = np.zeros(length_seq(hybrid_ref_path), dtype=int)
    recombination_distribution_10 = np.zeros(length_seq(hybrid_ref_path), dtype=int)
    tot_log_lik = 0

    with open(predictions_file) as file:
        tsv_file = csv.reader(file, delimiter="\t")
        for line in tsv_file:
            # read_name=line[0]
            mapping_start = int(line[1])
            # mapping_end=int(line[2])
            log_lik = float(line[3])

            compressed_prediction_array = retrive_compressed_array_from_str(line[4])
            prediction_array = decompress_array(compressed_prediction_array)

            tot_log_lik += log_lik

            pre_status = prediction_array[0]
            for i in range(len(prediction_array)):
                post_status = prediction_array[i]
                if pre_status != post_status:
                    recombination_distribution[mapping_start + i] += 1
                    if pre_status == 0:
                        recombination_distribution_01[mapping_start + i] += 1
                    else:
                        recombination_distribution_10[mapping_start + i] += 1
                pre_status = post_status

    return recombination_distribution, recombination_distribution_01, recombination_distribution_10, tot_log_lik


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="Makes a prediction on each evidence array and summarises the information in a single recombination array",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--predictions", help="path of the .tsv file containing the prediction arrays")
    parser.add_argument("--hybrid_ref", help="path of the hybrid reference")
    parser.add_argument("--out", help="output path of the .npz file containing the recombination array")
    parser.add_argument("--out_01", help="output path of the .npz file containing the recombination array from 0 to 1")
    parser.add_argument("--out_10", help="output path of the .npz file containing the recombination array from 1 to 0")

    args = parser.parse_args()
    predictions_file = args.predictions
    hybrid_ref_path = args.hybrid_ref
    output_path = args.out
    output_path_01 = args.out_01
    output_path_10 = args.out_10

    recombination_distribution, recombination_distribution_01, recombination_distribution_10, tot_log_lik = (
        get_recombination_array(predictions_file, hybrid_ref_path)
    )

    np.savez(output_path, recombination_distribution)
    np.savez(output_path_01, recombination_distribution_01)
    np.savez(output_path_10, recombination_distribution_10)
