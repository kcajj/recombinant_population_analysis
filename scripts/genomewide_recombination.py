import numpy as np
from handle_msa import length_msa
import csv
import sys

def get_recombination_array(predictions_file, refs_msa_path):

    csv.field_size_limit(sys.maxsize)

    recombination_distribution=np.zeros(length_msa(refs_msa_path),dtype=int)
    tot_log_lik=0

    with open(predictions_file) as file:
        tsv_file = csv.reader(file, delimiter="\t")
        for line in tsv_file:
            #read_name=line[0]
            mapping_start=int(line[1])
            #mapping_end=int(line[2])
            log_lik=float(line[3])
            prediction_array_str=line[4][1:-1].split(' ')
            prediction_array=np.array([int(e[:-1]) for e in prediction_array_str])

            tot_log_lik+=log_lik

            pre_status = prediction_array[0]
            for i in range(len(prediction_array)):
                post_status = prediction_array[i]
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
    parser.add_argument("--predictions", help="path of the .tsv file containing the prediction arrays")
    parser.add_argument("--msa_refs", help="path of the msa between the references")
    parser.add_argument("--out", help="output path of the .npz file containing the recombination array")

    args = parser.parse_args()
    predictions_file=args.predictions
    refs_msa_path=args.msa_refs
    output_path=args.out

    recombination_distribution, tot_log_lik = get_recombination_array(predictions_file, refs_msa_path)

    np.savez(output_path,recombination_distribution)