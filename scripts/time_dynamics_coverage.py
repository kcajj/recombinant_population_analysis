import numpy as np
import matplotlib.pyplot as plt
import csv
import sys
from handle_msa import length_seq
from array_compression import decompress_array, retrive_compressed_array_from_str

SMALL_SIZE = 20
MEDIUM_SIZE = 25
BIGGER_SIZE = 30

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

def get_references_coverage(predictions_file, hybrid_ref_path):

    csv.field_size_limit(sys.maxsize)

    coverage=np.zeros(length_seq(hybrid_ref_path),dtype=int)
    coverage0=np.zeros(length_seq(hybrid_ref_path),dtype=int)
    coverage1=np.zeros(length_seq(hybrid_ref_path),dtype=int)

    with open(predictions_file) as file:
        tsv_file = csv.reader(file, delimiter="\t")
        for line in tsv_file:
            #read_name=line[0]
            mapping_start=int(line[1])
            #mapping_end=int(line[2])
            #log_lik=float(line[3])
            
            compressed_prediction_array=retrive_compressed_array_from_str(line[4])
            prediction_array=decompress_array(compressed_prediction_array)

            for i in range(len(prediction_array)):
                coverage[mapping_start+i] += 1
                if prediction_array[i] == 0:
                    coverage0[mapping_start+i] += 1
                else:
                    coverage1[mapping_start+i] += 1

    return coverage, coverage0, coverage1

def plot_coverage_dynamics(timesteps,prediction_folder,hybrid_ref_path,references,coverage_threshold,output_path):
    figure, subplots = plt.subplots(len(timesteps), 1, figsize=(20, 15), sharex=True, sharey=True)

    for timestep in timesteps:
        predictions_file=f"{prediction_folder}/{timestep}.tsv"

        coverage, coverage0, coverage1 = get_references_coverage(predictions_file, hybrid_ref_path)

        normalised0 = np.divide(coverage0.astype(float), coverage.astype(float), out=np.zeros_like(coverage0.astype(float)), where=coverage.astype(float)!=0)
        normalised1 = np.divide(coverage1.astype(float), coverage.astype(float), out=np.zeros_like(coverage1.astype(float)), where=coverage.astype(float)!=0)
        
        # mask out the positions with coverage below the threshold
        mask = [i < coverage_threshold for i in coverage]
        masked_normalised0 = np.ma.masked_array(normalised0, mask)
        masked_normalised1 = np.ma.masked_array(normalised1, mask)

        # plot the distributions
        figure.tight_layout(pad=0)
        subplots[timesteps.index(timestep)].plot(masked_normalised0, alpha=1, linewidth=3)
        subplots[timesteps.index(timestep)].plot(masked_normalised1, alpha=1, linewidth=3)
        subplots[timesteps.index(timestep)].set_title(f"timestep {timestep}")
        subplots[timesteps.index(timestep)].set_ylabel("normalised \n coverage")
        figure.legend(references, loc='upper right')
        plt.xlabel("genome position")

    plt.savefig(output_path,dpi=600)

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="Plot the recombination array",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--hybrid_ref", help="path of the hybrid reference")
    parser.add_argument("--prediction", help="path of the folder containing prediction arrays")
    parser.add_argument("--timesteps", help="list of timesteps")
    parser.add_argument("--references", help="name of references")
    parser.add_argument("--coverage_threshold", help="coverage threshold")
    parser.add_argument("--out", help="output path of the plots")

    args = parser.parse_args()
    hybrid_ref_path=args.hybrid_ref
    prediction_folder=args.prediction
    timesteps=args.timesteps.split(",")[:-1]
    references=args.references.split(",")[:-1]
    coverage_threshold=int(args.coverage_threshold)
    output_path=args.out
    
    plot_coverage_dynamics(timesteps,prediction_folder,hybrid_ref_path,references,coverage_threshold,output_path)