import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import pysam
from Bio import SeqIO
import csv
import sys
from handle_msa import length_msa

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

def get_hyb_ref_map(bam_file):
    map_hyb_ref={}
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if not(read.is_secondary):
                alignment=read.get_aligned_pairs()
                for hyb,ref in alignment:
                    if hyb==None:
                        continue
                    if hyb in map_hyb_ref.keys():
                        if map_hyb_ref[hyb]==None:
                            map_hyb_ref[hyb]=ref
                        else:
                            continue
                    else:
                        map_hyb_ref[hyb]=ref
    return map_hyb_ref

def get_references_coverage(predictions_file, refs_msa_path):

    csv.field_size_limit(sys.maxsize)

    coverage=np.zeros(length_msa(refs_msa_path),dtype=int)
    coverage0=np.zeros(length_msa(refs_msa_path),dtype=int)
    coverage1=np.zeros(length_msa(refs_msa_path),dtype=int)

    with open(predictions_file) as file:
        tsv_file = csv.reader(file, delimiter="\t")
        for line in tsv_file:
            #read_name=line[0]
            mapping_start=int(line[1])
            #mapping_end=int(line[2])
            #log_lik=float(line[3])
            prediction_array_str=line[4][1:-1].split(' ')
            prediction_array=np.array([int(e) for e in prediction_array_str])

            for i in range(len(prediction_array)):
                coverage[mapping_start+i] += 1
                if prediction_array[i] == 0:
                    coverage0[mapping_start+i] += 1
                else:
                    coverage1[mapping_start+i] += 1

    return coverage, coverage0, coverage1

if __name__ == "__main__":
    populations=["P2"]
    timesteps=["1","3","5","7"]

    recombination_folder="results/genomewide_recombination"
    coverage_folder="results/coverage_arrays"
    refs_msa_path="results/msa/msa_refs.fasta"

    bas51_map=get_hyb_ref_map("thesis/alignments/hybrid_on_bas51.bam")

    for population in populations:
        figure, subplots = plt.subplots(len(timesteps), 1, figsize=(20, 15), sharex=True, sharey=True)
        output_path=f"thesis/plots/time_dynamics_coverage_{population}.png"

        for timestep in timesteps:
            predictions_file=f"results/prediction_arrays/{population}/{population}_{timestep}.tsv"

            coverage, coverage0, coverage1 = get_references_coverage(predictions_file, refs_msa_path)

            #convert genome coordinates
            converted_coverage = np.zeros_like(coverage)
            converted_coverage0 = np.zeros_like(coverage0)
            converted_coverage1 = np.zeros_like(coverage1)
            for i in range(len(coverage)):
                if bas51_map[i]==None: continue
                converted_coverage[bas51_map[i]]+=coverage[i]
                converted_coverage0[bas51_map[i]]+=coverage0[i]
                converted_coverage1[bas51_map[i]]+=coverage1[i]

            normalised_converted1 = np.divide(converted_coverage0.astype(float), converted_coverage.astype(float), out=np.zeros_like(converted_coverage0.astype(float)), where=converted_coverage.astype(float)!=0)
            normalised_converted0 = np.divide(converted_coverage1.astype(float), converted_coverage.astype(float), out=np.zeros_like(converted_coverage1.astype(float)), where=converted_coverage.astype(float)!=0)

            #normalised_converted0 = np.divide(coverage0.astype(float), coverage.astype(float), out=np.zeros_like(coverage0.astype(float)), where=coverage.astype(float)!=0)
            #normalised_converted1 = np.divide(coverage1.astype(float), coverage.astype(float), out=np.zeros_like(coverage1.astype(float)), where=coverage.astype(float)!=0)

            # plot the distributions
            figure.tight_layout(pad=0)
            subplots[timesteps.index(timestep)].plot(normalised_converted0, alpha=0.5, linewidth=3)
            subplots[timesteps.index(timestep)].plot(normalised_converted1, alpha=0.5, linewidth=3)
            subplots[timesteps.index(timestep)].set_title(f"timestep {timestep}")
            subplots[timesteps.index(timestep)].set_ylabel("covererage")
            figure.legend(["bas51","bas54"], loc='upper right')
            plt.xlabel("genome position")

        plt.savefig(output_path)