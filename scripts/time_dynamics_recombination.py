import matplotlib.pyplot as plt
import numpy as np

SMALL_SIZE = 20
MEDIUM_SIZE = 25
BIGGER_SIZE = 30

plt.rc("font", size=MEDIUM_SIZE)  # controls default text sizes
plt.rc("axes", titlesize=MEDIUM_SIZE)  # fontsize of the axes title
plt.rc("axes", labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc("ytick", labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc("legend", fontsize=MEDIUM_SIZE)  # legend fontsize
plt.rc("figure", titlesize=BIGGER_SIZE)  # fontsize of the figure title


def npz_extract(npz_file):
    npz = np.load(npz_file)
    lst = npz.files
    for item in lst:
        array = npz[item]
    return array


# to convert genome coordinates from hybrid to reference:
"""
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
"""


def plot_recombination_dynamics(
    timesteps, recombination_folder, coverage_folder, references, coverage_threshold, output_path
):
    figure, subplots = plt.subplots(len(timesteps), 1, figsize=(20, 15), sharex=True)

    # if you want to get converted genome coordinates from hybrid to reference
    # bas51_map=get_hyb_ref_map("thesis/alignments/hybrid_on_bas51.bam")

    for timestep in timesteps:
        recombination_path = f"{recombination_folder}/{timestep}.npz"
        recombination_01_path = f"{recombination_folder}/{timestep}_01.npz"
        recombination_10_path = f"{recombination_folder}/{timestep}_10.npz"
        coverage_array_path = f"{coverage_folder}/{timestep}.npz"

        recombination_distribution = npz_extract(recombination_path)
        recombination_distribution_01 = npz_extract(recombination_01_path)
        recombination_distribution_10 = npz_extract(recombination_10_path)

        coverage = npz_extract(coverage_array_path)

        # convert genome coordinates
        """
        converted_recombination = np.zeros_like(recombination_distribution)
        converted_recombination_01 = np.zeros_like(recombination_distribution_01)
        converted_recombination_10 = np.zeros_like(recombination_distribution_10)
        converted_coverage = np.zeros_like(coverage)
        for i in range(len(recombination_distribution)):
            if bas51_map[i]==None: continue
            converted_recombination[bas51_map[i]]+=recombination_distribution[i]
            converted_recombination_01[bas51_map[i]]+=recombination_distribution_01[i]
            converted_recombination_10[bas51_map[i]]+=recombination_distribution_10[i]
            converted_coverage[bas51_map[i]]+=coverage[i]

        normalised_converted = np.divide(converted_recombination.astype(float), converted_coverage.astype(float), out=np.zeros_like(converted_recombination.astype(float)), where=converted_coverage.astype(float)!=0)
        normalised_converted_01 = np.divide(converted_recombination_01.astype(float), converted_coverage.astype(float), out=np.zeros_like(converted_recombination_01.astype(float)), where=converted_coverage.astype(float)!=0)
        normalised_converted_10 = np.divide(converted_recombination_10.astype(float), converted_coverage.astype(float), out=np.zeros_like(converted_recombination_10.astype(float)), where=converted_coverage.astype(float)!=0)
        """
        normalised = np.divide(
            recombination_distribution.astype(float),
            coverage.astype(float),
            out=np.zeros_like(recombination_distribution.astype(float)),
            where=coverage.astype(float) != 0,
        )
        normalised_01 = np.divide(
            recombination_distribution_01.astype(float),
            coverage.astype(float),
            out=np.zeros_like(recombination_distribution_01.astype(float)),
            where=coverage.astype(float) != 0,
        )
        normalised_10 = np.divide(
            recombination_distribution_10.astype(float),
            coverage.astype(float),
            out=np.zeros_like(recombination_distribution_10.astype(float)),
            where=coverage.astype(float) != 0,
        )

        # get the highest values: recombination hotspots
        """
        threshold_frequency=0.01
        with open(f"thesis/plots/hot_sites_{population}.txt","a") as file:
            sorted_indices = np.argsort(normalised_converted)[::-1]
            to_slice=0
            for i in range(len(normalised_converted[sorted_indices])):
                if normalised_converted[sorted_indices[i]]>threshold_frequency: to_slice=i
                else: break
            highest_values = normalised_converted[sorted_indices][:to_slice]
            highest_indices = sorted_indices[:to_slice]
            file.write(f"{population} {timestep}\n")
            file.write(f"Values: {np.around(highest_values,decimals=4)}\n")
            file.write(f"Indices: {highest_indices}\n")
        """

        # mask out the positions with coverage below the threshold
        mask = [i < coverage_threshold for i in coverage]
        masked_normalised_01 = np.ma.masked_array(normalised_01, mask)
        masked_normalised_10 = np.ma.masked_array(normalised_10, mask)

        # plot the distributions
        figure.tight_layout(pad=0)
        subplots[timesteps.index(timestep)].plot(masked_normalised_01, alpha=1, linewidth=3)
        subplots[timesteps.index(timestep)].plot(masked_normalised_10, alpha=1, linewidth=3)
        subplots[timesteps.index(timestep)].set_title(f"timestep {timestep}")
        subplots[timesteps.index(timestep)].set_ylabel("recombination \n frequency")
        figure.legend(
            [f"from {references[0]} to {references[1]}", f"from {references[1]} to {references[0]}"], loc="upper right"
        )
        plt.xlabel("genome position")

    figure.savefig(output_path, format="pdf")


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="Plot the recombination array",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--recombination", help="path to the folder containing the recombination arrays")
    parser.add_argument("--coverage", help="path to the folder containing the coverage arrays")
    parser.add_argument("--timesteps", help="list of timesteps")
    parser.add_argument("--references", help="name of references")
    parser.add_argument("--coverage_threshold", help="coverage threshold")
    parser.add_argument("--out", help="output path of the plots")

    args = parser.parse_args()
    recombination_folder = args.recombination
    coverage_folder = args.coverage
    timesteps = args.timesteps.split(",")[:-1]
    references = args.references.split(",")[:-1]
    coverage_threshold = int(args.coverage_threshold)
    output_path = args.out

    plot_recombination_dynamics(
        timesteps, recombination_folder, coverage_folder, references, coverage_threshold, output_path
    )
