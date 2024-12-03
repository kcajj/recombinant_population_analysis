import matplotlib.pyplot as plt
import numpy as np
from array_compression import npz_extract

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


def plot_coverage_dynamics(
    timesteps,
    recombination_folder,
    coverage_folder,
    hybrid_ref_path,
    references,
    coverage_threshold,
    output_path,
):
    """
    plots the normalized coverage of the two phages and the recombination distribution.

    parameters:
    - timesteps: list of timesteps
    - recombination_folder: path to the folder containing the arrays with the recombination distribution
    - coverage_folder: path to the folder containing the coverage arrays
    - references: name of the references
    - coverage_threshold: coverage threshold below which the positions are masked
    - output_path: path to the output file where the plot will be saved
    """
    figure, subplots = plt.subplots(len(timesteps), 1, figsize=(20, 15), sharex=True, sharey=True)

    for timestep in timesteps:
        coverage_path = f"{coverage_folder}/{timestep}.npz"
        coverage_0_path = f"{coverage_folder}/{timestep}_0.npz"
        coverage_1_path = f"{coverage_folder}/{timestep}_1.npz"

        coverage = npz_extract(coverage_path)
        coverage0 = npz_extract(coverage_0_path)
        coverage1 = npz_extract(coverage_1_path)

        recombination_path = f"{recombination_folder}/{timestep}.npz"

        recombination_distribution = npz_extract(recombination_path)

        normalised_coverage = np.divide(
            coverage0.astype(float),
            coverage.astype(float),
            out=np.zeros_like(coverage0.astype(float)),
            where=coverage.astype(float) != 0,
        )

        normalised_recombination = np.divide(
            recombination_distribution.astype(float),
            coverage.astype(float),
            out=np.zeros_like(recombination_distribution.astype(float)),
            where=coverage.astype(float) != 0,
        )

        # mask out the positions with coverage below the threshold
        mask = [i < coverage_threshold for i in coverage]
        masked_normalised_coverage = np.ma.masked_array(normalised_coverage, mask)
        masked_normalised_recombination = np.ma.masked_array(normalised_recombination, mask)

        # plot the distributions
        figure.tight_layout(pad=0)

        x = np.arange(0, len(coverage))
        subplots[timesteps.index(timestep)].fill_between(x, masked_normalised_coverage)
        subplots[timesteps.index(timestep)].fill_between(x, masked_normalised_coverage, 1)
        subplots[timesteps.index(timestep)].plot(masked_normalised_recombination, alpha=1, linewidth=3, color="black")

        subplots[timesteps.index(timestep)].set_title(f"timestep {timestep}")
        subplots[timesteps.index(timestep)].set_ylabel("normalised \n coverage")
        figure.legend(references, loc="upper right")
        plt.xlabel("genome position")

    plt.savefig(output_path, format="pdf")


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="Plot the recombination array",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--hybrid_ref", help="path of the hybrid reference")
    parser.add_argument("--recombination", help="path to the folder containing the recombination arrays")
    parser.add_argument("--coverage", help="path to the folder containing the coverage arrays")
    parser.add_argument("--timesteps", help="list of timesteps")
    parser.add_argument("--references", help="name of references")
    parser.add_argument("--coverage_threshold", help="coverage threshold")
    parser.add_argument("--out", help="output path of the plots")

    args = parser.parse_args()
    hybrid_ref_path = args.hybrid_ref
    recombination_folder = args.recombination
    coverage_folder = args.coverage
    timesteps = args.timesteps.split(",")[:-1]
    references = args.references.split(",")[:-1]
    coverage_threshold = int(args.coverage_threshold)
    output_path = args.out

    plot_coverage_dynamics(
        timesteps,
        recombination_folder,
        coverage_folder,
        hybrid_ref_path,
        references,
        coverage_threshold,
        output_path,
    )
