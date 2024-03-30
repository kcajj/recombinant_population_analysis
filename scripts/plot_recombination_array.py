import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import pysam

if __name__ == "__main__":
    '''
    #bam_file = 'results/alignments/P2_7.bam'
    bam_file = 'data/test/hybrid_test_P2_7.bam'
    imported = pysam.AlignmentFile(bam_file, mode = 'rb')  # 'rb' ~ read bam
    coverage_arrays = imported.count_coverage("hybrid_ref")
    coverage=np.zeros(len(coverage_arrays[0]))
    for i in range(len(coverage_arrays[0])):
        for j in range(4):
            coverage[i]+=coverage_arrays[j][i]
    #np.savez("data/test/coverage_tot.npz",coverage)
    np.savez("data/test/coverage_test.npz",coverage)
    '''
    
    #npz=np.load("data/test/coverage_tot.npz")
    npz=np.load("data/test/coverage_test.npz")
    lst=npz.files
    for item in lst:
        coverage=npz[item]
    
    population="P2"
    timestep="7"

    #input_path=f"results/genomewide_recombination_arrays/{population}_{timestep}.npz"
    input_path=f"results/genomewide_recombination_arrays/test_{population}_{timestep}.npz"
    #input_path=f"results/genomewide_recombination_arrays/MAFFT_test_{population}_{timestep}.npz"

    #output_path="results/plots/genomewide_recombination/"
    output_path="results/plots/genomewide_recombination/test_"
    #output_path="results/plots/genomewide_recombination/MAFFT_test_"

    npz=np.load(input_path)
    lst=npz.files
    for item in lst:
        recombination_distribution=npz[item]

    norm_rec_dist = recombination_distribution/coverage
    plt.figure(figsize=(20, 5))
    plt.plot(norm_rec_dist)
    plt.xlabel("position")
    plt.ylabel("recombination rate")
    plt.savefig(f"{output_path}{population}_{timestep}_normalised_for_coverage.png")

    fig = px.line(x=range(len(norm_rec_dist)), y=norm_rec_dist)
    fig.write_html(f"{output_path}{population}_{timestep}_normalised_for_coverage.html")

    plt.figure(figsize=(20, 5))
    plt.plot(recombination_distribution)
    plt.xlabel("position")
    plt.ylabel("recombination rate")
    plt.savefig(f"{output_path}{population}_{timestep}.png")

    fig = px.line(x=range(len(recombination_distribution)), y=recombination_distribution)
    fig.write_html(f"{output_path}{population}_{timestep}.html")