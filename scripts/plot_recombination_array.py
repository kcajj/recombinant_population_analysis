import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px

if __name__ == "__main__":
    population="P2"
    timestep="7"

    input_path=f"results/genomewide_recombination_arrays/{population}_{timestep}.npz"
    #input_path=f"results/genomewide_recombination_arrays/MAFFT_test_{population}_{timestep}.npz"

    output_path="results/plots/genomewide_recombination/"
    #output_path="results/plots/genomewide_recombination/MAFFT_test_"

    npz=np.load(input_path)
    lst=npz.files
    for item in lst:
        recombination_distribution=npz[item]

    conv_rec_dist = np.convolve(recombination_distribution, np.ones(100), mode='valid')/100
    plt.figure(figsize=(20, 5))
    plt.plot(conv_rec_dist)
    plt.xlabel("position")
    plt.ylabel("recombination rate")
    plt.savefig(f"{output_path}{population}_{timestep}.png")

    fig = px.line(x=range(len(conv_rec_dist)), y=conv_rec_dist)
    fig.write_html(f"{output_path}{population}_{timestep}.html")