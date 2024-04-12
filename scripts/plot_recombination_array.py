import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px

def npz_extract(npz_file):
    npz=np.load(npz_file)
    lst=npz.files
    for item in lst:
        array=npz[item]
    return array

def plot(array, output_path, title, x_label, y_label):
    plt.figure(figsize=(20, 5))
    plt.plot(array)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.savefig(output_path)

    html_path = output_path[:-3]+"html"
    fig = px.line(x=range(len(array)), y=array)
    fig.write_html(html_path)

if __name__ == "__main__":
    
    population="P2"
    timestep="7"

    recombination_distribution_path=f"results/genomewide_recombination_arrays/{population}_{timestep}.npz"

    coverage_array_path=f"results/coverage_arrays/{population}_{timestep}.npz"

    output_path=f"results/plots/genomewide_recombination/{population}_{timestep}"

    import argparse

    parser = argparse.ArgumentParser(
        description="Plot the recombination array",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--recombination", help="path of the .npz recombination array")
    parser.add_argument("--coverage", help="path of the .npz coverage array")
    parser.add_argument("--out", help="output path of the plots")

    args = parser.parse_args()
    recombination_distribution_path=args.recombination
    coverage_array_path=args.coverage
    output_path=args.out
    
    recombination_distribution = npz_extract(recombination_distribution_path)

    coverage = npz_extract(coverage_array_path)

    id = output_path.split("/")[-1].split(".")[0]

    plot(recombination_distribution, output_path, f"Recombination distribution {id}", "position", "recombination rate")

    normalised = recombination_distribution/coverage
    normalised_output_path = output_path[:-4]+"_normalised.png"
    
    plot(normalised, normalised_output_path, f"Normalised recombination distribution {id}", "position", "recombination rate")