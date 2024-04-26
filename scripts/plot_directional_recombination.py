import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px

def npz_extract(npz_file):
    npz=np.load(npz_file)
    lst=npz.files
    for item in lst:
        array=npz[item]
    return array

def plot(array1, array2, output_path, title, x_label, y_label):
    plt.figure(figsize=(20, 5))
    plt.plot(array1, alpha=0.5)
    plt.plot(array2, alpha=0.5)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.savefig(output_path)

    html_path = output_path[:-3]+"html"
    fig = px.line(x=range(len(array1)), y=array1)
    fig.add_scatter(x=np.array(range(len(array2))), y=array2, mode='lines')
    fig.write_html(html_path)


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="Plot the recombination array",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--recombination_01", help="path of the .npz recombination array from 0 to 1")
    parser.add_argument("--recombination_10", help="path of the .npz recombination array from 1 to 0")
    parser.add_argument("--coverage", help="path of the .npz coverage array")
    parser.add_argument("--out", help="output path of the plots")

    args = parser.parse_args()
    recombination_01_path=args.recombination_01
    recombination_10_path=args.recombination_10
    coverage_array_path=args.coverage
    output_path=args.out

    recombination_distribution_01 = npz_extract(recombination_01_path)
    recombination_distribution_10 = npz_extract(recombination_10_path)

    coverage = npz_extract(coverage_array_path)

    id = output_path.split("/")[-1].split(".")[0]

    plot(recombination_distribution_01, recombination_distribution_10, output_path, f"Recombination distribution {id}", "position", "recombination rate")
    
    normalised_01 = np.divide(recombination_distribution_01.astype(float), coverage.astype(float), out=np.zeros_like(recombination_distribution_01.astype(float)), where=coverage.astype(float)!=0)
    normalised_10 = np.divide(recombination_distribution_10.astype(float), coverage.astype(float), out=np.zeros_like(recombination_distribution_10.astype(float)), where=coverage.astype(float)!=0)

    normalised_output_path = output_path[:-4]+"_normalised.png"
    
    plot(normalised_01, normalised_10, normalised_output_path, f"Normalised recombination distribution {id}", "position", "recombination rate")