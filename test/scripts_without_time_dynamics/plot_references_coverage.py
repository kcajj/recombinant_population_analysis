import numpy as np
from handle_msa import length_msa, extract_references_names
import csv
import sys
import matplotlib.pyplot as plt
import plotly.express as px

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

def plot_references_coverage(array1, array2, output_path, title, x_label, y_label, refs_msa_path):
    plt.figure(figsize=(20, 5))
    plt.plot(array1,alpha=0.5)
    plt.plot(array2,alpha=0.5)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.legend(extract_references_names(refs_msa_path))
    plt.savefig(output_path)

    html_path = output_path[:-3]+"html"
    fig = px.line(x=range(len(array1)), y=array1)
    fig.add_scatter(x=np.array(range(len(array2))), y=array2, mode='lines')
    fig.write_html(html_path)

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="Makes a prediction on each evidence array and summarises the information in a single recombination array",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--predictions", help="path of the .tsv file containing the prediction arrays")
    parser.add_argument("--msa_refs", help="path of the msa between the references")
    parser.add_argument("--out", help="output path of the .npz file containing the coverage plot")

    args = parser.parse_args()
    predictions_file=args.predictions
    refs_msa_path=args.msa_refs
    output_path=args.out

    coverage, coverage0, coverage1 = get_references_coverage(predictions_file, refs_msa_path)

    id = output_path.split("/")[-1].split(".")[0]

    plot_references_coverage(coverage0, coverage1, output_path, f"References coverage {id}", "position", "coverage", refs_msa_path)

    normalised0 = np.divide(coverage0.astype(float), coverage.astype(float), out=np.zeros_like(coverage0.astype(float)), where=coverage.astype(float)!=0)
    normalised1 = np.divide(coverage1.astype(float), coverage.astype(float), out=np.zeros_like(coverage1.astype(float)), where=coverage.astype(float)!=0)

    normalised_output_path = output_path[:-4]+"_normalised.png"

    plot_references_coverage(normalised0, normalised1, normalised_output_path, f"Normalised references coverage {id}", "position", "coverage rate", refs_msa_path)