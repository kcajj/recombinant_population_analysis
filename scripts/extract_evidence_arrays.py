import numpy as np
from handle_msa import get_evidences_distributions, add_to_msa
import pysam
import time
import csv
from array_compression import compress_array


def write_evidence_arrays(bam_file, refs_msa_path, output_evidences_path, length_threshold):
    """
    Extracts MSA recombination evidences from alignments of recombinant reads to a hybrid reference.
    Writes the evidence data to a specified output file.

    Parameters:
    - bam_file: Path to the BAM file containing alignments.
    - refs_msa_path: Path to the MSA between the references.
    - output_evidences_path: Path to the output file where evidence arrays will be written.
    - length_threshold: Minimum length of the alignment to be considered.
    """
    time_spent_per_read = []
    tot_time_start = time.time()
    c_tot_alignments = 0
    c_useful_alignments = 0
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        with open(output_evidences_path, "w", newline="") as tsvfile:
            writer = csv.writer(tsvfile, delimiter="\t", lineterminator="\n")
            for read in bam.fetch():
                if not (read.is_secondary):

                    l_alignment = read.query_alignment_length
                    if l_alignment > length_threshold:
                        start_time = time.time()

                        read_seq = read.query_sequence
                        read_msa_seq = ""

                        mapping_start = read.reference_start
                        mapping_end = read.reference_end

                        alignment_array = read.get_aligned_pairs()

                        for read_pos, ref_pos in alignment_array:
                            if ref_pos != None:  # clips and gaps in the reference are not considered
                                if read_pos == None:
                                    read_msa_seq += "-"  # gap in the read
                                else:
                                    read_msa_seq += read_seq[read_pos].lower()

                        msa_matrix = add_to_msa(refs_msa_path, read_msa_seq, mapping_start, mapping_end)

                        e_distribution_to_plot = get_evidences_distributions(msa_matrix)

                        e_distribution = np.where(
                            e_distribution_to_plot > 0, e_distribution_to_plot - 1, e_distribution_to_plot
                        )
                        compressed_e_distribution = compress_array(e_distribution)

                        np.set_printoptions(threshold=np.inf, linewidth=np.inf)
                        writer.writerow([read.query_name, mapping_start, mapping_end, compressed_e_distribution])

                        end_time = time.time()
                        time_spent_per_read.append(end_time - start_time)

                        c_useful_alignments += 1

                c_tot_alignments += 1

    output_stats_path = output_evidences_path[:-4] + "_stats.txt"
    with open(output_stats_path, "w") as f:
        f.write("evidence arrays extraction run of " + bam_file + "\n")
        f.write("mean time spent per read " + str(np.mean(time_spent_per_read)) + "\n")
        f.write("total time " + str(time.time() - tot_time_start) + "\n")
        f.write("total reads " + str(c_tot_alignments) + "\n")
        f.write("reads used " + str(c_useful_alignments) + "\n")


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="Extract all the evidence arrays from the alignments of a bam file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--bam", help="path of the bam file")
    parser.add_argument("--msa_refs", help="path of the msa between the references")
    parser.add_argument("--evidences_out", help="output path of the .tsv file containing the evidence arrays")
    parser.add_argument("--length_threshold", help="minimum length of the alignment to be considered", type=int)

    args = parser.parse_args()
    bam_file = args.bam
    refs_msa_path = args.msa_refs
    output_evidences_path = args.evidences_out
    length_threshold = args.length_threshold

    write_evidence_arrays(bam_file, refs_msa_path, output_evidences_path, length_threshold)
