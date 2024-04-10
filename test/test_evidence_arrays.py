import numpy as np
from handle_msa import get_evidences_distributions, add_to_msa
import pysam
import time
import csv

if __name__ == "__main__":

    population='P2'
    timestep='7'

    refs_msa_path="results/msa/refs_msa.fasta"

    bam_file=f"data/test/hybrid_test_{population}_{timestep}.bam" #for test dataset

    output_path=f"test/results/evidence_arrays/test_{population}_{timestep}.tsv" #for test dataset

    c_tot_alignments=0
    c_useful_alignments=0
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        with open(output_path, 'w', newline='') as tsvfile:
            writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
            for read in bam.fetch():
                if not(read.is_secondary):
                    read_seq = read.query_sequence
                    read_msa_seq = ''

                    l_read=read.infer_read_length()
                    l_alignment=read.query_alignment_length

                    mapping_start=read.reference_start
                    mapping_end=read.reference_end

                    start_time=time.time()

                    alignment_array = read.get_aligned_pairs()

                    for (read_pos,ref_pos) in alignment_array:
                        if ref_pos!=None: #clips and gaps ar
                            if read_pos==None:
                                read_msa_seq+='-'
                            else:
                                read_msa_seq+=read_seq[read_pos].lower()

                    msa_matrix = add_to_msa(refs_msa_path, read_msa_seq, mapping_start, mapping_end)

                    e_distribution_to_plot = get_evidences_distributions(msa_matrix)

                    e_distribution = np.where(e_distribution_to_plot > 0, e_distribution_to_plot-1, e_distribution_to_plot)
                    
                    np.set_printoptions(threshold=np.inf)
                    writer.writerow([mapping_start, mapping_end, e_distribution])

                    c_useful_alignments+=1
                    print(c_tot_alignments)

                c_tot_alignments+=1

    print("total reads", c_tot_alignments)
    print("reads used", c_useful_alignments)