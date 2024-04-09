import numpy as np
from handle_msa import get_evidences_distributions, add_to_msa
import pysam
import time
import csv

if __name__ == "__main__":

    population='P2'
    timestep='7'

    refs_msa_path="results/msa/refs_msa.fasta"

    bam_file=f"results/alignments/{population}_{timestep}.bam" #for test dataset

    output_path=f"results/evidence_arrays/{population}_{timestep}.tsv" #for test dataset

    length_treshold=5000

    time_spent_per_read=[]
    tot_time_start=time.time()
    c_tot_alignments=0
    c_useful_alignments=0
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        with open(output_path, 'w', newline='') as tsvfile:
            writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
            for read in bam.fetch():
                if not(read.is_secondary):
                    
                    l_alignment=read.query_alignment_length
                    if l_alignment>length_treshold:
                        start_time=time.time()

                        read_seq = read.query_sequence
                        read_msa_seq = ''

                        mapping_start=read.reference_start
                        mapping_end=read.reference_end

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

                        end_time=time.time()
                        time_spent_per_read.append(end_time-start_time)

                        c_useful_alignments+=1
                        print(c_tot_alignments)

                c_tot_alignments+=1

    print("mean time spent per read", np.mean(time_spent_per_read))
    print("total time", time.time()-tot_time_start)
    print("total reads", c_tot_alignments)
    print("reads used", c_useful_alignments)