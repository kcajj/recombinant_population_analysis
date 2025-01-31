from Bio import AlignIO
import numpy as np
import subprocess
from handle_msa import read_msa, get_evidences_distributions, map_refcoord_msacoord
import time
from collections import defaultdict
import pysam

phages={'EM11':0,'EM60':1}

refs_msa_path="results/msa/refs_msa.fasta"

null_prob=defaultdict(list)
a_prob=defaultdict(list)
b_prob=defaultdict(list)
time_spent=defaultdict(list)

for phage in phages:

    bam_file=f"data/test/test_{phage}_new_chemistry.bam" #for test dataset
    
    map_ref_msa=map_refcoord_msacoord(f"data/references/{phage}_assembly.fasta",refs_msa_path,i_ref_in_msa=phages[phage])
    
    temp_fasta_path = f"test/temp/{phage}_read.fasta"
    temp_output_path = f"test/temp/{phage}_read_msa.fasta"
    temp_refs_msa_path = f"test/temp/{phage}_refs_msa.fasta"
    
    c = 0
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if not(read.is_secondary):
                
                start_time=time.time()

                read_sequence=read.query_sequence

                #find the start and end of mapping. index the ref seuqneces with start and end indexes.
                mapping_start=read.reference_start
                mapping_end=read.reference_end-1

                alignment=AlignIO.read(open(refs_msa_path), "fasta")
                cut_alignment=alignment[:,map_ref_msa[mapping_start]:map_ref_msa[mapping_end]]
                AlignIO.write(cut_alignment, temp_refs_msa_path, "fasta")

                # Create a temporary fasta file with the id and sequence of the read
                with open(temp_fasta_path, "w") as temp_fasta:
                    temp_fasta.write(f">{read.query_name}\n{read_sequence}\n")
                    temp_fasta.close()

                # create the alignment
                msa_command = f"mafft --retree 1 --maxiterate 0 --add {temp_fasta_path} --keeplength {temp_refs_msa_path} > {temp_output_path}"
                subprocess.run(msa_command, shell=True)

                msa_matrix = read_msa(temp_output_path)

                e_distribution_to_plot = get_evidences_distributions(msa_matrix)

                e_distribution = np.where(e_distribution_to_plot > 0, e_distribution_to_plot-1, e_distribution_to_plot)

                l = len(e_distribution)
                null_prob[phage].append(np.count_nonzero(e_distribution == 0)/l)
                a_prob[phage].append(np.count_nonzero(e_distribution == 1)/l)
                b_prob[phage].append(np.count_nonzero(e_distribution == 2)/l)

                #remove temporary files
                rm_command = f"rm {temp_fasta_path} {temp_output_path} {temp_refs_msa_path}"
                subprocess.run(rm_command, shell=True)

                end_time=time.time()
                time_spent[phage].append(end_time-start_time)
                
                c+=1
                print()
                print()
                print()
                print(read.query_name)
                print(mapping_start,mapping_end)
                print(map_ref_msa[mapping_start],map_ref_msa[mapping_end])
                print(c)
                print()
                print()
                print()

print("mean null probability")
for k,v in null_prob.items():
    print(k," ",np.mean(v))
print("")

print("mean a probability")
for k,v in a_prob.items():
    print(k," ",np.mean(v))
print("")

print("mean b probability")
for k,v in b_prob.items():
    print(k," ",np.mean(v))
print("")

print("mean time spent per read")
for k,v in time_spent.items():
    print(k," ",np.mean(v))
print("")
