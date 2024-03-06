import gzip
from Bio import SeqIO,bgzf
from Bio.Seq import Seq
import numpy as np
import matplotlib.pyplot as plt
import subprocess
from handle_msa import read_msa, get_evidences_distributions
import time
from collections import defaultdict
import pysam

'''

this script gives some statistics on the fastq files of the nanopore run,
it is not important for the recombination analysis

'''

phages=['EM11','EM60']

refs_msa_path="results/msa/refs_msa.fasta"

null_prob=defaultdict(list)
a_prob=defaultdict(list)
b_prob=defaultdict(list)
time_spent=defaultdict(list)

for phage in phages:
    
    temp_fasta_path = f"results/temp/{phage}_read.fasta"
    temp_output_path = f"results/temp/{phage}_read_msa.fasta"
    
    c = 0
    with pysam.AlignmentFile(f"data/new_chemistry_{phage}.bam", "rb") as bam:
        for read in bam.fetch():
            if not(read.is_secondary) and not(read.is_supplementary):
                if read.is_reverse:
                    read_sequence=Seq(read.query_sequence).reverse_complement()
                else:
                    read_sequence=Seq(read.query_sequence)
                # we need to align the read to the refs_msa.fasta file and extract the evidences
                start_time=time.time()

                # Create a temporary fasta file with the id and sequence of the read
                with open(temp_fasta_path, "w") as temp_fasta:
                    temp_fasta.write(f">{read.query_name}\n{read_sequence}\n")
                    temp_fasta.close()

                # create the alignment
                msa_command = f"mafft --auto --addfragments {temp_fasta_path} {refs_msa_path} > {temp_output_path}"
                subprocess.run(msa_command, shell=True)

                msa_matrix = read_msa(temp_output_path)

                e_distribution = get_evidences_distributions(msa_matrix,i_ref1=0,i_ref2=1,i_extra=2)

                l = len(e_distribution)
                null_prob[phage].append(np.count_nonzero(e_distribution == 0)/l)
                a_prob[phage].append(np.count_nonzero(e_distribution == 1)/l)
                b_prob[phage].append(np.count_nonzero(e_distribution == 2)/l)
                end_time=time.time()
                time_spent[phage].append(end_time-start_time)

                plt.scatter(range(len(e_distribution)), e_distribution, c=e_distribution, marker='|', alpha=0.5)
                plt.title('evidence distribution (purple no evidence, green evidence for A, yellow evidence for B)')

                plt.show()

                #remove temporary files
                rm_command = f"rm {temp_fasta_path} {temp_output_path}"
                subprocess.run(rm_command, shell=True)
                
                c+=1
                if c==10:break

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

print("mean time spent")
for k,v in time_spent.items():
    print(k," ",np.mean(v))
print("")
