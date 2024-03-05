import gzip
from Bio import SeqIO,bgzf
import numpy as np
import matplotlib.pyplot as plt
import subprocess

'''

this script gives some statistics on the fastq files of the nanopore run,
it is not important for the recombination analysis

'''

phages=['EM11','EM60']

refs_msa_path="results/msa/refs_msa.fasta"

for phage in phages:
    
    temp_fasta_path = f"results/temp/{phage}_read.fasta"
    temp_output_path = f"results/temp/{phage}_read_msa.fasta"

    c=0
    with gzip.open(f"data/{phage}_new_chemistry.fastq.gz", "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            
            # we need to align the read to the refs_msa.fasta file and extract the evidences

            
            # mafft --auto 
            #       --addfragments results/seq_for_msa/${population}/${timepoint}/${population}_${timepoint}_${read}.fasta
            #       --keeplength results/msa/refs_msa.fasta
            #       > results/msa/${population}/${timepoint}/${population}_${timepoint}_${read}_msa.fasta

            # Create a temporary fasta file with the id and sequence of the read
            with open(temp_fasta_path, "w") as temp_fasta:
                temp_fasta.write(f">{record.id}\n{record.seq}\n")
                temp_fasta.close()

            # create the alignment
            msa_command = f"mafft --auto --addfragments {temp_fasta_path} {refs_msa_path} > {temp_output_path}"
            subprocess.run(msa_command, shell=True)

            c+=1
            if c==2: break

            #remove temporary files
            #rm_command = f"rm {temp_fasta_path} {temp_output_path}"
            #subprocess.run(rm_command, shell=True)
            
            #break   
            
