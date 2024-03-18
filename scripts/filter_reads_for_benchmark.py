'''

take the longest 100 reads from a fastq file and write them in a new file.
the parameters are chosen by hand for specific fastq file. it's not a smart script,
but it allows to create quickly a file that contains a subset of reads that will be useful for testing

'''
import gzip
from Bio import SeqIO,bgzf

import numpy as np

#reads="data/pure_reads/EM11_new_chemistry.fastq.gz" #longer than 40 kb
#threshold=40000 #reads longer than 40 kb, they are 111
#out_reads="data/test/test_EM11_new_chemistry.fastq.bgz"

reads="data/pure_reads/EM60_new_chemistry.fastq.gz" #longer than 40 kb
threshold=70000 #reads longer than 70 kb, they are 131
out_reads="data/test/test_EM60_new_chemistry.fastq.bgz"

#reads="data/population_reads/P2_7.fastq.gz"
#threshold=55000 #reads longer than 55 kb, they are 113
#out_reads="data/test/test_P2_7.fastq.bgz"

long_reads=[]

with gzip.open(reads, "rt") as handle:
    for record in SeqIO.parse(handle, "fastq"):
        l=len(record.seq)

        if l>threshold:
            long_reads.append(record)

print(len(long_reads))

with bgzf.BgzfWriter(out_reads, "wb") as outgz:
    SeqIO.write(sequences=long_reads, handle=outgz, format="fastq")
