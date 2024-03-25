# SPEEDING UP THE PROCESS

Using MAFFT for each read is too slow. We need to use another approach.

# Hybrid reference approach

0. create a benchmark dataset
1. build a hybrid reference from the MSA
2. align the recombinant reads to the hybrid reference with minimap2
3. remove the insertions of each read from the alignment (no gap in the reference)
4. add the resulting sequence to the msa of the two references, in the mapping region indicated by minimap
5. treat the new msa as the msa built by MAFFT
6. compare the results with the results obtained by MAFFT

## 0. create a test dataset (benchmark and run the MAFFT scripts on it)

1. filter the reads to obtain a subset of them that will be used as common testing dataset: test_EM11_new_chemistry.fastq.bgz, test_EM60_new_chemistry.fastq.bgz, test_P2_7.fastq.bgz

2. align these reads to the respective references and obtain .sam, .bam and .bai files:

<pre>

minimap2 -ax map-ont references/EM11_assembly.fasta test/test_EM11_new_chemistry.fastq.bgz > test/test_EM11_new_chemistry.sam
samtools sort -@ 4 -o test/test_EM11_new_chemistry.bam test/test_EM11_new_chemistry.sam
samtools index test/test_EM11_new_chemistry.bam test/test_EM11_new_chemistry.bam.bai

minimap2 -ax map-ont references/EM60_assembly.fasta test/test_EM60_new_chemistry.fastq.bgz > test/test_EM60_new_chemistry.sam
samtools sort -@ 4 -o test/test_EM60_new_chemistry.bam test/test_EM60_new_chemistry.sam
samtools index test/test_EM60_new_chemistry.bam test/test_EM60_new_chemistry.bam.bai

minimap2 -ax map-ont test/P2_references.fa test/test_P2_7.fastq.bgz > test/test_P2_7.sam
samtools sort -@ 4 -o test/test_P2_7.bam test/test_P2_7.sam
samtools index test/test_P2_7.bam test/test_P2_7.bam.bai

</pre>

parameters estimation with test dataset:

<pre>

mean null probability
EM11   0.967892723956481
EM60   0.9679468350646898

mean a probability
EM11   0.027694027601184522
EM60   0.0008853449723399029

mean b probability
EM11   0.0044132484423344415
EM60   0.031167819962970477

sum
EM11: 1.0
EM60: 1.0000000000000002

mean time spent
EM11   2.401213432613172
EM60   5.300837469427553

</pre>

recombination reads run on test datasett:

<pre>

mean time spent
P2   5.116986224410731

</pre>

BIG BIG BIG BIG BIG PROBLEM: THIS TIME CONSIDERS ALSO THE TIME USED BY THE VITERBI ALGORITHM

## 1. build a hibyrid reference from the MSA

1. build the msa:

    <pre>

    mafft --auto results/msa/refs.fasta > results/msa/refs_msa.fasta

    </pre>

2. build the hybrid reference

    [script](../scripts/hybrid_reference.py)

## 2. align the recombinant reads to the hybrid reference with minimap2

<pre>

minimap2 -ax map-ont results/msa/hybrid_ref.fasta data/test/test_P2_7.fastq.bgz > data/test/hybrid_test_P2_7.sam
samtools sort -@ 4 -o data/test/hybrid_test_P2_7.bam data/test/hybrid_test_P2_7.sam
samtools index data/test/hybrid_test_P2_7.bam data/test/hybrid_test_P2_7.bam.bai

</pre>

## 3. process bam

1. get the alignment, its an array of tuples, each tuple has the index fo the read and the index of the reference. if there is a gap in either one, the value is None.

2. we build the sequence of the aligned read, putting "-" when the alignment is None on the read (gaps) and ignoring the read bases when the alignment is None on the reference (insertions)

3. we cut the msa of the two references in correspondence of the start and end of the mapping of the read and we add the sequence directly to the msa, pasting it under the references

4. we analyse the msa the same as for the other cases

the time required to analyse each read of the test set alignment is:

<pre>

mean time spent
P2   1.8665060363374315

</pre>

that's a good improvement, i wonder how much is the mean time in a random dataset (i tried and it is 0.15)

## 6. compare the results with the results obtained by mafft

We will compare the genome wide recombination:

<pre>
MAFFT:
mean time spent
P2   5.116986224410731
</pre>

img:

![MAFFT_rec](../results/plots/genomewide_recombination/MAFFT_test_P2_7.png)

<pre>
hybrid ref:
mean time spent
P2  1.8665060363374315
</pre>

img:

![rec](../results/plots/genomewide_recombination/test_P2_7.png)

### we can do the same thing for parameter estimation

<pre>

minimap2 -ax map-ont results/msa/hybrid_ref.fasta data/test/test_EM11_new_chemistry.fastq.bgz > data/test/hybrid_test_EM11_new_chemistry.sam
samtools sort -@ 4 -o data/test/hybrid_test_EM11_new_chemistry.bam data/test/hybrid_test_EM11_new_chemistry.sam
samtools index data/test/hybrid_test_EM11_new_chemistry.bam data/test/hybrid_test_EM11_new_chemistry.bam.bai

minimap2 -ax map-ont results/msa/hybrid_ref.fasta data/test/test_EM60_new_chemistry.fastq.bgz > data/test/hybrid_test_EM60_new_chemistry.sam
samtools sort -@ 4 -o data/test/hybrid_test_EM60_new_chemistry.bam data/test/hybrid_test_EM60_new_chemistry.sam
samtools index data/test/hybrid_test_EM60_new_chemistry.bam data/test/hybrid_test_EM60_new_chemistry.bam.bai

</pre>

<pre>
hybrid reference                    MAFFT
mean null probability
EM11   0.9662429815268376           0.967892723956481
EM60   0.9679095862036405           0.9679468350646898

mean a probability
EM11   0.026737786705818084         0.027694027601184522
EM60   0.0011348214186898934        0.0008853449723399029

mean b probability
EM11   0.007019231767344301         0.0044132484423344415
EM60   0.03095559237766953          0.031167819962970477

mean time spent
EM11   0.14966010513590344          2.401213432613172
EM60   0.20331921707205222          5.300837469427553

</pre>

# whole dataset

align all the reads on the hybrid reference.

<pre>

minimap2 -ax map-ont results/msa/hybrid_ref.fasta data/population_reads/P2_7.fastq.gz > results/alignments/P2_7.sam
samtools sort -@ 4 -o results/alignments/P2_7.bam results/alignments/P2_7.sam
samtools index results/alignments/P2_7.bam results/alignments/P2_7.bam.bai

minimap2 -ax map-ont results/msa/hybrid_ref.fasta data/pure_reads/EM11_new_chemistry.fastq.gz > data/test/hybrid_test_EM11_new_chemistry.sam
samtools sort -@ 4 -o dacta/test/hybrid_test_EM11_new_chemistry.bam data/test/hybrid_test_EM11_new_chemistry.sam
samtools index data/test/hybrid_test_EM11_new_chemistry.bam data/test/hybrid_test_EM11_new_chemistry.bam.bai

minimap2 -ax map-ont results/msa/hybrid_ref.fasta data/pure_reads/EM60_new_chemistry.fastq.gz > data/test/hybrid_test_EM60_new_chemistry.sam
samtools sort -@ 4 -o data/test/hybrid_test_EM60_new_chemistry.bam data/test/hybrid_test_EM60_new_chemistry.sam
samtools index data/test/hybrid_test_EM60_new_chemistry.bam data/test/hybrid_test_EM60_new_chemistry.bam.bai

</pre>

we run hmm_rec_reads.py on the whole dataset. i want to see if the result is cool. then we will decide if and how to speed up the process.

i will put a threshold of 5kb of read length. this means that we will analyse 30% of reads and 70% of the total information in terms of mapped bases.

this is the result:

<pre>

mean time spent (per read and per base)
P2   0.42962472400840623
P2   3.931504231066367e-05
total reads 250355
reads used 68314
total time 29375.085064888

</pre>

img:

![rec_whole_dataset](../results/plots/genomewide_recombination/P2_7.png)

