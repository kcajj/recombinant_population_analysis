# SPEEDING UP THE PROCESS

Using MAFFT for each read is too slow. We need to use another approach.

# Hybrid reference approach

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
EM11   0.968531100625019
EM60   0.9675525481255135

mean a probability
EM11   0.02994434031592349
EM60   0.0004993323494807642

mean b probability
EM11   0.0015245590590574375
EM60   0.03194811952500564

sum
EM11: 1.0
EM60: 0.9999999999999999

mean time spent
EM11   3.396820112105903
EM60   5.120701194719504

</pre>

recombination reads run on test datasett:

<pre>

mean time spent
P2   14.039401678102356

</pre>

BIG BIG BIG BIG BIG PROBLEM: THIS TIME CONSIDERS ALSO THE TIME USED BY THE VITERBI ALGORITHM

the actual time required for MAFFT is 5.748953744769096


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

<pre>

mean time spent
P2   5.394195290548461

</pre>

that's a good improvement, i wonder how much is the mean time in a random dataset

## 6. compare the results with the results obtained by mafft


# whole dataset

minimap2 -ax map-ont results/msa/hybrid_ref.fasta data/population_reads/P2_7.fastq.gz > results/alignments/P2_7.sam
samtools sort -@ 4 -o results/alignments/P2_7.bam results/alignments/P2_7.sam
samtools index results/alignments/P2_7.bam results/alignments/P2_7.bam.bai