# SPEEDING UP THE PROCESS

Using MAFFT for each read is too slow. We need to use another approach.

# Hybrid reference approach

1. build a hybrid reference from the MSA
2. align the recombinant reads to the hybrid reference with minimap2
3. remove the insertions of each read from the alignment (no gap in the reference)
4. add the resulting sequence to the msa of the two references, in the mapping region indicated by minimap
5. treat the new msa as the msa built by MAFFT
6. compare the results with the results obtained by MAFFT

## 1. build a hibyrid reference from the MSA

1. build the msa:

    <pre>

    </pre>

2. build the hybrid reference:

## 6. compare the results with the results obtained by mafft

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