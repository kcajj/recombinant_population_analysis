# Pipeline

1. The MSA between the two ancestral phage sequences is produced
2. The MSA is used for creating a hybrid reference.
3. The ONT reads from the recombination experiment are mapped on the hybrid reference.
4. Every mapped read is analysed and the evidences of similarity to the two ancestral sequences are extracted.

# Scripts

## Hybrid reference

Extracts the two sequences from the MSA and iterates over each position. If one of the sequences has a gap at a given position, the nucleotide from the other sequence is used. In cases where both sequences have valid nucleotides, one is chosen at random with 50% probability.

## Extract evidence arrays

From the BAM file of the recombinant reads on the hybrid reference the MSA of each read with the two references is reconstructed (handle msa/add_to_msa). From this MSA the distribution of evidences is extracted (handle_msa/extract_evidence_arrays).

## Handle MSA


## Viterbi
