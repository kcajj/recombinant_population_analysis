# SPEEDING UP THE PROCESS

Using MAFFT for each read is too slow. We need to use another approach.

# Hybrid reference approach

1. build a hybrid reference from the MSA
2. align the recombinant reads to the hybrid reference with minimap2
3. remove the insertions of each read from the alignment (no gap in the reference)
4. add the resulting sequence to the msa of the two references, in the mapping region indicated by minimap
5. treat the new msa as the msa built by MAFFT
6. compare the results with the results obtained by MAFFT

