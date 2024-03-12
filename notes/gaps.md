# HOW TO HANDLE GAPS IN HMM PREDICTION OF RECOMBINANT SEQUENCES?

## Skip any kind of gap

![no_gaps](../results/plots/clones/P2_C2_msa.png)

## Consider any kind of gap

![all_gaps](../results/plots/clones/wg_P2_C2_msa.png)

## Skip gaps in references

![no_ref_gaps](../results/plots/clones/cg_P2_C2_msa.png)

## More in detail

- ref1,ref2,read
- x,x,-: no evidence (gap)
- -,-,-: no evidence
- -,-,x: no evidence (insertion)
- x,-,x: evidence 1
- -,x,x: evidence 2
- x,-,-: ?
- -,x,-: ?
