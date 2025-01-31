# FINDING RECOMBINANT READS

The problem we are facing and the strategy to solve it are well explained in [this pdf](plan.pdf).

For the following notes I will use the following notation:
- A: first phage genome
- B: second phage genome
- a: evidence on the read for the 

## 1. Create a functioning HMM script

The first step is to create a working HMM that can get in input the type of data we are dealing with (an array of recombination evidences, a/./b), and that provides in output what we expect (the series of hidden states corresponing to the genome that provided the particular evidence, A/B).

I will try to use the following implementation of the viterbi algorithm: https://medium.com/@zhe.feng0018/coding-viterbi-algorithm-for-hmm-from-scratch-ca59c9203964

### generate data on the basis of given probabilities

To test that it is working i have to define some reasonable probability matrices. we need three matrices: initial probability matrix, transition probability matrix and emission probability matrix.

- initial probability matrix

|    |      |
|----|------|
|A   |0.5   |
|B   |0.5   |

- transition probability matrix

|    |A       |B       |
|----|--------|--------|
|A   |0.99    |0.01    |
|B   |0.01    |0.99    |

- emission probability matrix

|    |.        |a        |b        |
|----|---------|---------|---------|
|A   |0.59     |0.40     |0.01     |
|B   |0.59     |0.01     |0.40     |

Now i will generate sequences with these probabilities

- hidden states: "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"

- emissions:     "bb...b......bb.b.b.b...bb.bb....b.baa...a.a..aaabaaaa.aaa....aa.aa....a.aa..a......a.a.aaa.........aa"

both the results seem reasonable.

Of course in the case of the experiment the amount of a and b will be a lot less.

### make predictions on the generated data

i think that to use the viterbi script i have to convert every title of the matrices to integers (indexes)

- initial probability matrix

|    |      |
|----|------|
|0   |0.5   |
|1   |0.5   |

- transition probability matrix

|    |0       |1       |
|----|--------|--------|
|0   |0.99    |0.01    |
|1   |0.01    |0.99    |

- emission probability matrix

|    |0        |1        |2        |
|----|---------|---------|---------|
|0   |0.59     |0.40     |0.01     |
|1   |0.59     |0.01     |0.40     |

for this reason the previously produced sequence has to be translated to numbers

this is the array of emitted states:

[1, 1, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 1, 1, 2, 1, 2, 1, 2, 1, 2, 2, 2, 1, 1, 2, 1, 1, 2, 2, 2, 2, 1, 2, 1, 0, 0, 2, 2, 2, 0, 2, 0, 2, 2, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 2, 2, 2, 2, 0, 0, 2, 0, 0, 2, 2, 2, 2, 0, 2, 0, 0, 2, 2, 0, 2, 2, 2, 2, 2, 2, 0, 2, 0, 2, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0]

this is the array of hidden states:

[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

this is the predicted array of hidden states:

[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

the prediction is 100% correct. I will run the script multiple times with random sequences

the accuracy is always around 99%. for now we are happy with this result.

## 2. Test on Aionostat data

### phage isolates sequences

The recombination evidences are extracted from the msa of the two reference sequences and the assembly of the phage clone. This is done through the handle_msa.py script.

In hmm.py it is defined the Viterbi algorithm and the probability matrices. We take the recombination evidences of the clones and we give them in input to the Viterbi algorithm.

for the clones i used the following matrices:

- initial probability matrix

|    |      |
|----|------|
|0   |0.5   |
|1   |0.5   |

- transition probability matrix

|    |0       |1       |
|----|--------|--------|
|0   |0.999   |0.001   |
|1   |0.001   |0.999   |

- emission probability matrix

|    |0         |1        |2        |
|----|----------|---------|---------|
|0   |0.949     |0.05     |0.001    |
|1   |0.949     |0.001    |0.05     |

to obtain sensible results from i modified the viterbi algorithm to do additions between logarithms instead of multiplications between probabilities. this is because otherwise the probability would become too small.

this is the result:
![example_clone_prediction](../results/plots/clones/P2_C2_msa.png)

WHAT TO DO WITH GAPS? we try to give an answer [here](gaps.md)

### test on reads of phage population

the same procedure is applied to the reads of the population. in this case the array of evidences is much more noisy:

![noisy_read_evidences](images/noisy_read_evidences.png)

i tried to tweak a bit the parameters but the resutls have many problems:

![noisy_read_prediction](images/noisy_read_prediction.png)

these problems were just due to a bug in the code, the actual situation is much nicer!!!

![nice_read_evidences](images/nice_read_evidences.png)

at this point it would be nice to infer the correct parameters of the probability matrices.

## 3. Parameters estimation

### emission probabilities of reads

to estimate the emission probabilities we can use reads of a sequencing run of a pure phage. by aligning these reads to the msa of references we can see with which frequency ., a and b evidences occur. we will just count the occurrence of each visible state.

The creation of the msa between the two references and a read takes more or less 30 seconds.

1. TIME PROBLEM

    - we could run the same script in parallel for multiple reads.

    - we could align multiple reads at the same time.

    WAIT! the huge amount of time is taken only by the reads that map badly in the msa. i.e. reads that map for the whole length of the references instead of in a specific region. if we fix this problem probably also the msa production time will be reduced a lot.

    the problem could be:

    1. reads that map badly cannot produce a good msa: probably wrong, the first read of EM11 maps really well with minimap2 but not with mafft
    2. maybe the reads are reversed or complemented: I will explore the .sam mapping of the first reads of the .fastq file to see if there is a correlation between a bad msa and some characteristic of the mapping of the read on the reference genome.
        turns out this is true, i already dealt with this problem in rec_genome_analysis and i took the reverse complement of the reads that were mapped in the in the other direction with minimap2. maybe i can try to use some options of mafft before implementing all of this again.
        
        i will try to use the --adjustdirection option, let's hope it is enough. it doesn't work, i will try with --adjustdirectionaccurately. it works!!! still it's a bit slow. I will try to run 100 reads for each phage to have an idea of the time needed.

        these are the results:
        <pre>
        mean time spent
        EM11   7.477514550685883
        EM60   10.225129489898682
        </pre>
        we can see if by giving in input directly the complemented reads and not using adjustidrectionaccurately we obtain a better result in terms of time. i will do this by using read.query_sequence, according to pysam documentation it gives the read sequence already complemented if the read is binding reversely.
        <pre>
        mean time spent
        EM11   2.8776113176345826
        EM60   6.12469485282898
        </pre>
        This time is really not acceptable.
2. BORDERS OF THE READ PROBLEM

    right now the script takes the frequencies from the whole length of the genome. maybe we can limit somehow the region in which the read is mapping.

In order to solve these problems we will try another approach: 

CUTTING THE MSA AROUND THE MAPPING REGION OF THE READ INDIVIDUATED WITH MINIMAP

since we are iterating through the bam file of the alignment between the reads of the sequencing run for a pure phage and the genome of the phage, we already have alignment information.

We can take the start and end positions for the alignment found by minimap2 and use these indexes to cut the MSA between the two references to have a shorter one. We will then add to the shorter msa the read, this will for sure speed up the process.

there are two main problems in this approach:

1. We have the mapping information for the read on just one reference, can we use this information to index both references. we will discuss about this problem in [this notebook](../notebooks/refs_plots.ipynb).

2. We have to create a map from the index of the plain references to their msa. to do this i created the [map_dictionary](scripts/handle_msa.py) function.

Now that we have solved these problems we can apply the approach. I will run the script on 100 reads taken more or less randomly along the genome.

<pre>

mean null probability
EM11   0.9825947907071777
EM60   0.9772102791664006

mean a probability
EM11   0.008471984172513073
EM60   0.0009509440754558615

mean b probability
EM11   0.008933225120309123
EM60   0.021838776758143563

mean time spent
EM11   0.5517718839645386
EM60   0.3548362874984741

</pre>

The time is becoming acceptable. We obtained this time by using the fastest msa method of MAFFT

The occurrences of evidences are weirdly different betweeen different phages: there is a mistake.

THE PROBLEM IS THAT I AM WORKING WITH THE PURE PHAGE SEQUENCING RUN ALIGNMENT MADE ON THE FORWARD SEQUENCE OF THE PHAGE AND WITH THE RECOMBINANT ANALYSIS + MSA MADE ON THE REVERSE COMPLEMENT.

I will realign the reads:

<pre>
minimap2 -ax map-ont references/EM11_assembly.fasta EM11_new_chemistry.fastq.gz > new_chemistry_EM11.sam
samtools sort -@ 4 -o new_chemistry_EM11.bam new_chemistry_EM11.sam
samtools index new_chemistry_EM11.bam new_chemistry_EM11.bam.bai
</pre>

now the results are the follwing:

mean null probability
EM11   0.9747077404528025
EM60   0.9772102791664006

mean a probability
EM11   0.024958402540660604
EM60   0.0009509440754558615

mean b probability
EM11   0.0003338570065367876
EM60   0.021838776758143563

sum
EM11: 0.9999999999999999
EM60: 1.0

mean time spent
EM11   0.38000538110733034
EM60   0.4888264489173889

much nicer!! now we can use these parameters for the model!

## 4. prediction on reads of phage population

now that we have an idea of the parameters of our model, we can start to set up the prediction framework.

We will use the same approach as for parameter estimation.

in this case we go through the reads mapped to both references, for each read:
- we check that it is not a secondary or supplementary alignment
- we check to which reference it is aligned
- we extract the start and end of alignmeeent
- we cut the msa on the basis of the index_map between the phage reference of the read and the msa of both references
- we add the read to the msa. we give to mafft the whole read sequence (even if it maps less) and we ask mafft to do the msa by keeping the length of the two references. (DO WE HAVE TO CHANGE THIS?)
- we do the prediction on msa of the 3 sequences
- we plot the result

The analysis takes roughly: 0.7569327437600424 seconds per read