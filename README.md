# Recombinant Population Analysis

This repository represents the largest portion of my [Bachelor thesis](/thesis/Thesis_Giacomo_Castagnetti.pdf). I graduated in Genomics at University of Bologna. This project started in summer 2023 during the Biozentrum research summer project, at [NeherLab](https://neherlab.org/). I graduated in July 2024, with a final score of 110/110 cum laude, along with a honourable mention.

The main use of this repository processes the data produced by an [Aionostat](https://edoc.unibas.ch/96360/) experiment, a machine that allows automatic experimental evolution of phages. You can have details on the specific experiment discussed in the thesis by exploring its conents. Right now, this repository aims to be a general tool that can be applied to analyse any heterogeneous population of recombinant molecular entities sequenced with ONT. The experimental requirements are the following:

- The recombinant population has to arise from just two ancestral species that mixed.
- The two ancestral species have to be fairly similar, allowing homologous recombination.
- Only homologous recombination is detected

The pipeline follows the following schematic workflow:



The two references corresponding to the ancestral phages are combined in a hybrid reference. This reference can be used to align the reads of the recombinant population with minimap2. For each read, the obtained alignment is approximated as being the MSA of the 2 references + the recombinant read. From the MSA the evidences of the read belonging to ancestral sequence 1 or 2 are extracted and feeded to the HMM model. To have more details on the HMM model see [here](documentation/hmm.md)

# Configuration

You can run the pipeline by properly setting up the [run_config.yml](/run_config.yml) file and by creating a folder with the input data.

## Input folder

The input folder should have the following structure:

data/
    reads/
        [replicate_code]_[timestep_code].fastq.gz
        ...
    references.fasta

Each fastq file should be named with two codes, one identifying the experimental replicate and one progressively numbering successive timestep (in case of a time series analysis).

The two reference genomes should be included in the same fasta file named "references.fasta".

## run_config.yml

The run_config file has 4 sections:

- run_config: describes the file configuration of the pipeline run. Write down the name of the two references and of the replicates and timesteps that have to be analyzed.

- alignments: set the length threshold below which the reads will be ignored and not aligned to the hybrid reference.

- HMM: define the HMM parameters. To have more details see [here](documentation/hmm.md)

- plots: set the coverage threshold below which the inferences carried out on the site will not be shown in the plots.

# Running the pipeline

## Local execution

<pre>
snakemake --profile local --configfile run_config.yml
</pre>

## HPC execution 

<pre>
snakemake --profile cluster --configfile run_config.yml
</pre>

# Ouput

Two plots are produced by the pipeline:

## Coverage plot

After gathering the inference carried out on all reads, for each site of the hybrid reference the fraction of reads assigned to ancestral sequence 1 and 2 is plotted.

Example:


## Recombination plot

After gathering the inference carried out on all reads, all the recombination events (i.e. position of a recombinant read where it is inferred the switch from a reference to the other) are plotted on the hybrid reference genome, normalised for the total amount of reads mapped on each position.

Example:
