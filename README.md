# Recombinant Population Analysis

This repository represents the largest portion of my [Bachelor thesis](/thesis/Thesis_Giacomo_Castagnetti.pdf). I graduated in Genomics at University of Bologna. This project started in summer 2023 during the Biozentrum research summer project, at [NeherLab](https://neherlab.org/). I graduated in July 2024, with a final score of 110/110 cum laude, along with a honourable mention.

The main use of this repository processes the data produced by an [Aionostat](https://edoc.unibas.ch/96360/) experiment, a machine that allows automatic experimental evolution of phages. You can have details on the specific experiment discussed in the thesis by exploring its conents. Right now, this repository aims to be a general tool that can be applied to analyse any heterogeneous population of recombinant molecular entities sequenced with ONT. The experimental requirements are the following:

- The recombinant population has to arise from just two ancestral species that mixed.
- The two ancestral species have to be fairly similar, allowing homologous recombination.
- Only homologous recombination is detected

# Configuration

You can run the pipeline by properly setting up the [run_config.yml](/run_config.yml) file and by creating a folder with the input data.

instructions...

# Ouput

We produced a pipeline that carries out the whole analysis, yielding the plots of the coverage of the genomes of the two phages for each timestep of the experiment.