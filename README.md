# recombinant_population_analysis

This repository represents the largest portion of my [Bachelor thesis](/thesis/Thesis_Giacomo_Castagnetti.pdf) in Genomics at University of Bologna. This project started in summer 2023 during the Biozentrum research summer project, at [NeherLab](https://neherlab.org/). I graduated in July 2024, with a final score of 110/110 cum laude, along with a honourable mention.

This repository processes the data produced by an experiment performed with the Aionostat, a machine that allows automatic experimental evolution of phages.

The experiment consists in evolving 3 phages with a bacterium that they do not infect well initially. Through an initial [exploratory analysis](https://github.com/kcajj/rec_genome_analysis) that I performed, we know that two phages recombined and one got extinct. The objective of this repository is to create a more precise and high througput procedure to measure recombination in population sequencing data.

We produced a pipeline that carries out the whole analysis, yielding the plots of the coverage of the genomes of the two phages for each timestep of the experiment.

test