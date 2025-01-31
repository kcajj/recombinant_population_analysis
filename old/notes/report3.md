# General version of the pipeline

I am still working on the pipeline because we will publish a paper

## Modifications

- creation of a proper run_config file that can accept various configurations and experimental setting. I will use as a template the run_config file of evo-genome-analysis pipeline from Marco Molari.

- implementation of time dynamics plots directly in the pipeline. this can be done by splitting the pipeline.

- implementation of a coverage threshold on the plots, we don't want to plot the normalised HMM prediction if there are few reads.

- implementation of a system to compress and decompress arrays. this will make the pipeline more space efficient.

- prettyfying the pipeline to make it usable from more people, hiding thesis-specific stuff and exploratory scripts.

## Unique plot

We want to add a combination of the two plots to put it in the paper.

1. we plot only one coverage line and color above and below it

2. we plot the recombination frequency on top of the coverage