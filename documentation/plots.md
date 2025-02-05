# Ouput

The main output of the pipeline is a plot that shows the fraction of the population belonging to each of the two reference sequences for each position of the genome. On top of this, a black line, representing the fraction of breakpoints recorded at each position, is present. The plot is created for each timestep at which the population was analyzed.

This is an example:

![unique_plot_example](assets/unique_P1.png)

This plot is the sum of the two plots described below.

## Coverage plot

After gathering the inference carried out on all reads, for each site of the hybrid reference the fraction of reads assigned to ancestral sequence 1 and 2 is plotted.

Example:


## Recombination plot

After gathering the inference carried out on all reads, all the recombination events (i.e. position of a recombinant read where it is inferred the switch from a reference to the other) are plotted on the hybrid reference genome, normalised for the total amount of reads mapped on each position.

Example:

## hyperparameter optimization

The pipeline can also perform the optimization of the hyperparameter θ, corresponding to the recombination frequency in the dataset. The optimization is carried out by maximum likelihood, for details see [here](scripts.md)

![hyperparameter_optimization_example](assets/hyperparameter_optimization.png)