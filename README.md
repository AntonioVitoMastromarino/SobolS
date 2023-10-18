# OT_based_Sobol_indexes
The main program computes Sobol indexes on a step-wise stochastic simulator accepting two inputs:
  former input (f) is a set of n simulators, namely the steps:
    each simulator has to be a function accepting the output of the previous steps as input together with an optional seed parameter;
  latter input (M) is a list of n integers that specify how many times the step shall run per input;
  an optional input shall include n arrays of legths specified by M containing the seeds to use
    possibly the optional input can be more general;
  the output is a set of 2^n confidence intervals (each represent the variance due to the interaction of a set of uncertainty sources).
The motivating example is a SAIRS model on random network composed by three steps:
  a simple sampler for epistemic parameters that specifies the distribution we are using for these;
  a random network generator like... ;
  a specialized Gillespie... ;
The main program shall support be parallel computing (multiprocessing + OS), consider using Snakemake for this in developing.
