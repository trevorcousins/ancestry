# Files

## msprime_models.py

Defines models in msprime. Each model should have a short descriptor above that describes the history specified. 

## msprime_models_check.py

Checks that models as specified in msprime_models.py are valid (that they can run without error).

# msprime_tms.py

With a given msprime history, it calcuates the transition matrix and plots it. Can use to calculate KL divergence and other things.

# msprime_utils.py

Different functions that are useful for various things regarding an msprime simulation. For example, returning the coalescent data.

# heatmaps_generate.py

Function to generate a heatmap.

# vcf_mhs.py

Functions to write (to disc) a vcf and a muhs (multihetsep) file, given an msprime simulation. 

# msprime_generate.py

Simulates a msprime model, then saves the coalescent data to disc. Is a bit different in the way that the msprime model can be well defined from the command line, making it useful for HPC submissions.

# abinitio_tm.py

Returns a transition matrix of expected probabilities, given a piecewise constant lambda array.
