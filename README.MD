# Files

## msprime_models.py

Defines models in msprime. Each model should have a short descriptor above that describes the history specified. 

##vmsprime_models_check.py

Checks that models as specified in msprime_models.py are valid (that they can run without error).
Usage: python msprime_models_check.py [model] [print] [tree], where model specifies the model, print instructs whether to print info about demographic history and number of segregating sites, and tree instructs whether to draw the coalescent tree with info. 
Example: python msprime_models_check.py m0001 --print --tree
The above shows the tree and the info. Omit --tree to not see tree (similarly with print)

#
This is an edit
