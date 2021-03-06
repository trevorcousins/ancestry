# This script checks that the models as specified in msprime_models.py are functional. 
# example usage:
# $ python msprime_models_check.py m0001 --print --tree

import argparse
from msprime_models import * 
import sys
from vcf_mhs import *
import pdb

parser = argparse.ArgumentParser()
parser.add_argument("model",help="Specify the model to be used, as defined in msprime_models.")
# parser.add_argument("print_",help="Print the DemographicDebugger (from msprime) and the number of segratating sites",nargs='?',default=False,type=bool)
# parser.add_argument("tree",help="Draw trees for this genealogy and info per sesgment?",nargs='?',default=False,type=bool)
parser.add_argument("--print",help="Print the DemographicDebugger (from msprime) and the number of segratating sites",action="store_true")
parser.add_argument("--tree",help="Draw trees for this genealogy and info per segment",action="store_true")
parser.add_argument("--mhs",help="Write a multihetsep file, which can be fed to MSMC.",action="store_true")
parser.add_argument("-suffix","--output_suffix",help="Give a suffix to the mhs and vcf files that will be written",default='',type=str)
args = parser.parse_args()
# TODO: add functionality to change path of vcf and mhs 

#grab args
model = args.model
if args.print:
    print_ = True
else: 
    print_ = False

print('Running simulation of model {}'.format(model))
exec('sim = ' + model + '(' + ')' )
print('Simulation finished')
    
# count number of segregating sites
# not entirely necessary
def seg_sites_info():
    count = 0
    show_variants_info = False # do you want to see the information about each segrating site
    for variant in sim.variants():
        if show_variants_info:
            print("\n\nCount is: ", count)
            print(
                variant.site.id, variant.site.position,
                variant.alleles, variant.genotypes, sep="\t")
        count += 1
    return count

def tree_info(sim):
    for tree in sim.trees():
        print("tree {}: interval = {}: TMRCA {}".format(tree.index, tree.interval, tree.time(tree.parent(0))))
        print(tree.draw(format="unicode"))
    return None

def draw_tree(sim):
    for tree in sim.trees():
        print("tree {}: interval = {}: TMRCA {}".format(tree.index, tree.interval, tree.time(tree.parent(0))))
        print(tree.draw(format="unicode"))
    return None

if print_:
    #seg_sites = seg_sites_info()
    print("Number of segregating sites is: {}".format(sim.num_sites))

if args.tree:
    draw_tree(sim)
    tree_info(sim)

if args.mhs:
    print('mhs is true. Writing vcf')
    # write vcf and mhs
    sim_to_mhs(sim,suffix = '_' + args.output_suffix)

    