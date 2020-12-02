# create msprime simulations and then save mhs and vcfs for there

import os
from datetime import datetime
import argparse
from msprime_models import * 
# from abinitio_tm import scaled_time_intervals
# from msprime_utils import scaled_time_intervals, get_het, round_coal_times, tm_counts, get_coal_data, round_bin_coal_data, normalise
from msprime_utils import *
from vcf_mhs import *
import numpy as np
import pdb
import math
import pandas as pd
from scipy.stats import entropy

# parse arguments
#  instant_struct0001(N_0,mig_prop,t_1,t_2,seq_length,recomb_rate,mut_rate,print)
parser = argparse.ArgumentParser()
parser.add_argument("model",help="Specify the model to be used, as defined in msprime_models.")
parser.add_argument("-N","--N_0",help="Initial population size (Default is 1e+04",default=1e+04,type=int)
parser.add_argument("-L","--seq_length",help="Length of the sequence to be simulated (default 3e+07)",default=int(3e+07),type=int)
parser.add_argument("-mig","--migration_prop",help="Proportion of migrating population (default 0.3)",default=0.3,type=float)
parser.add_argument("-t","--time_splits",nargs=2,help="Time of splits, in generations (default t_1 = 2e+04, t_2 = 4e+04)",default = [int(2e+04),int(4e+04)])
parser.add_argument("-bin","--bin_length",help="The length (in bases) of the how wide you want each bin. (Default is 100)",default=100,type=int)
parser.add_argument("-o_coaldir","--output_coaldir",help="Output directory of coalescent data",default=os.getcwd() + '/coal_data/',type=str)
parser.add_argument("-o_coalname","--output_coalname",help="Output name for coal_data filename",default='',type=str)
parser.add_argument("-o_mhsdir","--output_mhsdir",help="Output dir for vcf and mhs data",default=os.getcwd() + '/vcf_mhs/',type=str)
parser.add_argument("-o_mhsname","--output_mhsname",help="Output name for vcf and mhs data",default = '',type=str)
parser.add_argument("-rho","--recomb_rate",help="Rate of recombination per bp per generation",default=2e-08,type=float)
parser.add_argument("-mew","--mut_rate",help="Rate of mutation per bp per generation",default=2e-08,type=float)
parser.add_argument("--suffix_time",help="Boolean, whether to write the time in to the files (useful for multi runs of the same simulation",action="store_true")
parser.add_argument("--print",help="Print the DemographicDebugger (from msprime) and the number of segratating sites",action="store_true")
parser.add_argument("--tree",help="Draw trees and info per segment",action="store_true")
args = parser.parse_args()

model = args.model
print_ = args.print
if args.suffix_time:
    suffix = '_' + datetime.now().strftime("%H%M%S")
else:
    suffix = ''

print('Running simulation of model {}'.format(model))
command_line = 'sim = ' + model + '(' + str(print_) + ')'
exec('sim = ' + command_line)

tmrca_data = np.array(get_coal_data(sim, args)) # get coalescent data, as np array
header = "tree index, upper sequence interval,tMCRA"
np.savetxt(args.output_coaldir + args.output_coalname + model + suffix + '.txt', tmrca_data, header=header, delimiter="\t") #write to disc
print('coal_data written to {}.'.format(args.output_coaldir + args.output_coalname + suffix + '.txt'))


sim_to_mhs(sim,vcf_path = args.output_mhsdir + args.output_mhsname,mhs_path = args.output_mhsdir + args.output_mhsname, suffix = '_' + model + suffix )