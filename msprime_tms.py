# msprime_TMs
# This script takes an msprime model as specified in msprime_models.py and calculates the TM for this.
# and optionally plots it.
# example usage:
# $ python msprime_tms.py m0001 --print --tree -bin 1000 -N_T 32

# 1. simulate history and save data 
# 2. read coalescent times into an array
# 2.5 partition sequence into bins of defined length
# 3. create time intervals, then round this array of coalescent times into these intervals
# 4. iterate along these rounded time intervals and calculate the probability of transitioning from one state (time interval) to another, recording info in a matrix
# 5. (optional) print this matrix

# ref 1 is https://journals.plos.org/plosgenetics/article/file?id=10.1371/journal.pgen.1008552&type=printable

import argparse
from msprime_models import * 
# from abinitio_tm import scaled_time_intervals
from heatmaps_generate import heatmaps_seq, heatmaps_div
from abinitio_tm import abinitio
# from msprime_utils import scaled_time_intervals, get_het, round_coal_times, tm_counts, get_coal_data, round_bin_coal_data, normalise
from msprime_utils import *
import numpy as np
import pdb
import math
import pandas as pd
from scipy.stats import entropy

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("model",help="Specify the model to be used, as defined in msprime_models.")
parser.add_argument("--print",help="Print the DemographicDebugger (from msprime) and the number of segratating sites",action="store_true")
parser.add_argument("--tree",help="Draw trees and info per segment",action="store_true")
parser.add_argument("-bin","--bin_length",help="The length (in bases) of the how wide you want each bin. (Default is 100)",default=100,type=int)
parser.add_argument("-N_T","--Number_of_states",help="The number of states you want, for the time intervals (Default is 32)",default=32,type=int)
args = parser.parse_args()

# save args
model = args.model
N_T = args.Number_of_states
print('args.Number_of_states is {}'.format(args.Number_of_states))
print('args.bin_length is {}'.format(args.bin_length))
if args.print:
    print_ = True
else: 
    print_ = False

# run simulation
print('Running simulation of model {}'.format(model))
exec('sim = ' + model + '(' + str(print_) + ')' )
print('Simulation finished')
  
tmrca_data = get_coal_data(sim, args) # get true coalescent 
coal_times = round_bin_coal_data(sim.sequence_length,tmrca_data, args) # partition coalescent data into bins  
T_scaled, T, N_0 = scaled_time_intervals(sim,N_T = args.Number_of_states) # get scaled time intervals
coal_times_intervals = round_coal_times(coal_times,T_scaled,N_T=args.Number_of_states) # round binned coalescent times into their respected interval
tm = tm_counts(coal_times_intervals,N_T=args.Number_of_states) # generate a matrix of counts, where each element counts the number of times state j transitions to state i for matrix [i,j]

tm_norm = normalise(tm,'colsum') # normalise this such that it represents a probability distribution

# remove diagonals. let nd stand for 'nodiagonals'
tm_nd = np.copy(tm)
tm_nd[range(0,len(tm)),range(0,len(tm))] = 0
tm_nd_norm = normalise(tm_nd,'colsum') # normalise over max of each columnn (return a probability distribution over columns)

# TODO left-most and upper-most column look strange - try and fix. Update, think this is ok
heatmaps_seq(normalise(tm_nd,'colmax'),title='Normalised true TM') # show heatmap

pop_size = [10000 for i in range(0,50)] # set array of population sizes (for the lambda array)
ab_tm_nd, ab_tm = abinitio(pop_size,T,N_0,N_T) # generate theoretical probabilities of this matrix 
heatmaps_seq(normalise(ab_tm_nd,'colmax'),title='abintio TM') # normalise this matrix and then plot heatmap
ab_tm_norm = normalise(ab_tm,'colsum') # normalise probabilities
ab_tm_nd_norm = normalise(ab_tm_nd, 'colsum')

# compare probability distributions (With diagonals)
diff_all = tm_norm - ab_tm_norm
# compare without diagonals
diff_nd = tm_nd_norm - ab_tm_nd_norm

heatmaps_div(diff_nd,'Difference between true & abinitio TM')

# calculate KL divergence for each column
KL_divergence = [entropy(tm_nd_norm[:,i],ab_tm_nd_norm[:,i]) for i in range(0,N_T)]

pdb.set_trace()
sys.exit()

