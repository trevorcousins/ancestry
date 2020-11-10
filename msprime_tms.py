# msprime_TMs
# This script takes an msprime model as specified in msprime_models.py and calculates the TM for this.
# and optionally plots it.

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
parser.add_argument("-bin","--bin_length",help="The length (in bases) of the how wide you want each bin. (Default is 1000)",default=100,type=int)
parser.add_argument("-N_T","--Number_of_states",help="The number of states you want, for the time intervals (Default is 32)",default=32,type=int)
args = parser.parse_args()

model = args.model
N_T = args.Number_of_states
print('args.Number_of_states is {}'.format(args.Number_of_states))
if args.print:
    print_ = True
else: 
    print_ = False

# run simulation
print('Running simulation of model {}'.format(model))
exec('sim = ' + model + '(' + str(print_) + ')' )
print('Simulation finished')

# generate time intervals, scaling by 2 *N_0 so that they are in generations
# as per eq 13 in ref1
def scaled_time_intervals(sim, N_T = args.Number_of_states,mu = 2e-08, alpha = 0.1, Tmax = 15):
    # parameters: sim is the msprime simulation; mu is the per bp per gen mutation rate
    # ...alpha and Tmax are parameters defining the spacing of the intervals
    
    # heterozygosity is defined as the fraction of loci within an individual that are heterozygous.
    # het = sim.get_pairwise_diversity() / sim.get_sequence_length()
    het = get_het(sim)
    # N_0 (long term effective population size) is defined as heterozygosity / 4mu
    N_0 = het/(4*mu)
    # define time boundaries
    T = [0] # first lower boundary is 0 
    for i in range(0,N_T): # TODO: should T have upper interval
        T.append( alpha*math.exp( (i/N_T)*math.log(1 + Tmax/alpha) - 1))
    T_np = np.array(T) # convert to numpy
    # scale these to generations with *2*N_0, as per Schiffels' instruction
    T_scaled_np = T_np * 2 * N_0
    return T_scaled_np, T_np, N_0

# get(population) heterozygosity from simulation
def get_het(sim):
    het = sim.get_pairwise_diversity() / sim.get_sequence_length()
    return het

# round the coalescent times into their respective time intervals
def round_coal_times(coal_times,T_scaled_np):
    coal_times_intervals = [] 
    for i in range(0,len(coal_times)):
        diff = coal_times[i] - T_scaled_np
        diff_pos = diff[diff >0]
        where = np.argmin(diff_pos) # find closest interval
        if where == N_T: #TODO find a better fix for this. Have extended T to include upper interval
            where = where - 1
        coal_times_intervals.append(where) 
    return coal_times_intervals

# generate counts of transitioning from one state to another
def tm_counts(coal_times_intervals):
    tm = np.zeros(shape=(N_T,N_T))
    # tally the number of transitions from one interval to another (distributing over the rows)
    for i in range(0,len(coal_times_intervals)-1):
        col = coal_times_intervals[i+1]
        row  = coal_times_intervals[i]
        tm[row,col] += 1
    return tm

# get true coalescent times per the sequence returned in msprime; save this as a long array
def get_coal_times(simulation = sim,args = args):
    split_times = []
    coal_times = np.zeros(shape=(int(sim.sequence_length/args.bin_length))) # create empty array which will store the binned coalescent times
    ind = 0 #index for binning
    for tree in sim.trees():
        if args.tree: #draw tree with info if that is required
            print("tree {}: interval = {}".format(tree.index, tree.interval))
            print("TMRCA: ",tree.time(tree.parent(0)))
            print(tree.draw(format="unicode"))
        split_times.append(tree.time(tree.parent(0))) # append coalescent time to list
        coal_times[ind:math.floor(tree.interval[1]/args.bin_length)] = tree.time(tree.parent(0))
        ind = math.floor(tree.interval[1]/args.bin_length)
    return coal_times

def normalise(matrix,how = 'colsum'):
    options = ['rowmax', 'rowsum', 'colmax','colsum']
    if how in options:
        # actions for either row or column
        possible_actions = ['row','col']
        actions = ['[i,]','[:,i]']
        # calcs for either sum or max
        calcs = ['max','sum']
        normalised_tm = np.array(matrix)
        action = actions[possible_actions.index(how[0:3])]
        calc = calcs[calcs.index(how[3:6])]
        # check square
        if matrix.shape[0] != matrix.shape[1]:
            print('Warning! Matrix is not square. Continuing.') 
        for i in range(0,matrix.shape[1]):
            if max(matrix[:,i]) != 0.0:
                execute_me = 'normalised_tm' + action + ' = normalised_tm' + action + ' / ' + calc + '(normalised_tm' + action + ')'
                # print(execute_me)
                exec(execute_me)
        matrix = np.array(normalised_tm)
    return matrix
    
coal_times = get_coal_times(simulation = sim, args = args) # get true coalescent times
T_scaled, T, N_0 = scaled_time_intervals(sim) # get scaled time intervals
coal_times_intervals = round_coal_times(coal_times,T_scaled) # round coalescent times into their respected interval
tm = tm_counts(coal_times_intervals) # generate a matrix of counts, where each element counts the number of times state j transitions to state i for matrix [i,j]

tm_norm = normalise(tm,'colsum') # normalise this such that it represents a probability distribution

# remove diagonals. let nd stand for 'nodiagonals'
tm_nd = np.copy(tm)
tm_nd[range(0,len(tm)),range(0,len(tm))] = 0
tm_nd_norm = normalise(tm_nd,'colsum') # normalise over max of each columnn (return a probability distribution over columns)

# TODO left-most and upper-most column look strange - try and fix. Update, think this is ok
heatmaps_seq(normalise(tm_nd,'colmax'),title='Normalised true TM') # show heatmap

# TODO: check alpha and beta in abinitio, it seems matrix is coming out wrong way round
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
