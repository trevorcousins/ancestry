# This script hosts different functions that are useful for msprime simulations 
import math
import numpy as np

# generate time intervals as per Wang et. al. (2020), scaling them by the N_0 (long term effective pop size)
def scaled_time_intervals(sim, N_T,mu = 2e-08, alpha = 0.1, Tmax = 15):
    # parameters: sim is the msprime simulation; mu is the per bp per gen mutation rate. N_T is number of states
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

# function to return heterozygosity of an msprime simulation
def get_het(sim):
    het = sim.get_pairwise_diversity() / sim.get_sequence_length()
    return het

# function to round a sequence of coalescent times into their closest time interval
def round_coal_times(coal_times,T_scaled_np,N_T):
    # coal_times are the (exact) coalescent times along a sequence, as given by ms prime
    # T_scaled_np are the scaled time intervals (in a np aray)
    # N_T is the number of hidden states
    coal_times_intervals = [] 
    for i in range(0,len(coal_times)):
        diff = coal_times[i] - T_scaled_np
        diff_pos = diff[diff >0]
        where = np.argmin(diff_pos) # find closest interval
        if where == N_T: #TODO find a better fix for this. Have extended T to include upper interval
            where = where - 1
        coal_times_intervals.append(where) 
    return coal_times_intervals

# function to count the number of transitions between different states in a sequences of tMRCAs 
def tm_counts(coal_times_intervals,N_T):
    tm = np.zeros(shape=(N_T,N_T))
    # tally the number of transitions from one interval to another (distributing over the rows)
    for i in range(0,len(coal_times_intervals)-1):
        col = coal_times_intervals[i+1]
        row  = coal_times_intervals[i]
        tm[row,col] += 1
    return tm

# function to return an array of coalescent data (tree index, intervals separated by recombination, tMRCA time) from an msprime simulation.
def get_coal_data(sim,args):
    tmrca_data = [] # first tuple is index, second tuple is sequence interval, third tuple is tMRCA
    ind = 0 #index for binning
    for tree in sim.trees():
        if args.tree: #draw tree with info if that is required
            print("tree {}: interval = {}".format(tree.index, tree.interval))
            print("TMRCA: ",tree.time(tree.parent(0)))
            print(tree.draw(format="unicode"))
        tmrca_data.append((tree.index,tree.interval[1],tree.time(tree.parent(0))))
        ind = math.floor(tree.interval[1]/args.bin_length)
    return tmrca_data

# function to divide the true coalescent data into bins of defined size
def round_bin_coal_data(sim,tmrca_data,args):
    coal_times = np.zeros(shape=(int(sim.sequence_length/args.bin_length))) # create empty array which will store the binned coalescent times
    ind = 0 #index for binning
    for i in range(0,len(tmrca_data)):
        coal_times[ind:math.floor(tmrca_data[i][1]/args.bin_length)] = tmrca_data[i][2]
        ind = math.floor(tmrca_data[i][1]/args.bin_length)
    return coal_times

# function to normalise a matrix. Can be normalised over either rows or columns, using either the sum or max
# use sum for a probability distribution; use max for a good visual interpretation
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

def num_occcurences(array,time_intervals):
    counts = np.zeros(shape=len(time_intervals))
    for i in range(len(array)):
        counts[array[i]] += 1
    return counts

