# This script considers a structured population history, then creates a history of population size increases and decreases
# where the coalescent rates match. It models the scenario where we have constant population size for 0<=t<s1, then some amount of instantaneous structure
# between s1 <= t <s2, then back to regular population size with no structure at s2 <= t


# take as input migration amount and time of migration
# generate expected pop size in this interval
# generate msprime pop size changes in this
# simulate model


import os
import math
import matplotlib.pyplot as plt
import numpy as np
import pdb
import msprime
import argparse
from msprime_utils import *
from vcf_mhs import *
from datetime import datetime

parser = argparse.ArgumentParser()
parser.add_argument("-N","--N_0",help="Initial population size (Default is 1e+04",default=1e+04,type=int)
parser.add_argument("-L","--seq_length",help="Length of the sequence to be simulated (default 3e+07)",default=int(3e+07),type=int)
parser.add_argument("-mig","--migration_prop",help="Proportion of migrating population (gamma) (default 0.3)",default=0.3,type=float)
parser.add_argument("-t","--time_splits",nargs=2,help="Time of splits, in generations (default t_1 = 2e+04, t_2 = 4e+04)",default = [int(2e+04),int(4e+04)])
parser.add_argument("-time_ints","--time_intervals",help="Number of intervals to take in approximation of structured CR (default 50)",default=50,type=int)
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
parser.add_argument("--visualise",help="Show plot of ICR",action="store_true")
args = parser.parse_args()

def pdf(t,s1,gamma):
    # probability density function for lineages coalescing
    bigger_s1 = np.exp(-s1) # probability lineages didn't coalesce in 0 <= t < s1
    f = np.array(
        bigger_s1 * (
        ((1-gamma)**2)*(np.exp(-(t-s1)/(1-gamma))) + (gamma**2)*np.exp(-(t-s1)/gamma)
    )
    )
    return f
def CDF(t,s1,gamma):
# cumulative density function for probability of lineages coalescing
    F = np.array(
        1 - np.exp(-s1) + np.exp(-s1)*( ((gamma-1)**3)*np.exp(-(t-s1)/(1-gamma)) - gamma**3*np.exp(-(t-s1)/gamma) - (gamma-1)**3 + gamma**3)
    )
    return F
def ICR(t,s1,gamma):
# inverse coalescent rate 
    ICR = (1-CDF(t,s1,gamma))/(pdf(t,s1,gamma))
    return ICR

def pop_size(t,s1,s2,gamma,N_0):
    # function returns an array of population sizes

    
    # inverse coalescent rate returns Ne relative to starting size, so multiply by N_0
    pop_size = ICR(t,s1,gamma) * N_0
    return pop_size

def matching_CR(N_0,seq_length,mut_rate,recomb_rate,demographic_events_input,print_=True):
    # msprime model to implement population size changes
    # demographic_events_input is the array of population sizes that you are implementing
    # check these with the model you are matching
    N_0 = 10000
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=2, initial_size=N_0,growth_rate=0)
    ]
    migration_matrix = [[0]]
    demographic_events = demographic_events_input
    dd = msprime.DemographyDebugger(
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events)
    if print_:
        dd.print_history()
    sim = msprime.simulate(population_configurations=population_configurations,
                           migration_matrix=migration_matrix,
                           demographic_events=demographic_events, length=seq_length, recombination_rate=recomb_rate,mutation_rate=mut_rate)
    return sim

def model_check():
    # model check
    print_=True
    model = 'matching_CR'
    print('Running simulation of model {}'.format(model))
    exec('sim = ' + model + '(' + str(print_) + ')' )
    print('Simulation finished')
    return None

def pophist(N_0,s1,s2,gamma,intervals):
    # returns a string of commands which represent the population size changes (should be fed into msprime)
    # TODO: make inputs which can decide N_0, number of intervals, s1,s2
    N_0 = 10000 # TOOO write this into function, so it can be specified
    pophistlist = [] # list of strings of commands
    dem_events_list = [] # list of msprime commands 
    intervals = 50 # number of changes to make
    times_scaled = np.arange(s1,s2,(s2-s1)/intervals) # scaled times at which changes occur
    # times = [int(time) for time in times_scaled*10000] # times in generations
    times_float = np.array(times_scaled)*N_0 # times in generations TODO update to include N_0
    times = times_float.astype(int) # convert to int
    pop_sizes = pop_size(times_scaled,s1,s2,gamma,N_0) # array of population sizes
    for i in range(len(times)):
        # create a list of strings which are msprime commands
        string = 'msprime.PopulationParametersChange(time=' + str(times[i]) + ', initial_size=' + str(pop_sizes[i]) + ',growth_rate=0)'
        pophistlist.append(string)
        exec("dem_events_list.append(" + pophistlist[i] + ")")
    # return population size to N_0
    dem_events_list.append(msprime.PopulationParametersChange(time=s2*N_0, initial_size=N_0,growth_rate=0)) 
    return pophistlist, dem_events_list

def visualise(s1,s2,N_0,gamma):
    t2 = np.arange(s1,s2,(s2-s1)/1000)
    plt.plot(t2,N_0*ICR(t2,s1,gamma),'k')
    plt.show()

# default suffix. If specified, use datetime now. If not, nothing
if args.suffix_time:
    suffix = '_' + datetime.now().strftime("%H%M%S")
else:
    suffix = ''

model = 'matchingPSC'
if args.output_coalname == '': # by default, save name as Ymd + model + migprop
    coal_output = args.output_coaldir + datetime.now().strftime("%Y%m%d") + '_' + model + '_mig' + str(int(10*args.migration_prop))
else:
    coal_output = args.output_coaldir + datetime.now().strftime("%Y%m%d") + args.output_coalname 
if args.output_mhsname == '':
    mhs_output = args.output_mhsdir +  datetime.now().strftime("%Y%m%d") + '_' + model + '_mig' + str(int(10*args.migration_prop))
else: 
    mhs_output = args.output_mhsdir +  datetime.now().strftime("%Y%m%d") + args.output_mhsname


N_0 = args.N_0
s1 = int(args.time_splits[0])/N_0
s2 = int(args.time_splits[1])/N_0
gamma = args.migration_prop
intervals = args.time_intervals
# N_0, s1, s2, mig_prop (gamma)

pophistlist, dem_events_list = pophist(N_0,s1,s2,gamma,intervals)

# this works
# demographic_events_input = [msprime.PopulationParametersChange(time=131723, initial_size=314159,growth_rate=0),msprime.PopulationParametersChange(time=293137, initial_size=21718,growth_rate=0)] # this works

print('Running simulation of model...')
sim = matching_CR(N_0,args.seq_length,args.mut_rate,args.recomb_rate,dem_events_list)
print('Simulation finished.')

if args.visualise:
    visualise(s1,s2,N_0,gamma)

pdb.set_trace()
# save coalescent data to disc
#TODO fix output and then run simulations
tmrca_data = np.array(get_coal_data(sim, args)) # get coalescent data, as np array
header = "tree index, upper sequence interval,tMCRA"
np.savetxt(coal_output + suffix + '.txt', tmrca_data, header=header, delimiter="\t") #write to disc
print('coal_data written to {}.'.format(coal_output + suffix + '.txt'))

# save vcf and mhs data to disc
sim_to_mhs(sim,vcf_path = args.output_mhsdir,mhs_path = args.output_mhsdir, suffix = '_' + model + '_mig' + str(int(10*args.migration_prop)  ) + suffix )

# load this file
# test = np.loadtxt(path/to/file,comments="#")


