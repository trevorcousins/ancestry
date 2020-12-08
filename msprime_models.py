# This script details different msprime models, which are saved as functions. 
# These can then be executed in other files with greater ease.
# There should be a short description above each history.
# Times should be given in generations

import math
import msprime
import sys
import pdb
import os
import numpy as np


# simple population history, constant population size with no growth or migration
# one diploid individual
def m0001(print_):
    N_A = 10000
    T_A1 = 20000 
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=2, initial_size=N_A,growth_rate=0)
    ]
    migration_matrix = [
        [0],
        ]
    demographic_events = [
    ]
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    dd = msprime.DemographyDebugger(
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events)
    if print_:
        dd.print_history()
    sim = msprime.simulate(population_configurations=population_configurations,
                           migration_matrix=migration_matrix,
                           demographic_events=demographic_events, length=150e+6, recombination_rate=2e-8,mutation_rate=2e-8)
    return sim

# population that double in size at some time
def m0002(print_):
    N_A = 10000
    T_A1 = 20000 
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=2, initial_size=N_A,growth_rate=0)
    ]
    migration_matrix = [
        [0],
        ]
    demographic_events = [
        msprime.PopulationParametersChange(time=T_A1, initial_size=N_A*2)
    ]
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    dd = msprime.DemographyDebugger(
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events)
    if print_:
        dd.print_history()
    sim = msprime.simulate(population_configurations=population_configurations,
                           migration_matrix=migration_matrix,
                           demographic_events=demographic_events, length=150e+6, recombination_rate=2e-8,mutation_rate=2e-8)
    return sim

# model with instantaneous migration 
# TODO set so it matches instant_structure0001, where N_B is dependent on migration rate
# TODO set so it migration rate can be specifed through command line (like instant_structure0001)
def struct0001(print_):
    N_A0 = 1e+04
    N_B0 = 1e+04
    T_A = 2e+04
    T_B = 3e+04
    # initially.
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=2, initial_size=N_A0, growth_rate=0),
        msprime.PopulationConfiguration(
            sample_size=0, initial_size=N_B0, growth_rate=0),
    ]
    migration_matrix = [[0,0],[0,0]]
    demographic_events = [
        msprime.MassMigration(
            time=T_A, source = 0, dest=1,proportion=1),
        msprime.MassMigration(
            time=T_B, source=1, dest=0, proportion=1)
    ]
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    dd = msprime.DemographyDebugger(
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events)
    if print_:
        dd.print_history()
    sim = msprime.simulate(population_configurations=population_configurations,
                           migration_matrix=migration_matrix,
                           demographic_events=demographic_events, length=3e+6, recombination_rate=2e-8,
                           mutation_rate=2e-8,random_seed=50)
    return sim


# model with instant migration which can be called through command line. See msprime_generate.py 
# set mig_prop to your desired proportion of instantaneous migration
def instant_struct0001(N_0,mig_prop,t_1,t_2,seq_length,recomb_rate,mut_rate,print_):
    
    if not 0 < mig_prop < 1:
        print('Error. mig_prop is {} which is not between 0 and 1. Aborting\n'.format(mig_prop))
        sys.exit()
    print('N_0 is {}'.format(N_0))
    print('mig_prop is {}'.format(mig_prop))
    print('t_1 is {} and t_2 is {}'.format(t_1,t_2))
    print('seq_length is {}'.format(seq_length))
    print('recomb_rate is {}'.format(recomb_rate))
    print('mut_rate is {}'.format(mut_rate))

    # proportion of instantaneous migration
    print('migration proportion is {}'.format(mig_prop))
    N_B0 = mig_prop * 1e+04
    # initially.
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=2, initial_size=N_0, growth_rate=0),
        msprime.PopulationConfiguration(
            sample_size=0, initial_size=N_B0, growth_rate=0),
    ]
    demographic_events = [
        msprime.MassMigration(time=t_1, source = 0, dest=1,proportion=mig_prop),
        msprime.MassMigration(time=t_2, source=1, dest=0, proportion=1)
    ]
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    dd = msprime.DemographyDebugger(
        population_configurations=population_configurations,
        demographic_events=demographic_events)
    if print_:
        print('Demographic history:\n')
        dd.print_history()
    sim = msprime.simulate(population_configurations=population_configurations,
                           demographic_events=demographic_events, length=seq_length, recombination_rate=recomb_rate,
                           mutation_rate=mut_rate)
    return sim

# const population which can be called through command line. See msprime_generate.py 
def const0001(N_0,mig_prop,t_1,t_2,seq_length,recomb_rate,mut_rate,print_):
# def instant_struct0001(N_0,mig_prop,t_1,t_2,seq_length,recomb_rate,mut_rate,print_):
    print('N_0 is {}'.format(N_0))
    print('mig_prop is {}, but is irrelevant here'.format(mig_prop))
    print('t_1 is {} and t_2 is {}, but these are irrelevant here'.format(t_1,t_2))
    print('seq_length is {}'.format(seq_length))
    print('recomb_rate is {}'.format(recomb_rate))
    print('mut_rate is {}'.format(mut_rate))
    N_B0 = mig_prop * 1e+04
    # initially.
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=2, initial_size=N_0, growth_rate=0),
    ]
    demographic_events = []
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    dd = msprime.DemographyDebugger(
        population_configurations=population_configurations,
        demographic_events=demographic_events)
    if print_:
        print('Demographic history:\n')
        dd.print_history()
    sim = msprime.simulate(population_configurations=population_configurations,
                           demographic_events=demographic_events, length=seq_length, recombination_rate=recomb_rate,
                           mutation_rate=mut_rate)
    return sim

def mazet2016_3b(print_):
    N_0 = 10000
    # initially.
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=2, initial_size=N_0, growth_rate=0),
    ]
    demographic_events = [
        msprime.PopulationParametersChange(time=10000, initial_size=N_0/10,growth_rate=0),
        msprime.PopulationParametersChange(time=60000, initial_size=N_0,growth_rate=0),
        msprime.PopulationParametersChange(time=200000, initial_size=N_0/2,growth_rate=0)
    ]
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    dd = msprime.DemographyDebugger(
        population_configurations=population_configurations,
        demographic_events=demographic_events)
    if print_:
        print('Demographic history:\n')
        dd.print_history()
    sim = msprime.simulate(population_configurations=population_configurations,
                           demographic_events=demographic_events, length=150e+6, recombination_rate=2e-08,
                           mutation_rate=2e-08)
    return sim

def matching_mig01(print_):
    # this history of population size change matches the coalescent rate as detected by PSMC
    # for a structured history, with mig rate =01 (see my_log.txt)
    pwd = os.getcwd()
    if 'trevor' in pwd:
        hpc = False
    else: 
        hpc = True
    
    # load file, different depending on hpc or Surface
    if not hpc:
        path = '/home/trevor/ancestry/msmc_out/mig01/mig01_all.final.txt'
    else: 
        path = '/home/tc557/ancestry/msmc_out/mig01/mig01_all.final.txt'
     # load data

    mu = 2e-08
    data = np.loadtxt(path,skiprows=1)
    # load coalescent rate column
    lambda_ = data[:,3]
    # convert this into pop size estimates and get time intervals
    # check mu matches simulation that was used to create!
    popsize_estimate = (1/data[:,3])/(2*mu)
    times = data[:,1]/mu
    

    N_0 = 10000
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=2, initial_size=N_0, growth_rate=0),
    ]
    demographic_events = [
        msprime.PopulationParametersChange(time=times[0], initial_size=popsize_estimate[0]),
        msprime.PopulationParametersChange(time=times[1], initial_size=popsize_estimate[1]),
        msprime.PopulationParametersChange(time=times[2], initial_size=popsize_estimate[2]),
        msprime.PopulationParametersChange(time=times[3], initial_size=popsize_estimate[3]),
        msprime.PopulationParametersChange(time=times[4], initial_size=popsize_estimate[4]),
        msprime.PopulationParametersChange(time=times[5], initial_size=popsize_estimate[5]),
        msprime.PopulationParametersChange(time=times[6], initial_size=popsize_estimate[6]),
        msprime.PopulationParametersChange(time=times[7], initial_size=popsize_estimate[7]),
        msprime.PopulationParametersChange(time=times[8], initial_size=popsize_estimate[8]),
        msprime.PopulationParametersChange(time=times[9], initial_size=popsize_estimate[9]),
        msprime.PopulationParametersChange(time=times[10], initial_size=popsize_estimate[10]),
        msprime.PopulationParametersChange(time=times[11], initial_size=popsize_estimate[11]),
        msprime.PopulationParametersChange(time=times[12], initial_size=popsize_estimate[12]),
        msprime.PopulationParametersChange(time=times[13], initial_size=popsize_estimate[13]),
        msprime.PopulationParametersChange(time=times[14], initial_size=popsize_estimate[14]),
        msprime.PopulationParametersChange(time=times[15], initial_size=popsize_estimate[15]),
        msprime.PopulationParametersChange(time=times[16], initial_size=popsize_estimate[16]),
        msprime.PopulationParametersChange(time=times[17], initial_size=popsize_estimate[17]),
        msprime.PopulationParametersChange(time=times[18], initial_size=popsize_estimate[18]),
        msprime.PopulationParametersChange(time=times[19], initial_size=popsize_estimate[19]),
        msprime.PopulationParametersChange(time=times[20], initial_size=popsize_estimate[20]),
        msprime.PopulationParametersChange(time=times[21], initial_size=popsize_estimate[21]),
        msprime.PopulationParametersChange(time=times[22], initial_size=popsize_estimate[22]),
        msprime.PopulationParametersChange(time=times[23], initial_size=popsize_estimate[23]),
        msprime.PopulationParametersChange(time=times[24], initial_size=popsize_estimate[24]),
        msprime.PopulationParametersChange(time=times[25], initial_size=popsize_estimate[25]),
        msprime.PopulationParametersChange(time=times[26], initial_size=popsize_estimate[26]),
        msprime.PopulationParametersChange(time=times[27], initial_size=popsize_estimate[27]),
        msprime.PopulationParametersChange(time=times[28], initial_size=popsize_estimate[28]),
        msprime.PopulationParametersChange(time=times[29], initial_size=popsize_estimate[29]),
        msprime.PopulationParametersChange(time=times[30], initial_size=popsize_estimate[30]),
        msprime.PopulationParametersChange(time=times[31], initial_size=popsize_estimate[31])
    ]
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    dd = msprime.DemographyDebugger(
        population_configurations=population_configurations,
        demographic_events=demographic_events)
    if print_:
        print('Demographic history:\n')
        dd.print_history()
    sim = msprime.simulate(population_configurations=population_configurations,
                           demographic_events=demographic_events, length=150e+6, recombination_rate=2e-08,
                           mutation_rate=mu)
    return sim

def matching_mig02(print_):
    # this history of population size change matches the coalescent rate as detected by PSMC
    # for a structured history, with mig rate =01 (see my_log.txt)
    pwd = os.getcwd()
    if 'trevor' in pwd:
        hpc = False
    else: 
        hpc = True
    
    # load file, different depending on hpc or Surface
    if not hpc:
        path = '/home/trevor/ancestry/msmc_out/mig02/mig02_all.final.txt'
    else: 
        path = '/home/tc557/ancestry/msmc_out/mig02/mig02_all.final.txt'
    # load data

    mu = 2e-08
    data = np.loadtxt(path,skiprows=1)
    # load coalescent rate column
    lambda_ = data[:,3]
    # convert this into pop size estimates and get time intervals
    # check mu matches simulation that was used to create!
    popsize_estimate = (1/data[:,3])/(2*mu)
    times = data[:,1]/mu
    

    N_0 = 10000
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=2, initial_size=N_0, growth_rate=0),
    ]
    demographic_events = [
        msprime.PopulationParametersChange(time=times[0], initial_size=popsize_estimate[0]),
        msprime.PopulationParametersChange(time=times[1], initial_size=popsize_estimate[1]),
        msprime.PopulationParametersChange(time=times[2], initial_size=popsize_estimate[2]),
        msprime.PopulationParametersChange(time=times[3], initial_size=popsize_estimate[3]),
        msprime.PopulationParametersChange(time=times[4], initial_size=popsize_estimate[4]),
        msprime.PopulationParametersChange(time=times[5], initial_size=popsize_estimate[5]),
        msprime.PopulationParametersChange(time=times[6], initial_size=popsize_estimate[6]),
        msprime.PopulationParametersChange(time=times[7], initial_size=popsize_estimate[7]),
        msprime.PopulationParametersChange(time=times[8], initial_size=popsize_estimate[8]),
        msprime.PopulationParametersChange(time=times[9], initial_size=popsize_estimate[9]),
        msprime.PopulationParametersChange(time=times[10], initial_size=popsize_estimate[10]),
        msprime.PopulationParametersChange(time=times[11], initial_size=popsize_estimate[11]),
        msprime.PopulationParametersChange(time=times[12], initial_size=popsize_estimate[12]),
        msprime.PopulationParametersChange(time=times[13], initial_size=popsize_estimate[13]),
        msprime.PopulationParametersChange(time=times[14], initial_size=popsize_estimate[14]),
        msprime.PopulationParametersChange(time=times[15], initial_size=popsize_estimate[15]),
        msprime.PopulationParametersChange(time=times[16], initial_size=popsize_estimate[16]),
        msprime.PopulationParametersChange(time=times[17], initial_size=popsize_estimate[17]),
        msprime.PopulationParametersChange(time=times[18], initial_size=popsize_estimate[18]),
        msprime.PopulationParametersChange(time=times[19], initial_size=popsize_estimate[19]),
        msprime.PopulationParametersChange(time=times[20], initial_size=popsize_estimate[20]),
        msprime.PopulationParametersChange(time=times[21], initial_size=popsize_estimate[21]),
        msprime.PopulationParametersChange(time=times[22], initial_size=popsize_estimate[22]),
        msprime.PopulationParametersChange(time=times[23], initial_size=popsize_estimate[23]),
        msprime.PopulationParametersChange(time=times[24], initial_size=popsize_estimate[24]),
        msprime.PopulationParametersChange(time=times[25], initial_size=popsize_estimate[25]),
        msprime.PopulationParametersChange(time=times[26], initial_size=popsize_estimate[26]),
        msprime.PopulationParametersChange(time=times[27], initial_size=popsize_estimate[27]),
        msprime.PopulationParametersChange(time=times[28], initial_size=popsize_estimate[28]),
        msprime.PopulationParametersChange(time=times[29], initial_size=popsize_estimate[29]),
        msprime.PopulationParametersChange(time=times[30], initial_size=popsize_estimate[30]),
        msprime.PopulationParametersChange(time=times[31], initial_size=popsize_estimate[31])
    ]
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    dd = msprime.DemographyDebugger(
        population_configurations=population_configurations,
        demographic_events=demographic_events)
    if print_:
        print('Demographic history:\n')
        dd.print_history()
    sim = msprime.simulate(population_configurations=population_configurations,
                           demographic_events=demographic_events, length=150e+6, recombination_rate=2e-08,
                           mutation_rate=mu)
    return sim

def matching_mig03(print_):
    # this history of population size change matches the coalescent rate as detected by PSMC
    # for a structured history, with mig rate =01 (see my_log.txt)
    pwd = os.getcwd()
    if 'trevor' in pwd:
        hpc = False
    else: 
        hpc = True
    
    # load file, different depending on hpc or Surface
    if not hpc:
        path = '/home/trevor/ancestry/msmc_out/mig03/mig03_all.final.txt'
    else: 
        path = '/home/tc557/ancestry/msmc_out/mig03/mig03_all.final.txt'
    # load data

    mu = 2e-08
    data = np.loadtxt(path,skiprows=1)
    # load coalescent rate column
    lambda_ = data[:,3]
    # convert this into pop size estimates and get time intervals
    # check mu matches simulation that was used to create!
    popsize_estimate = (1/data[:,3])/(2*mu)
    times = data[:,1]/mu
    

    N_0 = 10000
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=2, initial_size=N_0, growth_rate=0),
    ]
    demographic_events = [
        msprime.PopulationParametersChange(time=times[0], initial_size=popsize_estimate[0]),
        msprime.PopulationParametersChange(time=times[1], initial_size=popsize_estimate[1]),
        msprime.PopulationParametersChange(time=times[2], initial_size=popsize_estimate[2]),
        msprime.PopulationParametersChange(time=times[3], initial_size=popsize_estimate[3]),
        msprime.PopulationParametersChange(time=times[4], initial_size=popsize_estimate[4]),
        msprime.PopulationParametersChange(time=times[5], initial_size=popsize_estimate[5]),
        msprime.PopulationParametersChange(time=times[6], initial_size=popsize_estimate[6]),
        msprime.PopulationParametersChange(time=times[7], initial_size=popsize_estimate[7]),
        msprime.PopulationParametersChange(time=times[8], initial_size=popsize_estimate[8]),
        msprime.PopulationParametersChange(time=times[9], initial_size=popsize_estimate[9]),
        msprime.PopulationParametersChange(time=times[10], initial_size=popsize_estimate[10]),
        msprime.PopulationParametersChange(time=times[11], initial_size=popsize_estimate[11]),
        msprime.PopulationParametersChange(time=times[12], initial_size=popsize_estimate[12]),
        msprime.PopulationParametersChange(time=times[13], initial_size=popsize_estimate[13]),
        msprime.PopulationParametersChange(time=times[14], initial_size=popsize_estimate[14]),
        msprime.PopulationParametersChange(time=times[15], initial_size=popsize_estimate[15]),
        msprime.PopulationParametersChange(time=times[16], initial_size=popsize_estimate[16]),
        msprime.PopulationParametersChange(time=times[17], initial_size=popsize_estimate[17]),
        msprime.PopulationParametersChange(time=times[18], initial_size=popsize_estimate[18]),
        msprime.PopulationParametersChange(time=times[19], initial_size=popsize_estimate[19]),
        msprime.PopulationParametersChange(time=times[20], initial_size=popsize_estimate[20]),
        msprime.PopulationParametersChange(time=times[21], initial_size=popsize_estimate[21]),
        msprime.PopulationParametersChange(time=times[22], initial_size=popsize_estimate[22]),
        msprime.PopulationParametersChange(time=times[23], initial_size=popsize_estimate[23]),
        msprime.PopulationParametersChange(time=times[24], initial_size=popsize_estimate[24]),
        msprime.PopulationParametersChange(time=times[25], initial_size=popsize_estimate[25]),
        msprime.PopulationParametersChange(time=times[26], initial_size=popsize_estimate[26]),
        msprime.PopulationParametersChange(time=times[27], initial_size=popsize_estimate[27]),
        msprime.PopulationParametersChange(time=times[28], initial_size=popsize_estimate[28]),
        msprime.PopulationParametersChange(time=times[29], initial_size=popsize_estimate[29]),
        msprime.PopulationParametersChange(time=times[30], initial_size=popsize_estimate[30]),
        msprime.PopulationParametersChange(time=times[31], initial_size=popsize_estimate[31])
    ]
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    dd = msprime.DemographyDebugger(
        population_configurations=population_configurations,
        demographic_events=demographic_events)
    if print_:
        print('Demographic history:\n')
        dd.print_history()
    sim = msprime.simulate(population_configurations=population_configurations,
                           demographic_events=demographic_events, length=150e+6, recombination_rate=2e-08,
                           mutation_rate=mu)
    return sim

def matching_mig04(print_):
    # this history of population size change matches the coalescent rate as detected by PSMC
    # for a structured history, with mig rate =01 (see my_log.txt)
    pwd = os.getcwd()
    if 'trevor' in pwd:
        hpc = False
    else: 
        hpc = True
    
    # load file, different depending on hpc or Surface
    if not hpc:
        path = '/home/trevor/ancestry/msmc_out/mig04/mig04_all.final.txt'
    else: 
        path = '/home/tc557/ancestry/msmc_out/mig04/mig04_all.final.txt'
    # load data

    mu = 2e-08
    data = np.loadtxt(path,skiprows=1)
    # load coalescent rate column
    lambda_ = data[:,3]
    # convert this into pop size estimates and get time intervals
    # check mu matches simulation that was used to create!
    popsize_estimate = (1/data[:,3])/(2*mu)
    times = data[:,1]/mu
    

    N_0 = 10000
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=2, initial_size=N_0, growth_rate=0),
    ]
    demographic_events = [
        msprime.PopulationParametersChange(time=times[0], initial_size=popsize_estimate[0]),
        msprime.PopulationParametersChange(time=times[1], initial_size=popsize_estimate[1]),
        msprime.PopulationParametersChange(time=times[2], initial_size=popsize_estimate[2]),
        msprime.PopulationParametersChange(time=times[3], initial_size=popsize_estimate[3]),
        msprime.PopulationParametersChange(time=times[4], initial_size=popsize_estimate[4]),
        msprime.PopulationParametersChange(time=times[5], initial_size=popsize_estimate[5]),
        msprime.PopulationParametersChange(time=times[6], initial_size=popsize_estimate[6]),
        msprime.PopulationParametersChange(time=times[7], initial_size=popsize_estimate[7]),
        msprime.PopulationParametersChange(time=times[8], initial_size=popsize_estimate[8]),
        msprime.PopulationParametersChange(time=times[9], initial_size=popsize_estimate[9]),
        msprime.PopulationParametersChange(time=times[10], initial_size=popsize_estimate[10]),
        msprime.PopulationParametersChange(time=times[11], initial_size=popsize_estimate[11]),
        msprime.PopulationParametersChange(time=times[12], initial_size=popsize_estimate[12]),
        msprime.PopulationParametersChange(time=times[13], initial_size=popsize_estimate[13]),
        msprime.PopulationParametersChange(time=times[14], initial_size=popsize_estimate[14]),
        msprime.PopulationParametersChange(time=times[15], initial_size=popsize_estimate[15]),
        msprime.PopulationParametersChange(time=times[16], initial_size=popsize_estimate[16]),
        msprime.PopulationParametersChange(time=times[17], initial_size=popsize_estimate[17]),
        msprime.PopulationParametersChange(time=times[18], initial_size=popsize_estimate[18]),
        msprime.PopulationParametersChange(time=times[19], initial_size=popsize_estimate[19]),
        msprime.PopulationParametersChange(time=times[20], initial_size=popsize_estimate[20]),
        msprime.PopulationParametersChange(time=times[21], initial_size=popsize_estimate[21]),
        msprime.PopulationParametersChange(time=times[22], initial_size=popsize_estimate[22]),
        msprime.PopulationParametersChange(time=times[23], initial_size=popsize_estimate[23]),
        msprime.PopulationParametersChange(time=times[24], initial_size=popsize_estimate[24]),
        msprime.PopulationParametersChange(time=times[25], initial_size=popsize_estimate[25]),
        msprime.PopulationParametersChange(time=times[26], initial_size=popsize_estimate[26]),
        msprime.PopulationParametersChange(time=times[27], initial_size=popsize_estimate[27]),
        msprime.PopulationParametersChange(time=times[28], initial_size=popsize_estimate[28]),
        msprime.PopulationParametersChange(time=times[29], initial_size=popsize_estimate[29]),
        msprime.PopulationParametersChange(time=times[30], initial_size=popsize_estimate[30]),
        msprime.PopulationParametersChange(time=times[31], initial_size=popsize_estimate[31])
    ]
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    dd = msprime.DemographyDebugger(
        population_configurations=population_configurations,
        demographic_events=demographic_events)
    if print_:
        print('Demographic history:\n')
        dd.print_history()
    sim = msprime.simulate(population_configurations=population_configurations,
                           demographic_events=demographic_events, length=150e+6, recombination_rate=2e-08,
                           mutation_rate=mu)
    return sim

def matching_mig05(print_):
    # this history of population size change matches the coalescent rate as detected by PSMC
    # for a structured history, with mig rate =01 (see my_log.txt)
    pwd = os.getcwd()
    if 'trevor' in pwd:
        hpc = False
    else: 
        hpc = True
    
    # load file, different depending on hpc or Surface
    if not hpc:
        path = '/home/trevor/ancestry/msmc_out/mig05/mig05_all.final.txt'
    else: 
        path = '/home/tc557/ancestry/msmc_out/mig05/mig05_all.final.txt'
    # load data

    mu = 2e-08
    data = np.loadtxt(path,skiprows=1)
    # load coalescent rate column
    lambda_ = data[:,3]
    # convert this into pop size estimates and get time intervals
    # check mu matches simulation that was used to create!
    popsize_estimate = (1/data[:,3])/(2*mu)
    times = data[:,1]/mu
    

    N_0 = 10000
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=2, initial_size=N_0, growth_rate=0),
    ]
    demographic_events = [
        msprime.PopulationParametersChange(time=times[0], initial_size=popsize_estimate[0]),
        msprime.PopulationParametersChange(time=times[1], initial_size=popsize_estimate[1]),
        msprime.PopulationParametersChange(time=times[2], initial_size=popsize_estimate[2]),
        msprime.PopulationParametersChange(time=times[3], initial_size=popsize_estimate[3]),
        msprime.PopulationParametersChange(time=times[4], initial_size=popsize_estimate[4]),
        msprime.PopulationParametersChange(time=times[5], initial_size=popsize_estimate[5]),
        msprime.PopulationParametersChange(time=times[6], initial_size=popsize_estimate[6]),
        msprime.PopulationParametersChange(time=times[7], initial_size=popsize_estimate[7]),
        msprime.PopulationParametersChange(time=times[8], initial_size=popsize_estimate[8]),
        msprime.PopulationParametersChange(time=times[9], initial_size=popsize_estimate[9]),
        msprime.PopulationParametersChange(time=times[10], initial_size=popsize_estimate[10]),
        msprime.PopulationParametersChange(time=times[11], initial_size=popsize_estimate[11]),
        msprime.PopulationParametersChange(time=times[12], initial_size=popsize_estimate[12]),
        msprime.PopulationParametersChange(time=times[13], initial_size=popsize_estimate[13]),
        msprime.PopulationParametersChange(time=times[14], initial_size=popsize_estimate[14]),
        msprime.PopulationParametersChange(time=times[15], initial_size=popsize_estimate[15]),
        msprime.PopulationParametersChange(time=times[16], initial_size=popsize_estimate[16]),
        msprime.PopulationParametersChange(time=times[17], initial_size=popsize_estimate[17]),
        msprime.PopulationParametersChange(time=times[18], initial_size=popsize_estimate[18]),
        msprime.PopulationParametersChange(time=times[19], initial_size=popsize_estimate[19]),
        msprime.PopulationParametersChange(time=times[20], initial_size=popsize_estimate[20]),
        msprime.PopulationParametersChange(time=times[21], initial_size=popsize_estimate[21]),
        msprime.PopulationParametersChange(time=times[22], initial_size=popsize_estimate[22]),
        msprime.PopulationParametersChange(time=times[23], initial_size=popsize_estimate[23]),
        msprime.PopulationParametersChange(time=times[24], initial_size=popsize_estimate[24]),
        msprime.PopulationParametersChange(time=times[25], initial_size=popsize_estimate[25]),
        msprime.PopulationParametersChange(time=times[26], initial_size=popsize_estimate[26]),
        msprime.PopulationParametersChange(time=times[27], initial_size=popsize_estimate[27]),
        msprime.PopulationParametersChange(time=times[28], initial_size=popsize_estimate[28]),
        msprime.PopulationParametersChange(time=times[29], initial_size=popsize_estimate[29]),
        msprime.PopulationParametersChange(time=times[30], initial_size=popsize_estimate[30]),
        msprime.PopulationParametersChange(time=times[31], initial_size=popsize_estimate[31])
    ]
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    dd = msprime.DemographyDebugger(
        population_configurations=population_configurations,
        demographic_events=demographic_events)
    if print_:
        print('Demographic history:\n')
        dd.print_history()
    sim = msprime.simulate(population_configurations=population_configurations,
                           demographic_events=demographic_events, length=150e+6, recombination_rate=2e-08,
                           mutation_rate=mu)
    return sim

def const_mig0001(print_):
    N_A0 = 1e+04
    N_B0 =  1e+04
    T_1 = 2e+04
    T_2 = 2.5e+04
    m = 4e-03
    
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=2, initial_size=N_A0, growth_rate=0),
        msprime.PopulationConfiguration(
            sample_size=0, initial_size=N_B0, growth_rate=0)
    ]
    migration_matrix = [[0,0],[0,0]]
    demographic_events = [
        msprime.MigrationRateChange(time = T_1,rate = m, matrix_index=(0,1)),
        msprime.MigrationRateChange(time = T_1,rate = m, matrix_index=(1,0)),
        msprime.MigrationRateChange(time = T_2,rate = 0, matrix_index=(0,1)),
        msprime.MigrationRateChange(time = T_2,rate = 0, matrix_index=(1,0)),
        msprime.MassMigration(time=T_2, source =1, destination = 0, proportion = 1)
    ]
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    dd = msprime.DemographyDebugger(
        population_configurations=population_configurations,migration_matrix=migration_matrix,
        demographic_events=demographic_events)
    if print_:
        print('Demographic history:\n')
        dd.print_history()
    sim = msprime.simulate(population_configurations=population_configurations,
                           demographic_events=demographic_events, length=150e+06, recombination_rate=2e-08,
                           mutation_rate=2e-08)
    return sim

def const_mig0002(print_):
    N_A0 = 1e+04
    N_B0 =  1e+04
    T_1 = 2e+04
    T_2 = 2.5e+04
    m = 4e-03
    
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=2, initial_size=N_A0, growth_rate=0),
        msprime.PopulationConfiguration(
            sample_size=0, initial_size=N_B0, growth_rate=0)
    ]
    migration_matrix = [[0,0],[0,0]]
    demographic_events = [
        msprime.MigrationRateChange(time = T_1,rate = m, matrix_index=(0,1)),
        msprime.MigrationRateChange(time = T_1,rate = m, matrix_index=(1,0)),
        msprime.MigrationRateChange(time = T_2,rate = 0, matrix_index=(0,1)),
        msprime.MigrationRateChange(time = T_2,rate = 0, matrix_index=(1,0)),
        msprime.MassMigration(time=T_2, source =0, destination =1, proportion = 1)
    ]
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    dd = msprime.DemographyDebugger(
        population_configurations=population_configurations,migration_matrix=migration_matrix,
        demographic_events=demographic_events)
    if print_:
        print('Demographic history:\n')
        dd.print_history()
    sim = msprime.simulate(population_configurations=population_configurations,
                           demographic_events=demographic_events, length=150e+06, recombination_rate=2e-08,
                           mutation_rate=2e-08)
    return sim


def const_mig0003(print_):
    N_A0 = 1e+04
    N_B0 =  1e+04
    T_1 = 2e+04
    T_2 = 4e+04
    m = 4e-03
    
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=2, initial_size=N_A0, growth_rate=0),
        msprime.PopulationConfiguration(
            sample_size=0, initial_size=N_B0, growth_rate=0)
    ]
    migration_matrix = [[0,0],[0,0]]
    demographic_events = [
        msprime.MigrationRateChange(time = T_1,rate = m, matrix_index=(0,1)),
        msprime.MigrationRateChange(time = T_1,rate = m, matrix_index=(1,0)),
        msprime.MigrationRateChange(time = T_2,rate = 0, matrix_index=(0,1)),
        msprime.MigrationRateChange(time = T_2,rate = 0, matrix_index=(1,0)),
        msprime.MassMigration(time=T_2, source =0, destination =1, proportion = 1)
    ]
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    dd = msprime.DemographyDebugger(
        population_configurations=population_configurations,migration_matrix=migration_matrix,
        demographic_events=demographic_events)
    if print_:
        print('Demographic history:\n')
        dd.print_history()
    sim = msprime.simulate(population_configurations=population_configurations,
                           demographic_events=demographic_events, length=150e+06, recombination_rate=2e-08,
                           mutation_rate=2e-08)
    return sim
