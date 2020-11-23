# This script details different msprime models, which are saved as functions. 
# These can then be executed in other files with greater ease.
# There should be a short description above each history.
# Times should be given in generations

import math
import msprime
import sys
import pdb


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
                           demographic_events=demographic_events, length=3e3, recombination_rate=2e-8,mutation_rate=2e-8,random_seed=24)
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
