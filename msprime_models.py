# This script details different msprime models, which are saved as functions. 
# These can then be executed in other files with greater ease.
# There should be a short description above each history.
# Times should be given in generations

import math
import msprime
import sys


# simple population history, constant population size with no growth or migration
# one diploid individual
def m0001(print):
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
    if print:
        dd.print_history()
    sim = msprime.simulate(population_configurations=population_configurations,
                           migration_matrix=migration_matrix,
                           demographic_events=demographic_events, length=6e7, recombination_rate=2e-8,mutation_rate=2e-8,random_seed=24)
    return sim


# delete this
# model with instantaneous migration 
def model_struct6(print):
    N_A0 = 1e+04
    N_B0 = 1e+04
    T_A = 2e+04
    T_B = 8e+04
    # Times are provided in years, so we convert into generations.
    generation_time = 25
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
            time=T_A, source = 0, dest=1,proportion=0.5),
        msprime.MassMigration(
            time=T_B, source=1, dest=0, proportion=1)
    ]
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    dd = msprime.DemographyDebugger(
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events)
    #dd.print_history()
    sim = msprime.simulate(population_configurations=population_configurations,
                           migration_matrix=migration_matrix,
                           demographic_events=demographic_events, length=3e+7, recombination_rate=2e-8,
                           mutation_rate=2e-8,random_seed=50)
    return sim


