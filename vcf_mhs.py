# script hosts a function which writes a vcf and from this writes a multihetsep, from a given msprime simulation
# see http://evomics.org/learning/population-and-speciation-genomics/2020-population-and-speciation-genomics/demography-psmc-msmc/

# algorithm:
# 1. simulate history with msprime and save this as sim
# 2. save array of variants (necessary?)
# 3. write a vcf to disc
# 4. save multihetsep from this vcv iterate along 

import allel
import pandas as pd
import os
from datetime import datetime
import pdb
import time

def dir_check(vcf_path,mhs_path):
    print('checking directories')

    if not os.path.isdir(vcf_path): 
        print('{} does not exist. Creating it now.'.format(vcf_path))
        os.mkdir(vcf_path)
    if not os.path.isdir(mhs_path):
        print('{} does not exist. Creating it now.'.format(mhs_path))
        os.mkdir(mhs_path)
    
    return None


def sim_to_mhs(sim,vcf_path=os.getcwd() + '/vcf_mhs/',vcf_filename=datetime.now().strftime("%Y%m%d"),mhs_path = os.getcwd() + '/vcf_mhs/',mhs_filename = datetime.now().strftime("%Y%m%d") + '_mhs',suffix = '',verbose="True" ):
    # function: given an msprime simulation, write a vcf for it and then call multihetsep() to write mhs
    # sim: msprime simulation object
    # vcf_path: path where you want to save vcf. Default is current directory
    # vcf_filename: the name of the vcf file that you want to write. Default is current date and time
    # vcf_filename = vcf_name
    # mhs_path: where do you want to save this mhs
    # mhs_filename: where do you want to save this file
    # verbose: Whether to print progress or not

    # get array of genotype for each variant. TODO this isn't necessary for PSMC's purposes
    # vars = [varient.genotypes for varient in sim.variants()]

    # check given directories are valid
    dir_check(vcf_path ,mhs_path)

    #write vcf of genotypes
    if verbose: print('Writing vcf...')
    with open(vcf_path + vcf_filename + suffix + ".vcf", "w") as vcf_file:
        sim.write_vcf(vcf_file, 2)
    if verbose: print('vcf written to {}'.format(vcf_path))

    # write matrix of genotypes, maybe delete this
    gen_mat = sim.genotype_matrix()

    #read vcf as data frame
    sim_vcf = allel.vcf_to_dataframe(vcf_path + vcf_filename + suffix + ".vcf",fields='*')

    # generate and write mhs
    multihetsep(sim_vcf,mhs_path,mhs_filename,suffix,gen_mat,verbose)

    # return the vcf file as a dataframe
    return None

def multihetsep(sim_vcf,mhs_path ,mhs_filename,suffix, gen_mat,verbose):
    # function: writes a multihetsep file
    # sim_vcf: the dataframe of the vcf


    t0 = time.time()
    
    # initialise data frames
    index = range(0,len(sim_vcf))
    columns = ['CHROM','POS','SSPSS','HAP']
    data = pd.DataFrame(index = index,columns = columns)
    
    # array for location of where NA's occur
    na_location = []

    # iterate along every row of the vcf
    if verbose: print('generating mhs...')
    for i in range(len(sim_vcf)):
        if i == 0: # for the start of the file
            row = ['chr' + str(sim_vcf['CHROM'][i]),int(sim_vcf['POS'][i]),int(sim_vcf['POS'][i]),str(gen_mat[i][0]) + str(gen_mat[i][1])]
            data.loc[i] = row
        elif i > 0:
            #sometimes the VCF files duplicates a row. And hence the resulting multihetsep file is invalid (it has 0 in 3rd column). These lines ensure skipping over that
            SSPSS = sim_vcf['POS'][i] - sim_vcf['POS'][i - 1] # caluclate the sites sinces previous segregating site (SSPSS)
            if SSPSS ==0:
                na_location.append(i)
                continue
            row = ['chr' + str(sim_vcf['CHROM'][i]),int(sim_vcf['POS'][i]),int(SSPSS),str(gen_mat[i][0]) + str(gen_mat[i][1])]
            data.loc[i] = row
        
    if verbose: print('mhs generated.')

    # remove NA rows
    data = data.drop(na_location)
    
    # convert to ints
    # data['POS'] = data['POS'].astype(int)
    # data['SSPSS'] = data['SSPSS'].astype(int)

    #write as csv
    if verbose: print('writing mhs to disc...')
    data = data.to_csv(mhs_path + mhs_filename + suffix + '.txt',sep='\t',header=False,index=False,line_terminator='\n')
    if verbose: print('mhs written to {}'.format(mhs_path))

    t1 = time.time()
    total = t1-t0   
    print('total time taken {}'.format(total))
    return None