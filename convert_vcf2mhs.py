# make mhs from given vcf

import allel 
from vcf_mhs import multihetsep
import argparse
import os
import pdb
import sys
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-dir",help="path to dir of vcf that you want converted into mhs files",default=os.getcwd() +'/vcf_mhs/',type=str)
parser.add_argument("-file","--filename",help="path to vcf you want converted into a mhs",type=str)
args = parser.parse_args()

filename = args.dir + args.filename

# check is file
print('File exists? {}'.format(os.path.isfile(filename)))
if args.filename[-4:len(args.filename)] == '.vcf':
    print('File is {}'.format(filename))
    print('File is vcf.')
else: 
    print('File is {}'.format(filename))
    print('File is not vcf. Aborting.')
    sys.exit()

# load vcf
sim_vcf = allel.vcf_to_dataframe(filename,fields='*')

# read vcf with tsk column, for genotype matrix
# TODO quite sure this doesn't really matter. Could just put arbitrary for genotype column
columns = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','tsk_0']
vcf_tsk = pd.read_csv('/home/trevor/ancestry/vcf_mhs/20201129_matchingPSC_mig1_142357.vcf',names = columns,delimiter="\t",comment = '#')
genotypes = vcf_tsk['tsk_0']
gen_mat = [gen[0] + gen[2] for gen in genotypes]
mhs_filename = args.filename[0:args.filename.find('_')] + '_mhs_' + args.filename[(args.filename.find('_')+1):(len(args.filename)-4)] 
suffix = ''
verbose=True
multihetsep(sim_vcf,args.dir,mhs_filename,suffix,gen_mat,verbose)