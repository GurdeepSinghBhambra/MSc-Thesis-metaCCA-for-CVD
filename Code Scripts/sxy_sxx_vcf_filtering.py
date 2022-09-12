# SXY SXX finalizing
# Using dask python to handle large SXX file

import dask.dataframe as dd
import pandas as pd

import datetime


st = datetime.datetime.now()

# Getting Pruned SNP list
# getting pruned snps to be used to filter sxx and sxy
pruned_snp_list = pd.read_csv('./data/g1000p3_EUR/prune/g1000p3_EUR_common_GTs_sorted_clean_head_ld_bed_1000k_02.prune.in', 
							  delimiter='\t', header=None)

# the file only has snps in the first column
pruned_snp_list = pruned_snp_list[0].values # storing the numpy array object instead of series object

print("Pruned SNP list loaded with size: ", pruned_snp_list.shape, "\n")
print("time taken:",datetime.datetime.now()-st)

# Pruning SXY_study_maf1.txt

st = datetime.datetime.now()
sxy_study1 = pd.read_csv('./data/processed/sxy_study_1_maf1.txt', delimiter='\t')

print(sxy_study1.shape)

print(sxy_study1.head())

filtered_sxy = sxy_study1[sxy_study1['SNPID_UKB'].isin(pruned_snp_list)]

print("SNPs after filtering:", filtered_sxy.shape[0], "/", pruned_snp_list.shape[0])

filtered_sxy.to_csv('./data/processed/sxy_study_1_maf1_pruned.txt', sep='\t', index=False)
print("SXY_STUDY_1 pruned file saved at:", './data/processed/sxy_study_1_maf1_pruned.txt')
print("time taken:",datetime.datetime.now()-st)
del filtered_sxy, sxy_study1 # releasing memory

print("\n\n")


st = datetime.datetime.now()
# filtering maf1_data.vcf file as it will used for gene annotation
vcf = pd.read_csv('./data/processed/maf1_data.vcf', delimiter='\t', header=None)
# columns=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'NA1', 'NA2', 'NA3']

print(vcf.shape)

print(vcf.head())

filtered_vcf = vcf[vcf[2].isin(pruned_snp_list)]

print("SNPs after filtering:", filtered_vcf.shape[0], "/", pruned_snp_list.shape[0])

filtered_vcf.to_csv('./data/processed/maf1_data_pruned.vcf', sep='\t', index=False, header=False)
print("VCF pruned file saved at:", './data/processed/maf1_data_pruned.vcf')
print("time taken:",datetime.datetime.now()-st)
del filtered_vcf, vcf # releasing memory

print("\n\n")


## Filtering GT file
st = datetime.datetime.now()
# loading the large GT file using dask as it uses batch and parallel processing
gt_df = dd.read_csv('./data/g1000p3_EUR/g1000p3_EUR_common_GTs_sorted_nohead.txt', #filepath
					  blocksize='500MB', # allocating about 500MB of memory for partition size to be read in RAM for parallel execution, if facing any memory errors, comment this line or reduce size
					  delimiter='\t', header=None) # since file has no head and are tab seperated, these options are used

print("GT file ready, filtering and saving GT pruned file ...")

gt_pruned_filepath = './data/g1000p3_EUR/g1000p3_EUR_common_GTs_sorted_nohead_ld02_pruned.txt'

gt_df = gt_df[gt_df[2].isin(pruned_snp_list)].compute()

print("SNPs after filtering:", gt_df.shape[0], "/", pruned_snp_list.shape[0])

gt_df.to_csv(gt_pruned_filepath, sep='\t', index=False, header=False)

print("GT pruned file saved at:", gt_pruned_filepath)
print("time taken:",datetime.datetime.now()-st)
