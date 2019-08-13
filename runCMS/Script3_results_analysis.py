# -*- coding: utf-8 -*-
"""
Created on Tue May 23 12:24:32 2017

@author: apolline
"""

import pandas as pd
import numpy as np
import sys

[pipeline_directory, analysis_directory, bim_file, summary_file, n_blocks, th] = sys.argv[1:7]

sys.path.append(pipeline_directory + '/others')  

from Script3_functions import compile_results,significant_snps,add_genes_ucsc,best_snp_per_region,count_significant_regions,result_plot,global_qqplot,global_quadrant_plot

blocks_directory  = analysis_directory + '/results_per_block'
summary_directory = analysis_directory + '/results_summary'
covariates_file   = analysis_directory + '/covariates.txt'

results_file         = summary_directory + '/complete_results.txt'
summary_snps_file    = summary_directory + '/summary_significant_snps.txt'
summary_regions_file = summary_directory + '/summary_significant_regions.txt'
counts_regions_file  = summary_directory + '/counts_regions.txt'

regions_file = pipeline_directory + '/others/regions.bed'

th = float(th)
res_list = ['_betaMA','_pvalMA','_betaMC','_pvalMC','_Yrsq','_Ncov']

# Get list of outcomes
phenotypes_summary = pd.read_csv(summary_file,sep=',')
pheno_list = phenotypes_summary.loc[np.where(phenotypes_summary['Outcome'] == 1)]['Label'] # list of phenotypes
print('Outcomes:\n')
print(pheno_list)

#################################################################################################################################
# 1) Compile results in one file
#################################################################################################################################
compile_results(blocks_directory, summary_directory, int(n_blocks), pheno_list, res_list, bim_file, regions_file, results_file)
data = pd.read_csv(results_file, sep = '\t')

#################################################################################################################################
# 2) Summary files
#################################################################################################################################

# All significant associations (one line per Phenotype - SNP association)
results_snps = significant_snps(data, pheno_list, res_list, th)
results_snps = add_genes_ucsc(results_snps, summary_directory) # Use UCSC database to match SNPs with genes
results_snps.to_csv(summary_snps_file, sep = '\t', index = 0)

# Associations grouped by region, i.e. independent LD blocks (one line per Phenotype - Region association)
results_regions = best_snp_per_region(results_snps, pheno_list)
results_regions.to_csv(summary_regions_file, sep = '\t', index = 0)

# Count significant regions per phenotype and per test (STD or CMS)
counts_regions = count_significant_regions(results_regions, pheno_list, th)
counts_regions.to_csv(counts_regions_file, sep = '\t', index = 0)

#################################################################################################################################
# 3) Plots
#################################################################################################################################

# 4 plots per phenotype: 2 Manhattan plots (Standard + CMS), QQplot and quadrant plot
result_plot(data, pheno_list, th, summary_directory)

# Summary QQplot (all phenotypes in one plot)
global_qqplot(data, pheno_list, summary_directory)

# Summary quadrant (all phenotypes in one plot)
global_quadrant_plot(data, pheno_list, th, summary_directory)
