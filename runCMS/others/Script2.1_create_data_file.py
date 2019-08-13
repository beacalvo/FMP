
# -*- coding: utf-8 -*-
"""
Created on Wed May 17 17:20:09 2017

@author: apolline
"""

import sys
import os
from rpy2 import robjects
r = robjects.r

###############################################################################
# Get the arguments from the batch command
[pipeline_directory , analysis_directory, genotypes_file, summaryFile, i, SNP1, SNP2] = sys.argv[1:8]

phenotypes_file   = analysis_directory + '/phenotypes_for_analysis.txt'
covariates_file   = analysis_directory + '/covariates.txt'
temp_directory    = analysis_directory + '/tmp'
results_directory = analysis_directory + '/results_per_block'
log_file          = analysis_directory + '/log_Script2.txt'

# Import R script to run CMS
r.source(pipeline_directory + '/others/Script2.2_run_CMS.R')

# Define temporary data file and result file for this block
temp_genotypes_file = temp_directory + '/temp_genotypes_' + i + '.txt'
results_file = results_directory+'/results_' + i + '.txt'


###############################################################################
# Create temporary data file by cutting the block from the genotype file
k = int(SNP1) # start at SNP1
command = 'cut -f2' # get second column (ID)

while(k < (int(SNP2) + 1)): # end at SNP2
    command = command + ',' + str(k + 6) # get each SNP column between SNP1 and SNP2
    k = k + 1

command = command + ' -d " " ' + genotypes_file + ' > ' + temp_genotypes_file # create the temporary file
os.system(command)

###############################################################################
# Launch the R script to run CMS on the created temporary file
r.CMS_analysis(pipeline_directory, phenotypes_file, summaryFile, covariates_file, temp_genotypes_file, results_file)

os.remove(temp_genotypes_file)
