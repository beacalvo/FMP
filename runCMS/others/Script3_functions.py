# -*- coding: utf-8 -*-
"""
Created on Mon May 29 09:50:47 2017

@author: apolline
"""
import pandas as pd
import numpy as np
import re
import os
import os.path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec


#################################################################################################################################
# Function to cut last two characters of a string
# INPUT: string  e.g. 'rs1234_A'
# OUTPUT: string e.g. 'rs1234'
################################################################################################################################# 
def cut_last(snp):
    snp = snp[0:(len(snp) - 2)]
    return snp   

#################################################################################################################################
# Function to tranform a list into a string
# INPUT: list of strings e.g. ['a', 'b', 'c']
# OUTPUT: string         e.g. 'a, b, c'
#################################################################################################################################
def list_to_string(l):
    if (type(l) == np.ndarray):
        s = l[0]
        l = l[1::]
        while (len(l) > 0):
            s = s + ', ' + l[0]
            l = l[1::]
        return(s)

#################################################################################################################################
# Function to add genome region corresponding to chromose and position (1703 independent LD blocks)
# INPUT: dataframe, file with regions
# OUTPUT: None
#################################################################################################################################  
def fill_regions(df, regions_file):
    regions = pd.read_csv(regions_file, sep='\t')
    reg = pd.DataFrame(np.zeros(df.shape[0], dtype = int))
    for r in range(0, 1784):
        left = regions[' start '][r]    
        right = regions[' stop'][r]    
        chrnum = int(''.join(re.findall('\d+', regions['chr '][r])))        
        ind = np.where((df['CHR'] == chrnum) & (df['POS'] >= left ) & (df['POS'] <= right ))
        #print(left, right, chrnum, ind)
        reg.iloc[ind] = r + 1
    df['Region'] = reg

#################################################################################################################################
# Function to compile results in one file with one line per SNP
# INPUT: Folder with results per blocks, folder for summarized results, number of blocks, list of phenotypes, list of results, bim file, regions file, file for summarized results
# OUTPUT: None
#################################################################################################################################
def compile_results(blocks_directory, summary_directory, n_blocks, pheno_list, res_list, bim_file, regions_file, results_file):

    summary_file = summary_directory + '/results_all_SNPs.txt'    
        
    os.system('for i in `seq 1 '+ str(n_blocks) + '`; do cat ' + blocks_directory + '/results_${i}.txt; done >> ' + summary_file)
    
    results = pd.read_csv(summary_file, sep = '\t', header = None)
    results.columns = ['SNP'] + [pheno + res for pheno in pheno_list for res in res_list]
    
    # Replace rsxx_A/T/C/G with rsxx
    results['SNP'] = results['SNP'].apply(cut_last)
    
    # Merge with bim file to get chromosome and position
    bim = pd.read_csv(bim_file, sep = '\t', header = None)[[0, 1, 3, 4, 5, 6]]
    bim.columns = ['CHR', 'SNP', 'POS', 'A1', 'A2', 'MAF']
    data = pd.merge(results, bim, on = 'SNP', how = 'left')
    
    # Compute regions from chromosome and position
    fill_regions(data, regions_file)
    
    # Export 
    data.to_csv(results_file, sep = '\t', index = 0)
    os.system('rm ' + summary_file)

#################################################################################################################################
# Function to build a data frame with significant associations (one line per Phenotype - SNP association)
# INPUT: data frame with all results, list of phenotypes, list of results, significance threshold
# OUTPUT: data frame with one line per significant Phenotype - SNP association
#################################################################################################################################
def significant_snps(df, pheno_list, res_list, th):
    col = ['Phenotype', 'Region', 'CHR', 'POS','SNP', 'A1', 'A2', 'MAF', 'Standard_beta', 'Standard_pval', 'CMS_beta', 'CMS_pval', 'CMS_rsquared', 'CMS_ncovs']
    table_signals = pd.DataFrame(columns = col, index = range(len(df)))
    
    s = 0
    
    print('\nBuild file with significant SNP-Phenotype associations...')
    for pheno in pheno_list:
        
        temp_cms = df.iloc[np.where(df[pheno + '_pvalMC'] < th)]
        temp_std = df.iloc[np.where(df[pheno + '_pvalMA'] < th)] 
        temp_all = pd.merge(temp_cms, temp_std, how = 'outer')
        
        temp_all = temp_all[['Region','CHR', 'POS','SNP', 'A1', 'A2', 'MAF'] + [pheno + res for res in res_list]]
        temp_all.columns = col[1::]
        temp_all.index = range(s, s + len(temp_all))
        
        for i in range(s, s + len(temp_all)):
            table_signals.loc[i, 'Phenotype'] = pheno
            table_signals.loc[i, col[1::]] = temp_all.loc[i]
                     
        s = s + len(temp_all)
        
    table_signals = table_signals.dropna(how = 'all')
    return(table_signals)

#################################################################################################################################
# Function to get genes from UCSC data base
# INPUT: data frame with significant Phenotype - SNP associations, folder for summarized results
# OUTPUT: data frame with genes in a new column
#
# Tables :
#  - refFlat: A gene prediction with additional geneName field
#     *geneName = Gene name
#     *txStart = Transcription start position
#     *txEnd = Transcription end position
#
#  - snp150: Polymorphism data from dbSNP
#     *name = SNP id
#     *chromStart =  	Start position in chromosome
#     *chromEnd = End position in chromosome
#################################################################################################################################
def add_genes_ucsc(results_snps, summary_directory):
    
    tmp_file = summary_directory + '/tmp_ucsc'
    
    list_snps = np.unique(results_snps['SNP']) # list of significant snps
    print(list_snps)
    
    command = 'mysql --user=genome --host=genome-euro-mysql.soe.ucsc.edu -A -D hg19 -e \'select R.geneName, R.txStart, R.txEnd, S.name, S.chromStart, S.chromEnd from snp150 as S left join refFlat as R on (S.chrom=R.chrom and not(R.txEnd+100000<S.chromStart or S.chromEnd+100000<R.txStart)) where S.name in ('
    
    for snp in list_snps:
        command = command + '"' + snp + '",'
    
    command = command + '"' + list_snps[len(list_snps)-1] + '")\' > ' + tmp_file
    print(command)
    os.system(command)
    
    ucsc = pd.read_csv(tmp_file, sep = '\t')
    
    ucsc = ucsc.dropna()
    ucsc = ucsc.drop_duplicates()
    ucsc['dist'] = range(len(ucsc))
    
    for i in ucsc.index:
        gene_start = ucsc.loc[i, 'txStart']
        gene_end   = ucsc.loc[i, 'txEnd']
        snp_start  = ucsc.loc[i, 'chromStart']
        snp_end    = ucsc.loc[i, 'chromStart']
        if(snp_start > gene_end): #gene before SNP
            dist = snp_start - gene_end
        elif (snp_end < gene_start): # gene after SNP
            dist = gene_start - snp_end
        else: # SNP on gene
            dist = 0
        ucsc.loc[i, 'dist'] = dist
    
    snps = np.unique(ucsc['name'])
    genes = pd.DataFrame(columns = ['SNP', 'Gene', 'dist'], index = range(len(ucsc)))
    
    l = 0
    for s in snps:
        
        tmp = ucsc.iloc[np.where(ucsc['name'] == s)][['name', 'geneName', 'dist']]
        min_dist = min(tmp['dist'])
        
        closest = tmp.iloc[np.where(tmp['dist'] == min_dist)]
        n = len(closest)
        closest.index = range(l,l+n)
        closest.columns = ['SNP', 'Gene', 'dist']
        genes.loc[range(l,l+n)] = closest
        
        l = l+n
    
    genes = genes.dropna()
    genes = genes.drop_duplicates()
    
    recap = pd.DataFrame()
    recap['Genes'] = genes.groupby('SNP')['Gene'].apply(list).apply(np.unique)
    recap['SNP'] = recap.index
            
    results_snps = pd.merge(results_snps, recap, on = 'SNP', how = 'left')
    results_snps['Genes'] = results_snps['Genes'].apply(list_to_string)
    
    os.system('rm ' + tmp_file)
    return(results_snps)

#################################################################################################################################
# Function to build a data frame with best snps per region (one line per Phenotype - Region association)
# INPUT: data frame with significant Phenotype - SNP associations, list of phenotypes
# OUTPUT: data frame with one line per Phenotype - Region association
#################################################################################################################################
def best_snp_per_region(results_snps, pheno_list):
    table_recap   = pd.DataFrame(columns = results_snps.columns, index = range(len(results_snps)))
    
    r = 0
    
    print('\nBuild file with significant Loci-Phenotype associations...')
    for pheno in pheno_list:
        
        temp = results_snps.iloc[np.where(results_snps['Phenotype'] == pheno)]       
        list_regions = np.unique(temp['Region'])
        
        for i in range(len(list_regions)):
            reg = list_regions[i]
            best_pval = min(temp[['Region', 'Standard_pval']].groupby('Region').min()['Standard_pval'][reg],temp[['Region', 'CMS_pval']].groupby('Region').min()['CMS_pval'][reg])
            if (len(temp.iloc[np.where(temp['CMS_pval'] == best_pval)]) != 0):   
               line = temp.iloc[np.where(temp['CMS_pval'] == best_pval)]
            else:
                line = temp.iloc[np.where(temp['Standard_pval'] == best_pval)]
            line.index = range(r + i, r + i + len(line))
            table_recap.loc[r + i, ] = line.loc[r + i]
                
        r = r + len(list_regions)
        
    table_recap = table_recap.dropna(how = 'all')
    return(table_recap)
 
#################################################################################################################################
# Function to count significant regions per phenotype
# INPUT: data frame with significant Phenotype - Region associations, list of phenotypes, significance threshold
# OUTPUT: data frame with number of significant regions per phenotype and per test (STD or CMS)
#################################################################################################################################   
def count_significant_regions(results_regions, pheno_list, th):    
    table_counts = pd.DataFrame(columns = ['Phenotype', 'Standard_only', 'CMS_only', 'Standard_and_CMS', 'Total'], index = range(len(pheno_list)))
         
    p = 0
    
    print('\nCount significant Loci-Phenotype associations...')
    for pheno in pheno_list:
        
        temp = results_regions.iloc[np.where(results_regions['Phenotype'] == pheno)]
        
        std = temp.iloc[np.where(temp['Standard_pval'] < th)]
        cms = temp.iloc[np.where(temp['CMS_pval'] < th)]
        
        table_counts['Phenotype'][p] = pheno
        table_counts['Standard_only'][p] = len(std) - len(set(std.index) & set(cms.index))
        table_counts['CMS_only'][p] = len(cms) - len(set(std.index) & set(cms.index))
        table_counts['Standard_and_CMS'][p] = len(set(std.index) & set(cms.index))
        table_counts['Total'][p] = len(set(std.index) | set(cms.index))
        
        p = p + 1
    
    return(table_counts)
    
#################################################################################################################################
# Function to build a Manhattan plot
# INPUT: data frame with all results, name of phenotype, name of p-value column, significance threshold, colors, folder for summarized results
# OUTPUT: None
#################################################################################################################################
def manhattan_plot(df, pheno, pval, th, colors, directory):
    threshold = -np.log10(th)
    
    df['-log10(' + pval + ')'] = - np.log10(df[pval])
    df_grouped = df.groupby(('CHR'))
    
    fig = plt.figure(figsize = (18, 9))
    ax = fig.add_subplot(111)
    x_labels = []
    x_labels_pos = []
  
    m = 0
    i = 0
    
    for num, (name, group) in enumerate(df_grouped):
        group = pd.DataFrame(group)
        group['absolute_position'] = group['POS'] + m
        m = np.max(group['absolute_position'])
        group.plot(kind = 'scatter' , x = 'absolute_position', y = '-log10('+pval+')', color = colors[i], ax = ax)
        x_labels.append(name)
        x_labels_pos.append((group['absolute_position'].iloc[-1] - (group['absolute_position'].iloc[-1] - group['absolute_position'].iloc[0])/2))
        i = i + 1
        
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlim([0, m])
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('-log10(P)')
    plt.axhline(y = threshold, color = 'r', linestyle = '-')
    if (pval[len(pval)-1] == 'A'):
        plt.title('Standard test', size = 16)
        fig.savefig(directory + '/Manhattan_plot_' + pheno + '_Standard.png', dpi=120, bbox_inches='tight')
    else:
        plt.title('CMS test', size = 16)
        fig.savefig(directory + '/Manhattan_plot_' + pheno + '_CMS.png', dpi=120, bbox_inches='tight')   
    
    plt.close(fig)


#################################################################################################################################
# Function to build a QQplot
# INPUT: data frame with all results, name of phenotype, folder for summarized results
# OUTPUT: None
#################################################################################################################################
def qqplot(df, pheno, directory):
    n_snps = len(df)
    
    pvalue_std = -np.log10(df[pheno+'_pvalMA'])
    pvalue_cms = -np.log10(df[pheno+'_pvalMC'])
    
    x_unif = -np.log10(np.arange(1, n_snps + 1)/(n_snps + 2)) #Uniformed values between 0 and 1, 0 and 1 excluded
    y_std = np.sort(pvalue_std[:])
    y_cms = np.sort(pvalue_cms[:])
    
    fig = plt.figure(figsize = (9, 9))
    plt.scatter(x_unif[::-1], y_std, s = 5, c = 'blue')
    plt.scatter(x_unif[::-1], y_cms, s = 5, c = 'red')
    
    x = np.linspace(0, 6)
    plt.plot(x, x, c = 'black')
    
    lambda_value_std = np.nanmedian(y_std) / np.nanmedian(x_unif)
    lambda_value_cms = np.nanmedian(y_cms) / np.nanmedian(x_unif)
    
    plt.text(x = 0, y = 8, s = 'Standard test: lambda = ' + str(lambda_value_std), color = 'blue')
    plt.text(x = 0, y = 7, s = 'CMS test: lambda = '      + str(lambda_value_cms), color = 'red')
    plt.title('QQplot', size = 16)
    
    fig.savefig(directory + '/QQplot_' + pheno + '.png', dpi = 120, bbox_inches = 'tight')
    plt.close(fig)


#################################################################################################################################
# Function to build a Quadrant plot (distinction of STD or CMS significant regions)
# INPUT: data frame with all results, name of phenotype, significance threshold, folder for summarized results
# OUTPUT: None
#################################################################################################################################        
def quadrant_plot(df, pheno, th, directory):
    
    fig = plt.figure(figsize=(9, 9))
    
    n_stand = 0
    n_cms = 0
    n_both = 0
    n_others = 0
               
    std = df[[pheno + '_pvalMA', 'Region']].groupby(('Region')).min()[pheno + '_pvalMA']
    cms = df[[pheno + '_pvalMC', 'Region']].groupby(('Region')).min()[pheno + '_pvalMC']      
    
    for i in range(1, 1785):
        if (i in std.index):
            if(cms[i] < th and std[i] > th):
                plt.scatter(-np.log10(std[i]), -np.log10(cms[i]), c = 'lightgreen')
                n_cms = n_cms + 1
            elif(cms[i] < th and std[i] < th):
                plt.scatter(-np.log10(std[i]), -np.log10(cms[i]), c = 'turquoise')
                n_both = n_both + 1
            elif(cms[i] > th and std[i] < th):
                plt.scatter(-np.log10(std[i]), -np.log10(cms[i]), c = 'salmon')
                n_stand = n_stand + 1
            else:
                plt.scatter(-np.log10(std[i]), -np.log10(cms[i]), c = 'grey')
                n_others = n_others + 1
            
    plt.axvline(-np.log10(th), color = 'grey', linestyle = '--')
    plt.axhline(-np.log10(th), color = 'grey', linestyle = '--')
    
    plt.xlabel('-log10(P) for standard test', fontsize=18)
    plt.ylabel('-log10(P) for CMS test',      fontsize=16)
    
    red_patch   = mpatches.Patch(color='salmon',     label='Significant pvalues for standard test only ('     + str(n_stand) + ')')
    blue_patch  = mpatches.Patch(color='turquoise',  label='Significant pvalues for CMS and standard tests (' + str(n_both) + ')')
    green_patch = mpatches.Patch(color='lightgreen', label='Significant pvalues for CMS test only ('          + str(n_cms)+')')
    black_patch = mpatches.Patch(color='grey',       label='Non significant pvalues ('                        + str(n_others)+')')
    
    plt.legend(handles = [red_patch, blue_patch, green_patch, black_patch])
    plt.title('CMS vs Standard test', size = 16)
    
    fig.savefig(directory + '/Quadrant_plot_' + pheno + '.png', dpi = 120, bbox_inches = 'tight')
    plt.close(fig)
    

#################################################################################################################################
# Function to gather two Manhattan plot (STD + CMS), a QQplot and a quadrant plot on the same image for one phenotype
# INPUT: data frame with all results, list of phenotypes, significance threshold, folder for summarized results
# OUTPUT: None
#################################################################################################################################
def result_plot(df, pheno_list, th, directory):
    
    colors = [[ 0.23625862,  0.52860261,  0.78364119], [ 0.03608782,  0.55719131,  0.1664441 ], [ 0.94307407,  0.35026836,  0.68064485], [ 0.25533632,  0.18023988,  0.42811556], [ 0.83003783,  0.81582078,  0.37646543], [ 0.65316214,  0.13796559,  0.45594953], [ 0.47797385,  0.18014566,  0.09708977], [ 0.44313678,  0.59353625,  0.50071228], [ 0.5985085 ,  0.76463359,  0.75645524], [ 0.67242609,  0.13103319,  0.17105948], [ 0.05590727,  0.56748858,  0.47916381], [ 0.80333682,  0.72828047,  0.39257576], [ 0.46142155,  0.65768434,  0.98160204], [ 0.85915425,  0.80420234,  0.47793356], [ 0.16912189,  0.17209861,  0.4464207 ], [ 0.6452403 ,  0.39523205,  0.59043336], [ 0.56400114,  0.46104143,  0.71735811], [ 0.86413373,  0.26718486,  0.04390139], [ 0.76195355,  0.17223056,  0.38066256], [ 0.4721071 ,  0.18925624,  0.25739407], [ 0.62012875,  0.97591735,  0.21921685], [ 0.64509652,  0.54463331,  0.23358183], [0.035, 0.556, 0.541], [0.556, 0.392, 0.035]]
    
    print('\nBuild plots for each phenotype (Manhattan plots + QQplot + Quadrant plot)...')
    for pheno in pheno_list:
        print(pheno)
        manhattan_plot(df, pheno, pheno + '_pvalMA', th, colors, directory)
        manhattan_plot(df, pheno, pheno + '_pvalMC', th, colors, directory)
        qqplot(df, pheno, directory)
        quadrant_plot(df, pheno, th, directory)
    
        manMA = mpimg.imread(directory + '/Manhattan_plot_' + pheno + '_Standard.png')
        manMC = mpimg.imread(directory + '/Manhattan_plot_' + pheno + '_CMS.png')
        qq =    mpimg.imread(directory + '/QQplot_' + pheno + '.png') 
        quad =  mpimg.imread(directory + '/Quadrant_plot_' + pheno + '.png')
        
        fig = plt.figure(figsize=(27, 18))
        
        gs = gridspec.GridSpec(2, 2, width_ratios=[2, 1], height_ratios=[1, 1])
                       
        # First Manhattan plot (standard test)
        plt.subplot(gs[0])
        plt.imshow(manMA)
        plt.axis('off')
        
        # QQ plot
        plt.subplot(gs[1])
        plt.imshow(qq)
        plt.axis('off')
        
        # Second Manhattan plot (CMS test)
        plt.subplot(gs[2])
        plt.imshow(manMC)
        plt.axis('off')
        
        # Quadrant plot
        plt.subplot(gs[3])
        plt.imshow(quad)
        plt.axis('off')
    
        #plt.suptitle(pheno, size = 32)
        
        plt.tight_layout()
        
        fig.savefig(directory + '/Plots_' + pheno + '.png', dpi = 240,  bbox_inches = 'tight')
        plt.close(fig)
        
        os.system('rm ' + directory + '/Manhattan_plot_*')
        os.system('rm ' + directory + '/QQplot_*')
        os.system('rm ' + directory + '/Quadrant_plot_*')


#################################################################################################################################
# Function to build a QQplot with all phenotypes together
# INPUT: data frame with all results, list of phenotypes, folder for summarized results
# OUTPUT: None
#################################################################################################################################
def global_qqplot(df, pheno_list, directory):
    
    fig = plt.figure(figsize=(9, 9)) 
    
    n_points = len(df) * len(pheno_list)
    x_unif = -np.log10(np.arange(1, n_points + 1) / (n_points + 2))[::-1] #Uniformed values between 0 and 1, 0 and 1 excluded
    y_std = []
    y_cms = []
    
    print('\nBuild global QQplot with all phenotypes...')
    for pheno in pheno_list:        
        y_std = y_std + list(-np.log10(list(df[pheno + '_pvalMA'])))
        y_cms = y_cms + list(-np.log10(list(df[pheno + '_pvalMC'])))  
    
    y_std = np.sort(y_std)
    y_cms = np.sort(y_cms)
    
    plt.scatter(x_unif, y_std, s = 5, c = 'blue')
    plt.scatter(x_unif, y_cms, s = 5, c = 'red')
        
    x = np.linspace(0, 8)
    plt.plot(x, x, c = 'black')
        
    lambda_value_std = np.nanmedian(y_std) / np.nanmedian(x_unif)
    lambda_value_cms = np.nanmedian(y_cms) / np.nanmedian(x_unif)
        
    plt.text(x = 0, y = 6,  s = 'Standard test: lambda = ' + str(lambda_value_std), color = 'blue')
    plt.text(x = 0, y = 10, s = 'CMS test: lambda = '      + str(lambda_value_cms), color = 'red')
    plt.title('QQplot', size = 16)
        
    fig.savefig(directory + '/Global_QQplot.png', dpi = 120, bbox_inches = 'tight')
    plt.close(fig)
        
    fig = plt.figure(figsize = (9, 9))
    plt.scatter(x_unif[range(int(0.9999 * n_points))], y_std[range(int(0.9999 * n_points))], s = 5, c = 'blue')
    plt.scatter(x_unif[range(int(0.9999 * n_points))], y_cms[range(int(0.9999 * n_points))], s = 5, c = 'red')
    
    x = np.linspace(0,2)
    plt.plot(x, x, c = 'black')
    fig.savefig(directory + '/Global_QQplot_zoom.png', dpi = 120, bbox_inches = 'tight')
    plt.close(fig)

       
#################################################################################################################################
# Function to build a Quandrant plot with all phenotypes together
# INPUT: data frame with all results, list of phenotypes, significance threshold, folder for summarized results
# OUTPUT: None
#################################################################################################################################        
def global_quadrant_plot(df, pheno_list, th, directory):
    
    fig = plt.figure(figsize=(9, 9))
    
    n_stand = 0
    n_cms = 0
    n_both = 0
    n_others = 0
    
    print('\nBuild global Quadrant plot with all phenotypes...')
    for pheno in pheno_list:
        std = df[[pheno + '_pvalMA', 'Region']].groupby(('Region')).min()[pheno + '_pvalMA']
        cms = df[[pheno + '_pvalMC', 'Region']].groupby(('Region')).min()[pheno + '_pvalMC']       
        
        for i in range(1, 1785):
            if (i in std.index):
                if(cms[i] < th and std[i] > th):
                    plt.scatter(-np.log10(std[i]), -np.log10(cms[i]), c = 'lightgreen')
                    n_cms = n_cms + 1
                elif(cms[i] < th and std[i] < th):
                    plt.scatter(-np.log10(std[i]), -np.log10(cms[i]), c = 'turquoise')
                    n_both = n_both + 1
                elif(cms[i] > th and std[i] < th):
                    plt.scatter(-np.log10(std[i]), -np.log10(cms[i]), c = 'salmon')
                    n_stand = n_stand + 1
                else:
                    plt.scatter(-np.log10(std[i]), -np.log10(cms[i]), c = 'grey')
                    n_others = n_others + 1
            
    plt.axvline(-np.log10(th), color = 'grey', linestyle = '--')
    plt.axhline(-np.log10(th), color = 'grey', linestyle = '--')
    
    plt.xlabel('-log10(P) for standard test', fontsize = 18)
    plt.ylabel('-log10(P) for CMS test'     , fontsize=16)
    
    red_patch   = mpatches.Patch(color = 'salmon',     label='Significant pvalues for standard test only ('     + str(n_stand) + ')')
    blue_patch  = mpatches.Patch(color = 'turquoise',  label='Significant pvalues for CMS and standard tests (' + str(n_both) + ')')
    green_patch = mpatches.Patch(color = 'lightgreen', label='Significant pvalues for CMS test only ('          + str(n_cms) + ')')
    black_patch = mpatches.Patch(color = 'grey',       label='Non significant pvalues ('                        + str(n_others) + ')')
    
    plt.legend(handles = [red_patch, blue_patch, green_patch, black_patch])
    plt.title('CMS vs Standard test', size = 16)
    
    fig.savefig(directory + '/Global_Quadrant_plot.png', dpi = 120, bbox_inches = 'tight')
    plt.close(fig)  
