# Main code


The main code used throughout this project is available at this folder. The code can be in three different programming languages: Bash, Python and, mainly, R. The main code is divided in smaller chunks of code.


## 1. Code to summarize variant quality control results 



This script was used to count the number of variants excluded from the original dataset after the variant quality control, as well as to create some plots to compare the distribution of some parameters before and after the variant quality control. 



### 1.1. Counting excluded variants

On the one hand, this script indicates, for each chromosome, the number of variants removed owing to each quality control filter:

* Variant call rate < 95%
* Variants not in HWE (P<sub>HWE</sub> < 1.0x10<sup>-06</sup>) 
* MAF < 1% 

Finally, the script summarizes the results per chromosome into a single summary file containing the results of all chromosomes.

NOTE: This script was first aimed at analyzing imputed genotype data. For this reason, this script also outputs the number of variants falling below an imputation accuracy (R<sup>2</sup>) threshold. In fact, three different minimum R<sup>2</sup> thresholds were tested. The number of variants with an R<sup>2</sup> below 0.25, 0.50 or 0.80 were reported, in order to help establish a definitive R<sup>2</sup> cut-off.



### 1.2. Summarizing performance of QC by creating plots

On the other hand, this script creates plots that allow to visually compare the distribution of some parameters before and after variant quality control. In this case, R<sup>2</sup> and MAF before and after quality control are represented in histograms. Moreover, a plot correlating R<sup>2</sup> and MAF is also obtained.

(See scripts: 1.2.QC_summary_1.R, 1.2.QC_summary_2.sh and 1.2.QC_summary_3.R)





## 2. Code used for data processing




### 2.1. Metabolomic data processing

The metabolomic data was contained in two different datasets: one fore serum and one for urine. In both cases, metabolite levels were contained in ExpressionSet objects. The *Minfi* package was used to manage this data (1). This script was used to join both datasets and, once serum and urine metabolomic data had been merged, the resulting dataset was merged to the phenodata.

Moreover, metabolites names were changed to the reference names. Common metabolites in serum and urine were addressed by adding to each metabolite name the extension ".serum" and ".urine" to serum and urine metabolites, respectively.

Next, metabolite ratios were computed from the 44 urine metabolites.

(See script 2.1.Metabolomic_data_processing.R)

### 2.2. Metadata data processing

As mentioned in the methodology in Section 3.4., the metadata was used to identify those individuals for which the ancestry predicted by the *Peddy* program was not European. Then, these individuals were removed from the dataset. This way, the final dataset contained metadata, serum and urine metabolite levels, as well as urine ratios, of the 996 individuals with European predicted ancestry.

(See 2.2.Metadata_data_processing.R)

### 2.3. Genotype data processing

The genotype data was already in PLINK binary format (2). However, it contained variants located in the mitochondrial DNA and in the non-canonical PAR region of chromosome Y. These variants were excluded from further analysis. To do so, the following script was used.

(See 2.3.Genotype_data_processing.sh)

It was checked that genotyped variants in chromosome X were coded as 0 or 2 in males (representing a single copy of the X chromosome), whereas females variants in chromosome X were coded as 0, 1 or 2 (owing to the pair of X chromosomes). Instead, for the PAR regions (chromosome 25 in PLINK), both males and females variants were coded as diploid (0, 1 or 2).




## 3. Preparing input data for CMS method
 



The CMS scripts required the input data to be in a specific format, as well as the creation of other files with information of our data.





### 3.1. Code to prepare genotype data for statistical analysis










#### 3.1. Preparing genotype data for whole-genome analysis


*Script 2* required the genotype data to be recoded from the binary PLINK format to the `.raw` file, in order to use the `.raw` file as input. In this format genotypes are coded as a single allele dosage number (additive coding), which represents the number of minor alleles per person.

*Script 3* requires the `BIM` file as input, in order to include SNP-related information (chromosome, variant identifier, position, base-pair coordinate, allele 1 and allele 2) to the summary tables. In order for the MAF to be included in these tables, the `BIM` file was modified by adding a column containing the MAF for each SNP.

(See scripts 3.1.Preparing_genotype_data_1.sh and 3.1.Preparing_genotype_data_2.sh)



### 3.2. Code to prepare phenotype data for statistical analysis


The `phenotypes_file.txt` required as input in *Script 1*  was a text file, with columns separated by tabulations and a header line. It contained one line per individual and the first column was the individual ID. This file contained both confounding variables and serum and urine metabolite levels. The runCMS code did not allow categorical variables (i.e. sex and urine sampling type), so the categorical variables in the metadata were transformed to dummy variables. The following script shows how the sex and urine sampling type variables were transformed to dummy variables and, finally, the `phenotypes_file.txt` used as input in *Script 1* is created. The `fastDummies` (3) package was used.

(See 3.2.Code_prepare_phenotype_data.R)

### 3.3. Code to create other input files required to run CMS


The runCMS code required the creation of another file as input of *Scripts 1*, *2* and *3*, called `summary_file.csv`. This is a `.csv` file with columns separated by commas and a header line. This file aims at describing the role of each variable contained in the phenotypes file. For each selected variable, a label and a binary indicator for classification as confounding factors (i.e. variables systematically included as covariates), outcome (i.e. each single variable that will be treated as a primary outcome) and candidate covariates (i.e. variables that will be assessed by CMS for inclusion as a covariate) are provided. By default, all variables in "Covariates" column are included as covariates in each outcome analysis. However, the "Excluded" column can be used to exclude specific variables from being considered as covariates for a given outcome. In our study, we excluded the covariates that were hierarchical parent of the outcome under study or vice versa, in order to reduce bias. For example, "X3.aminoisobutyrate.urine" was excluded from being a possible covariate of "X3.hydroxybutyrate.3.aminoisobutyrate.urine" and vice versa, since "X3.hydroxybutyrate.3.aminoisobutyrate.urine" represents the concentrations of both "X3.aminoisobutyrate.urine" and "3.hydroxybutyrate" (i.e. the metabolite "X3.hydroxybutyrate.3.aminoisobutyrate.urine" is a hierarchical parent of "X3.aminoisobutyrate.urine", this is why both metabolites are excluded as possible covariates of the other metabolite).
Moreover, for the urine metabolites that were also measured in serum, the corresponding serum metabolite was added to the "Excluded" column, since there is a strong possibility that they share some genetic source of variability.

(See 3.3.Code_create_input_files_CMS.R)

## 4. Computing heritability

SNP heritability was computed using genome-wide complex trait analysis (GCTA), which estimates the additive contribution to a trait's heritability of a particular subset of SNPs (GREML analysis) (4). GCTA is a command-line tool with options similar to the commonly used PLINK software. Only autosomal SNPs and SNPs with MAF > 0.01 were considered to build the GRM.

First, the input files needed for the analysis were created (See 4.Computing_heritability_1.R).

Next, the GRM was computed by GCTA and GCTA-GREML analysis was carried out (See 4.Computing_heritability_2.sh, 4.Computing_heritability_3.R and 4.Computing_heritability_4.sh).

Finally, SNP-heritability results were summarized. Moreover, results were plotted in two circular plots. One represented variation from additive genetic effects (V(G)), variation from residual effects V(e) and the total phenotypic variation (Vp); while the other one showed the proportion of genotypic to phenotypic variation (V(G)/Vp). In order to create the plots, the R packages from CRAN *tidyverse* (5), *viridis* (6) and *ggplot2* (7) were used (See 4.Computing_heritability_5.R).


### 5.1. Modifications made to the CMS scripts


The runCMS pipeline offers the code used by the authors to run their analysis. Therefore, some modifications needed to be carried out in order to adapt the script to our analysis and the cluster used. 


* *Script 1* consists in an R code (`Script1_data_preparation.R`). This script was run without any modification.

* *Script 2* is a SLURM batch script (`Script2_launch_analysis_SLURM.sh`) used to launch a job array in a cluster, with one task per SNP block. This script was adapted to the cluster used by modifying the modules loading and the memory requested, as well as the output and error log directories. Once the job array was launched, each job called a Python script (`Script2.1_create_data_file.py`) which ran CMS on each SNP block, by calling an R script (`Script2.2_run_CMS.R`), and wrote results in one file per block, by calling another R script (`Script2.3_CMS.R`).

* *Script 3* is a Python script (`Script3_results_analysis.py`) which parses the results to build plots and summary tables. This script requires the `BIM` file as input, in order to include SNP-related information (chromosome, position, minor allele, etc.) to the summary tables. However, the MAF is not included in any of these tables. This is why the script was modified so that the MAF was reported in these summary tables. More concretely, the modification was done in file `Script3_functions.py`, in the `compile_results` and `significant_snps` functions. In addition, this script summarizes results by loci which correspond to approximately independent LD blocks based on a recombination map (`regions.bed`). However, this map only considers chromosomes 1 to 22. This is why it was modified so that it included chromosome X and PAR1 and PAR2:

    + Each pseudoautosomal region was considered as a locus, so two loci were added each one spanning the length of one pseudoautosomal region.
    + No results could be found in the literature regarding independent LD blocks or recombination pattern in chromosome X. For this reason, the X chromosome was divided randomly into regions of 2 Mb, which resulted in a total of 79 loci.

    After modification of the `regions.bed` file, *Script 3* was adapted to consider the newly incorporated loci. Modifications were made in file `Script3_functions.py`, in the `fill_regions`, `quadrant_plot` and `global_quadrant_plot` functions.


The resulting scripts and `regions.bed` file after modification were the ones used to carry out our analysis. These files, together with the rest of unmodified scripts used to run CMS, can be found at: https://github.com/beacalvo/FMP/tree/master/runCMS.

### 5.2. Computing the effective number of tests (ENT)


In order to account for all variants and all metabolites tested, the p-value threshold used to determine significant associations was calculated by dividing the standard genome-wide significance threshold of 5 x 10-8 by the ENT. The ENT is the virtual number of independent tests across the real number of tests performed and it is computed by taking into account the high degree of correlation between metabolites levels. This way, by using the ENT to correct the p-value threshold, overcorrection for multiple testing is prevented (10).

The R code used to estimate the ENT by this method is the following (10): see script 5.2.Computing_ENT.R .

### 5.3. Selecting confounding variables



All models in CMS were adjusted by age, sex, the 20 first principal components (PCs) and the urine sampling type, which refers to the type of urine sample (night only, morning only or pooled sample). In order to select these variables as confounding, correlation of the available variables and metabolite levels was tested so as to identify those variables explaining part of the metabolite variability. The ones explaining more phenotypic variability were the ones included in the model, in order to decrease metabolic variability and, thus, reduce bias such as the one due to batch effect.

Although the cohort variable was statistically significantly associated to metabolite levels (adjusted p-value < 0.05 for 39 of the 44 metabolites), it was decided to adjust the model by the 20 first PCs since there was a high collinearity between cohort and PCs (mean adjusted r2 = 0.6781519, median = 0.8887265) and because it was considered that the PCs better capture undesired structures in the data which can lead to false associations (the 20 first PCs explained 46.2% of the variability).

(See 5.3.Selecting_confounding_variables.R)


### 5.4. Code used to run CMS



As mentioned in Section 2.6., three different analysis were run. Two of them were one positive control and one negative control, in which only genetic variants located in two different loci (FADS and AGXT2) of the genome were tested. The third, and most relevant analysis, tested the association between the 44 urine metabolite levels and the whole dataset of genotype data. These three analysis were run very similarly: the same phenotypes and summary files were used, only the file containing the genotype data differed.


`runCMS` pipeline is composed of 3 main scripts -as stated in Section 2.6.-:

* *Script 1* was run locally as an R script. The arguments used were the folder path for pipeline, the folder for analysis, the phenotypes file (containing confounding variables, as well as urine and serum metabolite levels), the phenotypes summary file, the number of covariates (urine and serum metabolites) selected for each outcome (urine metabolites) and the threshold on outcome variance explained by covariates.

* *Script 2* was run in the cluster, which used SLURM as job scheduler. As previously mentioned, *Script 2* launches a job array in which each jobs runs CMS for a block of SNPs. The arguments used were: the folder path for pipeline, the folder for analysis, the genotypes file (in `.raw` format), the phenotypes summary file, the total number of SNPs analysed and the number of SNPs per block.

* *Script 3* was run locally as a Python script. The arguments used to run it were: the folder path for pipeline, the folder for analysis, the genotypes file (`BIM` file), the phenotypes summary file, the number of blocks of SNPs (i.e. the number of jobs in the array) and the significance threshold.



The following code shows how these scripts were called and which arguments were used, in each of the three analysis that were run.


#### 5.4.1. Positive control (AGXT2 locus)

(See 5.4.1.AGXT2_CMS.R)

#### 5.4.2. Negative control (FADS locus)

(See 5.4.2.FADS_CMS.R)

#### 5.4.3. Whole genome analysis

(See 5.4.3.Whole_genome_CMS.R)

## 6. Comparison of the identified metabQTLs with urinary metabQTLs in adults


Comparison of our results with findings from other studies with adults was done by two approaches: SNP-based and locus-based. In both cases, a reference dataset containing adult urinary metabQTLs identified in other studies was built for comparison. Locus-based comparison was carried out by comparing the LD block number to which the locus corresponded to, instead of directly comparing the locus in which the SNP was found. The corresponding LD block number was addded for each association from the created reference dataset based on the position from the SNP. Our results already had an LD block number for each SNP.

(See 6.Comparison_adults_1.R, 6.Comparison_adults_2.py and 6.Comparison_adults_3.R)


## 7. GxE interaction analysis



For statistically significant locus-metabolite associations, where the metabolite had previously been related to a dietary factor, GxE interactions were tested. Standard linear regression models (`lm()`) were used under the assumption of an additive genetic factor interacting with a dietary factor treated in tertiles. 

(See 7.GxE.R)




## References

1. Aryee MJ, Jaffe AE, Corrada-Bravo H, Ladd-Acosta C, Feinberg AP, Hansen KD, et al. Minfi: A flexible and comprehensive Bioconductor package for the analysis of Infinium DNA methylation microarrays. Bioinformatics. 2014.
2. Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender D, et al. PLINK: A Tool Set for Whole-Genome Association and Population-Based Linkage Analyses. Am J Hum Genet. 2007; 
3. Kaplan J. fastDummies: Fast Creation of Dummy (Binary) Columns and Rows from Categorical Variables [Internet]. 2019. Available from: https://cran.r-project.org/package=fastDummies
4. Yang J, Lee SH, Goddard ME, Visscher PM. GCTA: A tool for genome-wide complex trait analysis. Am J Hum Genet. 2011; 
5. Wickham H. tidyverse: Easily Install and Load “Tidyverse” Packages. R package version 1.0.0. 2016. 
6. Garnier S. viridis: Default Color Maps from “matplotlib.” R package version 0.5.1. 2018. 
7. Ginestet C. ggplot2: Elegant Graphics for Data Analysis. J R Stat Soc Ser A (Statistics Soc. 2011; 
8. Gallois A, Mefford J, Ko A, Vaysse A, Laakso M, Zaitlen N, et al. A comprehensive study of metabolite genetics reveals strong pleiotropy and heterogeneity across time and context. bioRxiv. 2018; 
9. Aschard H, Guillemot V, Vilhjalmsson B, Patel CJ, Skurnik D, Ye CJ, et al. Covariate selection for association screening in multiphenotype genetic studies. Nat Genet. 2017; 
10. Li MX, Yeung JMY, Cherny SS, Sham PC. Evaluating the effective numbers of independent tests and significant p-value thresholds in commercial genotyping arrays and public imputation reference datasets. Hum Genet. 2012; 
