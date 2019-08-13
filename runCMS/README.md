This pipeline is made to run CMS in large-scale analyses.

# Description

This package is composed of different elements to run a complete analysis with CMS.

- The script 1 aims at preparing the data.
- The script 2 aims at launching the analysis.
- The script 3 aims at parsing the results.
- The folder "test" contains examples of input files (simulated data).
- The folder "others" contains secondary scripts and files.


# Requirements

All the scripts are made to be run on a cluster.  
To use these scripts, you need R version 3 and Python version 3 with the package `Rpy2`.  

Each script can be run with a command line in a terminal.  
First, you need to download runCMS and put it in a dedicated folder.  
Then, you need to create an empty folder for all your analysis, where runCMS will generate intermediate files and results.

# Details

For more details on input formats and usage, see the [wiki](https://gitlab.pasteur.fr/statistical-genetics/runCMS/wikis/home).


# References

For the CMS algorithm:  
Aschard, H. et al. Covariate selection for association screening in multiphenotype genetic studies. *Nat Genet* **49**, 1789-1795 (2017).

For the LD blocks used to summarize results per loci:  
Berisa, T. & Pickrell, J.K. Approximately independent linkage disequilibrium blocks in human populations. *Bioinformatics* **32**, 283-5 (2016).


# License

This project is licensed under GNU General Public License version 2.


# Author

Apolline Gallois (C3BI, Institut Pasteur)