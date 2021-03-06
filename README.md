# autoLM
Gene-environment association test using spatial autocorreltation

## Summary

autoLM determines correlations between allele frequencies and environmental measures. It will use simple linear regrestion to determine strongest  GEA correlations. If genetic cluster and population coordinates are included, autoLM uses spatial autocorrelation and population structure as random variables to predict GEAs.

### Input

There are two primary input files and one optional file (see example files):
1) Allele frequencies for each population, tab deliminted <br>
    Format: First column is population name, following columns are allele frequencies for each SNP. 
2) Measures of environmental variables for each population <br>
    Format: Rows correspond to populations, in the same order as file 1. Do not include population names. Following columns are raw environmental measures for each population. 
3) Population coordinates, in meters (optional), tab deliminted  <br>
    Format: First column is population name, second column is genetic cluster (integer), third is east (meters), fourth is north (meters). 
    

### Running 
Use the provided Shell wrapper to run or physically open the R script and edit parameters starting on line 25.

FREQ_FILE= allele frequency file  <br>
ENV_FILE= environmental measure file  <br>
COORDS_FILE= coordinates and genetic cluster file  <br>
STANDARDIZE= if TRUE, all environmental variables will be normalized. If values are already normalized, set to false. 

### Output
autoGLM outputs R2, p-value, mean allele frequency, max difference in allele frequency, and a score for each SNP x environmental variable. The score is a value between 0-1 (low - high) thats weights R2 values by the distribution of and max differences of allele frequencies for each SNP. The goal of the score is to identify promising SNPs with correlations to environmental variables that also show raw variation between populations.


