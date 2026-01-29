# Pipeline for data analysis of bulk RNA-Seq

## Normalization of RNA expression matrix
Count matrix -> RPKM/TPM matrix  
__Required packages__: `dplyr`, `GenomicFeatures`  
``` r
source('https://github.com/KaWingLee9/in_house_tools/blob/main/RNANormalization/RNANormalization.R')

norm_matrix=NormalizeCount(count_matrix,length.type='transcript',method='tpm')
```
Parameters:  
+ `count_mat`: count matrix with sample x gene  
+ `species`: species of the genes, e.g. `human` (default), `mouse`  
+ `method`: `tpm` (default) or `rpkm`  
+ `length.type`: method to define gene length: `exon` or `transcript` (default)  

## Differential gene expression analysis
