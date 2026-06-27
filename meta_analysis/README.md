# Meta analysis

## Rank combination
Combine the rank of P-value.  
``` r
library(Seurat)
library(SeuratData)

InstallData('pbmcsca')
seurat_obj=LoadData('pbmcsca')

table(pbmcsca@meta.data[,c('Experiment','Method')])
#           Method
# Experiment 10x Chromium (v2) 10x Chromium (v2) A 10x Chromium (v2) B
#      pbmc1                 0                3222                3222
#      pbmc2              3362                   0                   0
#           Method
# Experiment 10x Chromium (v3) CEL-Seq2 Drop-seq inDrops Seq-Well Smart-seq2
#      pbmc1              3222      253     3222    3222     3222        253
#      pbmc2                 0      273     3362    3362      551        273

lapply


combined_result=CombRank_DFLs(test_result_ls,p_col='p_val',ES_col='avg_log2FC',
                              min_num=0,min_ratio=0.5,method='RankSum',quantile_est=1/4)
```
Parameters of the function:  
+ `l`: list of the test result of each samples/studie  
+ `p_col`: column name of p_val  
+ `ES_col`: column name of ES  
+ `min_num`: minimal study number with test result  
+ `min_ratio`: % of the results (ignoring NA) have consistent effect would be considered  
+ `method`: "RankSum", "RankProd"  
+ `ncores`: number of cores to use  
Return a data frame of meta analysis:  


'Combined rank - Median effect size' plot could be used for visualization, just as volcano plot:  

`signed_combined_rank` could be used for GSEA.  
__Citation__: Pan-cancer analysis of spatial transcriptomics reveals heterogeneous tumor spatial microenvironment  

