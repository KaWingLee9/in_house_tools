# Pipeline for bulk RNA-Seq data analysis

## Normalization of RNA expression matrix
Count matrix -> RPKM/TPM matrix  
__Required packages__: `dplyr`, `GenomicFeatures`  
``` r
source('https://github.com/KaWingLee9/in_house_tools/blob/main/BulkRNASeq/RNANormalization.R')

# remove non-coding genes
count_mat=KeepProteinGene(count_mat,species='human')
# count matrix -> tpm/fpkm matrix
tpm_mat=NormalizeCount(count_matrix,species='human',length.type='transcript',method='tpm')
```

Parameters:  
+ `count_mat`: count matrix with sample x gene  
+ `species`: species of the genes, e.g. `human` (default), `mouse`  
+ `method`: `tpm` (default) or `rpkm`  
+ `length.type`: method to define gene length: `exon` or `transcript` (default)  

## QC of samples
__Required packages__: `ggplot2`, `ComplexHeatmap`  
``` r
source('https://github.com/KaWingLee9/in_house_tools/blob/main/BulkRNASeq/DEAnalysis.R')

p=SampleQC_PCA(tpm_mat,group1,group2,PC=c('PC1','PC2'))
p

ht=SampleQC_cor(log2tpm_mat_scale,group1,group2,method='pearson')
ht
```
Parameters for `SampleQC_PCA`:  
+ `exp_mat`: tpm_mat or fpkm_mat with sample x gene  
+ `group1`, `group2` (optional): vectors to show sample group  
+ `PC`: two PC dimensions for visualization  
+ `plot_label`: whether to show sample label  

Parameters for `SampleQC_cor`:  
+ `exp_mat`: tpm_mat or fpkm matrix with sample x gene  
+ `group1`, `group2` (optional): vectors to show sample group  
+ `method`: method to calculate correlation, could be `pearson` or `spearman`  
+ `...`: other arguments passed on to `ComlexHeatmap::Heatmap`

## Differential gene expression analysis
__Required packages__: DESeq2, edgeR, limma, ggplot2, ggrepel, EnhancedVolcano  
```r
source('https://github.com/KaWingLee9/in_house_tools/blob/main/BulkRNASeq/DEAnalysis.R')

test_result=DiffExp(count_mat,mat_type='count',group=c(),group_ctrl='WT',method='DESeq2')
test_result=data.frame(test_result)
DrawVolcano(test_result,x='log2FoldChange',y='padj',
            pCutoff=0.05,FCcutoff=1.5)

test_result=DiffExp(count_mat,mat_type='count',group=c(),group_ctrl='WT',method='edgeR')
test_result=data.frame(test_result)
```
Parameters for `DiffExp`:
+ `exp_mat`: expression matrix (gene x sample)
+ `mat_type`: one of `count`, `tpm`, `fpkm`  
+ `group`: vector for sample groups  
+ `group_ctrl`: name of control group   
+ `method`: one of `DESeq2`, `edgeR`, `limma`,`wilcox`. __Note that: count_mat for DESeq2 and edgeR; and tpm_mat/fpkm_mat for limma and wilcox!!!__  
+ `test`, `fitType`, `sfType`, `lfc_shrinkage_method` (recommended: ashr): passed on to `DESeq2::DESeq`  
+ `useQL`: whether to use Quasi-likelihood (if FALSE, using v2; if TRUE, using v3/v4)  
+ `NormFactorMethod`,: passed on to `edgeR`  
+ `p_adj_method`: method for p value adjust, passed on to `p.adjust(method=...)`  

Parameters for `DrawVolcano`:  
+ `deg_result`: a data frame resulting fro DEG analysis  
+ `x`, `y`: column name of logfc and p-vale  
+ `FCcutoff`, `pCutoff`: the cutoff of logfc and p-value to select DEGs  
+ `col`: color vector diaplaying significantly upregulated, significantly downregulated, and non-significantly changed genes; should be length of three  
+ `selectLab`: genes that label on the plot    
+ `alpha`: the transparency of the points  
+ `...`: other arguments that passed on to `EnhancedVolcano::EnhancedVolcano`  