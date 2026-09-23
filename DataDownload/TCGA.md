### ==============================
### Ways to download TCGA data
1. GDCRNATools
2. RTCGA (not updated)
3. TCGABiolinks
4. GDC Data Portal: https://portal.gdc.cancer.gov/
5. UCSC Xena: https://xenabrowser.net/datapages/
6. Broad Institute GDAC: https://gdac.broadinstitute.org/
7. cBioPortal: https://www.cbioportal.org/datasets

### ==============================
### Records for data download
Last update since download: 2022-03-29 (Data Release 32.0)
Reference Genome: hg38
Mapping software: STAR
Transcriptome Profiling-Gene Expression Quantification: 2023-11-21 (mRNA)

### Codes for data download
### TCGAbiolinks
# Codes were from: https://mp.weixin.qq.com/s/m8w1L4N2aXAIers_ZJvp_g, https://mp.weixin.qq.com/s/wI0_GyVl5LiKAjX5C3f-NQ, https://zhuanlan.zhihu.com/p/556196846

```r
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinks") # version should be higher than 2.25.1

library(dplyr)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)

GDCquery('TCGA-BRCA')
--------------------------------------
o GDCquery: Searching in GDC database
--------------------------------------
Genome of reference: hg38

| file_count| case_count|data_category                |
|----------:|----------:|:----------------------------|
|      17337|       1098|Simple Nucleotide Variation  |
|       9281|       1098|Sequencing Reads             |
|       5316|       1098|Biospecimen                  |
|       2288|       1098|Clinical                     |
|      12292|       1098|Copy Number Variation        |
|       4876|       1097|Transcriptome Profiling      |
|       3714|       1097|DNA Methylation              |
|        919|        881|Proteome Profiling           |
|        226|        101|Somatic Structural Variation |
|       4924|       1095|Structural Variation         |
Error in checkDataCategoriesInput(project, data.category, legacy) :
  Please set a data.category argument from the column data_category above
```

### ==============================
### 1. WGS/WES

### Records for data download
Last update since download: Download day 2024.7.20                         2022-03-29 (Data Release 32.0)
Reference Genome: hg38
```r
projects=TCGAbiolinks::getGDCprojects()$project_id %>% .[grepl('^TCGA',.,perl=TRUE)]

sapply(projects, function(project){
  query=GDCquery(project = project,
                    data.category="Simple Nucleotide Variation",
                    data.type = "Masked Somatic Mutation",
                    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking") # legacy: TRUE-hg19 (GDC Legacy Archive), FALSE-hg38 (GDC harmonized database, default)
  GDCdownload(query)
  GDCprepare(query,save=T,save.filename=paste0("DNA/RData/",project,"_SNP.Rdata")) # 100 files download in each time
  data=GDCprepare(query,save=FALSE) 
})

|data_type in query          |
|:---------------------------|
|Aggregated Somatic Mutation |
|Annotated Somatic Mutation  |
|Masked Somatic Mutation     |
|Raw Simple Somatic Mutation |
|Simple Germline Variation   |

> unique(query$results[[1]]$analysis_workflow_type)
 [1] "MuTect2 Annotation"
 [2] "VarScan2 Annotation"
 [3] "VarScan2"
 [4] "MuTect2"
 [5] "Aliquot Ensemble Somatic Variant Merging and Masking"
 [6] "MuSE"
 [7] "Pindel"
 [8] "MuSE Annotation"
 [9] "CaVEMan"
[10] "Pindel Annotation"
[11] "Birdseed"
```

# MC3 pipeline
# Reference: Scalable Open Science Approach for Mutation Calling of Tumor Exomes Using Multiple Genomic Pipelines
```r
data=getMC3MAF()
```

### ==============================
### 2. RNA
### Records for data download
Last update since download: 2022-03-29 (Data Release 32.0)
Reference Genome: hg38
Mapping software: STAR
Transcriptome Profiling-Gene Expression Quantification: 2023-11-21 (mRNA)
```r
library(dplyr)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)

# get TCGA tumor names
projects=TCGAbiolinks::getGDCprojects()$project_id %>% .[grepl('^TCGA',.,perl=TRUE)]

sapply(projects, function(project){
  query=GDCquery(project = project,
                    data.category="Transcriptome Profiling",
                    data.type="Gene Expression Quantification",
                    workflow.type="STAR - Counts") # legacy: TRUE-hg19 (GDC Legacy Archive), FALSE-hg38 (GDC harmonized database, default)
  GDCdownload(query, files.per.chunk = 100)
  GDCprepare(query,save=T,save.filename=paste0("mRNA/RData/",project,"_mRNA.Rdata")) # 100 files download in each time
  data=GDCprepare(query,save=FALSE)
})

> unique(rowdata$gene_type)
 [1] "protein_coding"                     "transcribed_unitary_pseudogene"
 [3] "transcribed_unprocessed_pseudogene" "processed_pseudogene"
 [5] "lncRNA"                             "polymorphic_pseudogene"
 [7] "unprocessed_pseudogene"             "transcribed_processed_pseudogene"
 [9] "IG_V_pseudogene"                    "unitary_pseudogene"
[11] "TR_V_pseudogene"                    "IG_V_gene"
[13] "snRNA"                              "miRNA"
[15] "misc_RNA"                           "snoRNA"
[17] "rRNA_pseudogene"                    "rRNA"
[19] "TR_V_gene"                          "Mt_tRNA"
[21] "Mt_rRNA"                            "IG_C_gene"
[23] "IG_J_gene"                          "TR_J_gene"
[25] "TR_C_gene"                          "TR_J_pseudogene"
[27] "IG_D_gene"                          "ribozyme"
[29] "IG_C_pseudogene"                    "TR_D_gene"
[31] "TEC"                                "IG_J_pseudogene"
[33] "scRNA"                              "scaRNA"
[35] "translated_processed_pseudogene"    "vault_RNA"
[37] "sRNA"                               "translated_unprocessed_pseudogene"
[39] "pseudogene"                         "IG_pseudogene"
```

### ==============================
### 2A. mRNA

```r
sapply(projects, function(project){
  
  load(paste0("mRNA/RData/",project,"_mRNA.Rdata"))

  rowdata=rowData(data)

  # mRNA
  data_mrna=data[rowdata$gene_type=="protein_coding",]
  symbol_mrna=rowData(data_mrna)$gene_name

  # raw counts data
  expr_mat=assay(data_mrna,"unstranded")
  expr_mat=cbind(data.frame(symbol_mrna),as.data.frame(expr_mat))
  expr_mat=expr_mat %>% 
    as_tibble() %>% 
    mutate(meanrow = rowMeans(.[,-1]), .before=2) %>% 
    arrange(desc(meanrow)) %>% 
    distinct(symbol_mrna,.keep_all=T) %>% 
    select(-meanrow) %>% 
    column_to_rownames(var = "symbol_mrna") %>% 
    as.data.frame()
  write.table(expr_mat,paste0('mRNA/counts/',project,'_mRNA_counts.txt'),quote=FALSE)

  # TPM
  expr_mat=assay(data_mrna,"tpm_unstrand")
  expr_mat=cbind(data.frame(symbol_mrna),as.data.frame(expr_mat))
  expr_mat=expr_mat %>% 
    as_tibble() %>% 
    mutate(meanrow = rowMeans(.[,-1]), .before=2) %>% 
    arrange(desc(meanrow)) %>% 
    distinct(symbol_mrna,.keep_all=T) %>% 
    select(-meanrow) %>% 
    column_to_rownames(var = "symbol_mrna") %>% 
    as.data.frame()
  write.table(expr_mat,paste0('mRNA/TPM/',project,'_mRNA_TPM.txt'),quote=FALSE)

  # FPKM
  expr_mat=assay(data_mrna,"fpkm_unstrand")
  expr_mat=cbind(data.frame(symbol_mrna),as.data.frame(expr_mat))
  expr_mat=expr_mat %>% 
    as_tibble() %>% 
    mutate(meanrow = rowMeans(.[,-1]), .before=2) %>% 
    arrange(desc(meanrow)) %>% 
    distinct(symbol_mrna,.keep_all=T) %>% 
    select(-meanrow) %>% 
    column_to_rownames(var = "symbol_mrna") %>% 
    as.data.frame()
  write.table(expr_mat,paste0('mRNA/FPKM/',project,'_mRNA_FPKM.txt'),quote=FALSE)

  # FPKM-UQ
  # expr_mat=assay(data_mrna,"fpkm_uq_unstrand")
  # expr_mat=cbind(data.frame(symbol_mrna),as.data.frame(expr_mat))
  # expr_mat=expr_mat %>% 
  #   as_tibble() %>% 
  #   mutate(meanrow = rowMeans(.[,-1]), .before=2) %>% 
  #   arrange(desc(meanrow)) %>% 
  #   distinct(symbol_mrna,.keep_all=T) %>% 
  #   select(-meanrow) %>% 
  #   column_to_rownames(var = "symbol_mrna") %>% 
  #   as.data.frame()
  # write.table(expr_mat,paste0(project,'_mRNA_FPKM_UQ.txt'),quote=FALSE)

  # metadata
  metadata=colData(data) %>% as.data.frame()
  write.table(metadata,paste0('metadata/',project,'_metadata.txt'),quote=FALSE)

})
```

### ==============================
### N. Clinical data

projects=TCGAbiolinks::getGDCprojects()$project_id %>% .[grepl('^TCGA',.,perl=TRUE)]

clinical_file_str='clinical_patient'
sapply(projects, function(project){
  query=GDCquery(project = project,
    data.category = "Clinical",
    data.type = "Clinical Supplement", 
    data.format = "BCR Biotab")

  GDCdownload(query, files.per.chunk = 100)
  data=GDCprepare(query,save=FALSE)

  data_clin=data[[grep(clinical_file_str,names(data),value=TRUE)]] %>% as.data.frame()
  data_clin=data_clin[-c(1:2),]
  rownames(data_clin)=1:nrow(data_clin)

)}