# Phylogenic analysis to study tumor evolution from omics data  
## Phylogenic analysis from multi-regional analysis from WGS/WES data
Sample data are aquired from __Quantitative evidence for early metastatic seeding in colorectal cancer__ (Supplementary Table 8) and __Mapping the spreading routes of lymphatic metastases in human colorectal cancer__ ().  
### Heterogeneity measure between samples


### Phylogenetic reconstruction


### Metastatic dissemination pattern dissection


## Phylogenic analysis from scRNA-Seq/ST
Sample data are aquired from __Spatial multi-omics landscape of colorectal cancer macro- and micrometastases__ (PAT1 in Figure 2). Processed data acquired from GSE294385 (M-ST-13: CRC, M-ST-05: CLiM, M-ST-14: CLuM). Spot annotation could be found in https://github.com/yliuup/CRC_micromets_ST/blob/main/Meta_data/visium_meta_raw.rds.  
### Phylogenetic reconstruction from somatic mutation
#### Mutation calling from fastq


### Phylogenetic reconstruction from CNV
#### CNV
``` R
library(dplyr)
library(Seurat)
library(infercnv)

selected_cluster=c('Liver macrometastasis tumor','Liver micrometastasis tumor',
  'Lung macrometastasis tumor',
  'Non-neoplastic colon',
  'Tumor cells in mucosa and submucosa','Tumor cells at tumor invasive front/border','Tumor cells in muscularis propria','Tumor cells in muscularis propria')
seurat_obj=subset(seurat_obj,subset= ( Layer3 %in% selected_cluster ))

# run infercnv 
raw_count_mat=as.matrix(seurat_obj@assays[["RNA"]]@counts)
gene_ordering_table=read.table('./gene_ordering_file.txt',
                                  header=FALSE,sep=" ",check.names=FALSE) %>%
   .[! duplicated(.[,1],fromLast=TRUE,),] %>%
   data.frame(row.names=1,check.names=FALSE)
cell_type_annotation=seurat_obj@meta.data['Layer3']

infercnv_obj=infercnv::CreateInfercnvObject(raw_counts_matrix=raw_count_mat,
                                            annotations_file=cell_type_annotation, 
                                            gene_order_file=gene_ordering_table, 
                                            ref_group_names=c('Non-neoplastic colon'), 
                                            chr_exclude=c("chrX","chrY","chrM","KI270734.1"))
infercnv_obj=infercnv::run(infercnv_obj,
                           cutoff=0.1, 
                           out_dir='./infercnv_with_reference',
                           num_threads=15,
                           cluster_by_groups=FALSE,
                           denoise=TRUE,
                           HMM=FALSE)
```
Based on the CNV profiles, identify phylogenetic clones and construct of phylogenetic tree.  
``` r

```

__Reference__:  
[1] Erickson, A., He, M., Berglund, E., Marklund, M., Mirzazadeh, R., Schultz, N., Kvastad, L., Andersson, A., Bergenstråhle, L., Bergenstråhle, J., Larsson, L., Alonso Galicia, L., Shamikh, A., Basmaci, E., Díaz De Ståhl, T., Rajakumar, T., Doultsinos, D., Thrane, K., Ji, A. L., Khavari, P. A., … Lundeberg, J. (2022). Spatially resolved clonal copy number alterations in benign and malignant tissue. Nature, 608(7922), 360–367. https://doi.org/10.1038/s41586-022-05023-2
[2] Liu, Y., Jadhav, A. S., Pan, Y., Liao, J., Khanduri, I., Liu, Y., Katkhuda, R., Lu, W., Cho, K. S., Tong, Z., Zhou, T., Lin, K., Sun, B., Jiang, M., Hernandez, S. D., Lubo Julio, I. C., Brennan, P., Pei, G., Yu, K., Dai, Y., … Maru, D. (2026). Spatial multi-omics landscape of colorectal cancer macro- and micrometastases. Cancer cell, 44(8), 1568–1586.e11.  
[3] Pei, G., Min, J., Rajapakshe, K. I., Branchi, V., Liu, Y., Selvanesan, B. C., Thege, F., Sadeghian, D., Zhang, D., Cho, K. S., Chu, Y., Dai, E., Han, G., Li, M., Yee, C., Takahashi, K., Garg, B., Tiriac, H., Bernard, V., Semaan, A., … Maitra, A. (2025). Spatial mapping of transcriptomic plasticity in metastatic pancreatic cancer. Nature, 642(8066), 212–221.  
[4] Yu, K., Chen, J., Chu, Y. Y., Nair, S., Crupi, E., Hasanov, E., Mei, Y., Han, X., Liu, Y., Liu, Y., Pei, G., Peng, F., Liao, J., Dai, E., Chu, T., Cho, K. S., Jiang, J., Yan, X., Dai, Y., Wang, J., … Wang, L. (2026). A Spatial Atlas of Muscle-Invasive Bladder Cancer Reveals Lineage-Specific Vulnerabilities and Immune Architecture. Cancer discovery, 10.1158/2159-8290.CD-26-0099. https://doi.org/10.1158/2159-8290.CD-26-0099  

## Phylogenic analysis from scRNA-Seq/ST with paired WES/WGS  


## Phylogenic analysis from scDNA-Seq
Sample data are aquired from .


