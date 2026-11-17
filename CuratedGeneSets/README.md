# Curated gene sets for tumor biology
## Information of the gene sets
__A. Cancer and cancer cell related gene sets__  
+ [Oncogene and tumor suppresor genes](https://github.com/KaWingLee9/in_house_tools/blob/main/CuratedGeneSets/TumorGene_gs.csv)  
+ [Transcriptional programs of cancer cells across tumor types](https://github.com/KaWingLee9/in_house_tools/blob/main/CuratedGeneSets/PanCancer_gs.csv)  
+ [Cancer antigen genes](https://github.com/KaWingLee9/in_house_tools/blob/main/CuratedGeneSets/CancerAntigen.csv)  
+ [Housekeeping genes](https://github.com/KaWingLee9/in_house_tools/blob/main/CuratedGeneSets/HouseKeeping_gs.csv)  

__B. Cancer type specific gene sets__  
+ [Colorectal cancer](https://github.com/KaWingLee9/in_house_tools/blob/main/CuratedGeneSets/CRC.csv)  
+ [Central nervous system cancer and glioblastoma](https://github.com/KaWingLee9/in_house_tools/blob/main/CuratedGeneSets/CNSCancer.csv)  

__C. TME related gene sets__  
+ [TME evaluation, Immune subtype classification and TLS (Tertiary lymphoid structures)](https://github.com/KaWingLee9/in_house_tools/blob/main/CuratedGeneSets/ImmEva_gs.csv)
+ [T/NK](https://github.com/KaWingLee9/in_house_tools/blob/main/CuratedGeneSets/TNK_gs.csv)  


## Usage of the gene sets
__A. ClusterProfiler__  
``` R
# ORA enrichment
enrich_result=enricher(target_gene_ls,TERM2GENE=gs_table,pAdjustMethod='BH',pvalueCutoff=0.05,qvalueCutoff=0.05)
# GSEA
```
__B. irGSEA__  

``` R

```

__C. AUCell in python__  

``` python
import numpy as np
import pandas as pd
import scanpy as sc
from pyscenic import aucell

adata.X=adata.layer['normalized']
gs_df=pd.read_csv('TCell_gs.csv')
gs=np.unique(gs_df['gsName'])[0:2]
gene_sets={i:list(gs_df[ gs_df['gsName']==i ]['gsGene']) for i in gs}
```

__D. score_genes in scanpy__  

```
for score_name, gene_list in gene_sets.items():
    sc.tl.score_genes(adata,gene_list=gene_list,score_name=score_name)

sc.pl.embedding(adata,color=gs,basis='X_umap',cmap='bwr')

# using aucell (not finalized)
# exp_mat=pd.DataFrame(adata.layer['normalized'].toarray(),index=adata.obs_names,columns=adata.var_names)

# for score_name, gene_list in gene_sets.items():
#     aucell_df=aucell(exp_mat, gs, auc_threshold=gene_thres, **kwargs)
```
__Reference__: Itay Tirosh, Benjamin Izar, Sanjay M Prakadan, Marc H Wadsworth, Daniel Treacy, John J Trombetta, Asaf Rotem, Christopher Rodman, Christine Lian, George Murphy, and others. Dissecting the multicellular ecosystem of metastatic melanoma by single-cell rna-seq. Science, 352(6282):189–196, 2016.
