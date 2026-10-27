

## Intra-tumoral and inter-tumoral heterogeneity of cells 
### Intra-tumoral heterogeneity of cells
Conduct NMF within each sample and get consensus metaprograms.  
Python version:  
```python
import pandas as pd
import anndata as ad
import scanpy as sc
from intra_nmf import *

# data normalization
adata.layers['nmf']=scale_by_batch(adata,'sample_id',target_layer='normalized',max_value=10,min_value=0)
# get MPs
nmf_programs=run_sample_nmf(adata,'sample_id',layer='nmf',min_cells_per_sample=10,n_jobs=5)
MPs_BeforeFilter=derive_metaprograms(nmf_programs,require_different_k=False)
MPs_AfterFilter=filter_metaprograms(MPs_BeforeFilter)
# genes for each MP
MPs_genes=pd.DataFrame(MPs_AfterFilter.genes)
plot_metaprograms(MPs_AfterFilter,order_mode='mp',program_scope='robust')

# scoring and visualization
gene_sets={i:list( MPs_genes.loc[:,i] ) for i in MPs_genes.columns}
for score_name, gene_list in gene_sets.items():
    sc.tl.score_genes(adata,gene_list=gene_list,score_name=score_name)
sc.pl.embedding(adata,basis='X_umap',color=MPs_genes.columns,cmap='bwr',vmax=1,vmin=0)

# annotation using ORA enrichment analysis
enrichment_result=enrich_metaprograms(MPs_AfterFilter)
enrichment_result=pd.DataFrame(enrichment_result)

# 
```

R version:


__Reference__: 
[1] Tirosh I, Suva ML. Cancer cell states: Lessons from ten years of single-cell RNA-sequencing of human tumors. Cancer Cell. 2024 Sep 9;42(9):1497-1506. doi: 10.1016/j.ccell.2024.08.005. Epub 2024 Aug 29.  
[2] Gavish A, Tyler M, Greenwald AC, Hoefflin R, Simkin D, Tschernichovsky R, Galili Darnell N, Somech E, Barbolin C, Antman T, Kovarsky D, Barrett T, Gonzalez Castro LN, Halder D, Chanoch-Myers R, Laffy J, Mints M, Wider A, Tal R, Spitzer A, Hara T, Raitses-Gurevich M, Stossel C, Golan T, Tirosh A, Suvà ML, Puram SV, Tirosh I. Hallmarks of transcriptional intratumour heterogeneity across a thousand tumours. Nature. 2023 Jun;618(7965):598-606. doi: 10.1038/s41586-023-06130-4. Epub 2023 May 31.  
[3] Nomura M, Spitzer A, Johnson KC, Garofano L, Nehar-Belaid D, Galili Darnell N, Greenwald AC, Bussema L, Oh YT, Varn FS, D'Angelo F, Gritsch S, Anderson KJ, Migliozzi S, Gonzalez Castro LN, ChowdhFury T, Robine N, Reeves C, Park JB, Lipsa A, Hertel F, Golebiewska A, Niclou SP, Nusrat L, Kellet S, Das S, Moon HE, Paek SH, Bielle F, Laurenge A, Di Stefano AL, Mathon B, Picca A, Sanson M, Tanaka S, Saito N, Ashley DM, Keir ST, Ligon KL, Huse JT, Yung WKA, Lasorella A, Verhaak RGW, Iavarone A, Suvà ML, Tirosh I. The multilayered transcriptional architecture of glioblastoma ecosystems. Nat Genet. 2025 May;57(5):1155-1167. doi: 10.1038/s41588-025-02167-5. Epub 2025 May 9.  

### Inter-tumoral heterogeneity of cells




__Reference__:  
[1] Tirosh I, Suva ML. Cancer cell states: Lessons from ten years of single-cell RNA-sequencing of human tumors. Cancer Cell. 2024 Sep 9;42(9):1497-1506. doi: 10.1016/j.ccell.2024.08.005. Epub 2024 Aug 29.  
[2] Nomura M, Spitzer A, Johnson KC, Garofano L, Nehar-Belaid D, Galili Darnell N, Greenwald AC, Bussema L, Oh YT, Varn FS, D'Angelo F, Gritsch S, Anderson KJ, Migliozzi S, Gonzalez Castro LN, ChowdhFury T, Robine N, Reeves C, Park JB, Lipsa A, Hertel F, Golebiewska A, Niclou SP, Nusrat L, Kellet S, Das S, Moon HE, Paek SH, Bielle F, Laurenge A, Di Stefano AL, Mathon B, Picca A, Sanson M, Tanaka S, Saito N, Ashley DM, Keir ST, Ligon KL, Huse JT, Yung WKA, Lasorella A, Verhaak RGW, Iavarone A, Suvà ML, Tirosh I. The multilayered transcriptional architecture of glioblastoma ecosystems. Nat Genet. 2025 May;57(5):1155-1167. doi: 10.1038/s41588-025-02167-5. Epub 2025 May 9.  


