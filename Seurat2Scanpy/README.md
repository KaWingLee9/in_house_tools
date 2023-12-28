Functions to transform Seurat object to Scanpy object

### Tutorial
1. To use the script here, you should prepare:  
A) Python console with: `scanpy`, `numpy`, `pandas`, `rpy2`  
B) R console with: `Seurat`, `dplyr`  
3. Sample `.rds` files of ST datasets could be downloaded from  
4. Download the source file Seurat2Scanpy.py and import the functions:  
```
exec(open('Seurat2Scanpy.py','r').read())  
r_home="~/miniconda3/envs/R/lib/R" # path of R console 
```
4. For scRNA-Seq data:  
```
adata=sc_Seurat2Scanpy('seurat_obj.rds',r_home,
                    exp_mat_slot=['RNA','data'])
```
5. For spatial-omics data:  
```
adata=st_Seurat2Scanpy('seurat_obj.rds',r_home,
                    exp_mat_slot=['RNA','data'],res_type='lowres')
```
`exp_mat_slot` is a list with `assay` and `slot` name of expression matrix in seurat object to export to `adata.X`.

### Citation
Wang H, Li J, Jing S, et al. SOAPy: a Python package to dissect spatial architecture, dynamics and communication[J]. bioRxiv, 2023: 2023.12. 21.572725.
