Functions to transform Seurat object to Scanpy object  
__{For both single-cell and spatial-omics data!!!}__
Results are adapted for SOAPy  
Temporally for Seurat v4   

### Tutorial
1. To use the script here, you should prepare:  
A) Python console with: `scanpy`, `numpy`, `pandas`, `rpy2`. For Windows user to use rpy2, please refer to https://support.anaconda.com/hc/en-us/articles/360023857134-Setting-up-rpy2-on-Windows.  
B) R console with: `Seurat`, `dplyr`.  
2. Download the source file Seurat2Scanpy.py and import the functions:  
```
exec(open('Seurat2Scanpy.py','r').read())  
r_home="~/miniconda3/envs/R/lib/R" # path of R console 
```
3. For scRNA-Seq data:  
```
adata=sc_Seurat2Scanpy('seurat_obj.rds',r_home,
                    exp_mat_slot=['RNA','data'])
```
4. For spatial-omics data:  
```
adata=st_Seurat2Scanpy('seurat_obj.rds',r_home,
                    exp_mat_slot=['RNA','data'],res_type='lowres')
```
`exp_mat_slot` is a list with `assay` and `slot` name of expression matrix in seurat object to export to `adata.X`.

### Citation
Wang H, Li J, Jing S, et al. SOAPy: a Python package to dissect spatial architecture, dynamics and communication[J]. bioRxiv, 2023: 2023.12. 21.572725.
