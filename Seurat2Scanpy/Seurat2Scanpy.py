# function to transform 10X Visium data
def st_Seurat2Scanpy(seurat_file,r_home,
                     exp_mat_slot=['RNA','counts'],res_type='lowres'):
    
    
    # Python library preparation
    import anndata as ad
    import scanpy as sc
    import numpy as np

    # R console preparation
    import os
    os.environ['R_HOME']=r_home
    from rpy2.robjects import r
    import rpy2.robjects as robjects
    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri
    
    # assert r_home!=None 'R Console is not assigned'

    # Automatically transfer R-data.frame to Python-DataFrame
    pandas2ri.activate()

    # R fucntion preparation
    readRDS=robjects.r['readRDS']
    to_array=robjects.r['as.array']
    to_data_frame=robjects.r['as.data.frame']
    transpose=robjects.r['t']
    
    # required R package loading
    
    Seurat=importr('Seurat')
    SeuratObject=importr('SeuratObject')

    # get simple data as single cell 
    GetMetaData=r("""
    function(seurat_obj){
    library(Seurat)
    metadata=seurat_obj@meta.data
    metadata$Idents=Idents(seurat_obj)
    return(metadata)
    }
    """)

    GetMetaFeature=r("""
    function(seurat_obj,assay='RNA'){
    library(Seurat)
    library(dplyr)
    gene_feature=GetAssay(seurat_obj,assay=assay)@meta.features
    gene_feature$'genome'=1
    return(gene_feature)
    }
    """)

    GetAssayData=r("""
    function(seurat_obj,assay,slot){
    library(Seurat)
    exp_mat=GetAssayData(seurat_obj,assay=assay,slot=slot)
    return(exp_mat %>% as.data.frame())
    }
    """)

    # get spatial data specific specific for ST
    GetSpatial=r("""
    function(seurat_obj)
    return(seurat_obj@images)
    """)

    GetCoordinates_Visium=r("""
    function(image_obj){
    library(Seurat)
    coord_mat=image_obj@coordinates
    return(coord_mat[,c('imagerow','imagecol')])
    }
    """)
    
    GetCoordinates_SlideSeq=r("""
    function(image_obj){
    library(Seurat)
    coord_mat=image_obj@coordinates
    return(coord_mat[,c('x','y')])
    }
    """)
    
    GetCellCoordinates=r("""
    function(image_obj){
    library(Seurat)
    return(GetTissueCoordinates(image_obj[["centroids"]])[,c('x','y')])
    }
    """)
    
    GetSegmentationCoordinates=r("""
    function(image_obj){
    library(Seurat)
    return(GetTissueCoordinates(image_obj[["segmentation"]]))
    }
    """)
    
    GetMoleculeCoordinates=r("""
    function(image_obj,assay,slot){
    library(Seurat)
    library(dplyr)
    x=rownames(GetAssayData(seurat_obj,assay=assay,slot=slot))
    return(GetTissueCoordinates(image_obj[["molecules"]]) %>% filter(molecule %in% x))
    }
    """)

    GetArrayCoordinates=r("""
    function(image_obj){
    library(Seurat)
    coord_mat=image_obj@coordinates[,c('col','row')]
    colnames(coord_mat)=c('array_row','array_col')
    return(coord_mat)
    }
    """)

    GetImage=r("""
    function(image_obj){
    library(Seurat)
    return(image_obj@image)
    }
    """)

    GetScaleFactor=r("""
    function(image_obj){
    library(Seurat)
    return(image_obj@scale.factors)
    }
    """)

    GetImageType=r("""
    function(image_obj){
    library(Seurat)
    return(class(image_obj))
    }
    """)

    import pandas as pd
    import numpy as np
    
    seurat_obj=readRDS(seurat_file)
    
    # raw count matrix
    assay=exp_mat_slot[0]
    slot=exp_mat_slot[1]
    count_mat=SeuratObject.GetAssayData(seurat_obj,assay=assay,slot=slot)

    count_mat=to_array(count_mat)
    count_mat=count_mat.transpose()

    # get gene names
    metafeature=GetMetaFeature(seurat_obj,assay)

    # metadata
    metadata=GetMetaData(seurat_obj)
    
    # create anndata
    adata=ad.AnnData(X=count_mat,var=metafeature,obs=metadata)  
    coord_mat_ls=[]
    
    images=GetSpatial(seurat_obj)
    ImageType=GetImageType(images[0])[0]
    
    # 10X Visium
    if ImageType=='VisiumV1':
        array_coord_mat_ls=[]
        spatial_dir={}
        for i in range(len(images)):
            image_obj=images[i]
            # sample name
            sample_name=images.names[i]
            # coordinate matrix
            coord_mat=GetCoordinates_Visium(image_obj)
            coord_mat=np.array(coord_mat)
            coord_mat_ls.append(coord_mat)
            array_coord_mat=GetArrayCoordinates(image_obj)
            array_coord_mat_ls.append(array_coord_mat)
            # image array
            image_array=SeuratObject.GetImage(image_obj,mode='raw')
            # scale factor
            scalefactors=GetScaleFactor(image_obj)
            scalefactors={scalefactors.names[i]:scalefactors[i][0] for i in range(len(scalefactors))}
            scalefactors={'fiducial_diameter_fullres':scalefactors['fiducial'],
            'tissue_hires_scalef':scalefactors['hires'],
            'tissue_lowres_scalef':scalefactors['lowres']}
            scalefactors['spot_diameter_fullres']=0.6196*scalefactors['fiducial_diameter_fullres']
            spatial_dir[sample_name]={'images':{res_type:image_array},
                                      'scalefactors':scalefactors}
        # add spatial data
        adata.uns['spatial']=spatial_dir
        if len(array_coord_mat_ls)==1:
            array_coord_mat=array_coord_mat_ls[0]
        else:
            array_coord_mat=pd.concat(array_coord_mat_ls,axis=0)
        adata.obs=pd.concat((adata.obs,array_coord_mat),axis=1)
    
    # Slide-Seq
    if ImageType=='SlideSeq':
        for i in range(len(images)):
            image_obj=images[i]
            # sample name
            sample_name=images.names[i]
            # coordinate matrix
            coord_mat=GetCoordinates_SlideSeq(image_obj)
            coord_mat=np.array(coord_mat)
            coord_mat_ls.append(coord_mat)
    
    # STARmap
    # if ImageType=='STARmap':
        
        
    # MERSCOPe, Xenium, Nanostring CosMx, CODEX, etc.
    if ImageType=='FOV':
        for i in range(len(images)):
            image_obj=images[i]
            # sample name
            sample_name=images.names[i]
            # coordinate matrix
            coord_mat=GetCellCoordinates(image_obj)
            coord_mat=np.array(coord_mat)
            coord_mat_ls.append(coord_mat)
            # image_data=GetMoleculeCoordinates(image_obj,assay,slot)
            # segmentation_data=GetSegmentationCoordinates=(image_obj)      
    
    # add coordination matrix
    if len(coord_mat_ls)==1:
        coord_mat=coord_mat_ls[0]
    else:
        coord_mat=np.concatenate(coord_mat_ls,axis=0)
    coord_mat=coord_mat[:,[1,0]]
    adata.obsm['spatial']=coord_mat

    return(adata)

# function to transform scRNA-Seq data
def sc_Seurat2Scanpy(seurat_file,r_home,
                     exp_mat_slot=['RNA','counts']):
    
    
    # Python library preparation
    import anndata as ad
    import scanpy as sc
    import numpy as np

    # R console preparation
    import os
    os.environ['R_HOME']=r_home
    from rpy2.robjects import r
    import rpy2.robjects as robjects
    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri
    
    # assert r_home!=None 'R Console is not assigned'

    # Automatically transfer R-data.frame to Python-DataFrame
    pandas2ri.activate()

    # R fucntion preparation
    readRDS=robjects.r['readRDS']
    to_array=robjects.r['as.array']
    to_data_frame=robjects.r['as.data.frame']
    transpose=robjects.r['t']
    
    # required R package loading
    
    Seurat=importr('Seurat')
    SeuratObject=importr('SeuratObject')

    # get simple data as single cell 
    GetMetaData=r("""
    function(seurat_obj){
    library(Seurat)
    metadata=seurat_obj@meta.data
    return(metadata)
    }
    """)

    GetMetaFeature=r("""
    function(seurat_obj,assay='RNA'){
    library(Seurat)
    library(dplyr)
    gene_feature=GetAssay(seurat_obj,assay=assay)@meta.features
    gene_feature$'genome'=1
    return(gene_feature)
    }
    """)

    GetAssayData=r("""
    function(seurat_obj,assay,slot){
    library(Seurat)
    exp_mat=GetAssayData(seurat_obj,assay=assay,slot=slot)
    return(exp_mat %>% as.data.frame())
    }
    """)


    import pandas as pd
    import numpy as np
    
    seurat_obj=readRDS(seurat_file)
    
    # raw count matrix
    assay=exp_mat_slot[0]
    slot=exp_mat_slot[1]
    count_mat=SeuratObject.GetAssayData(seurat_obj,assay=assay,slot=slot)

    count_mat=to_array(count_mat)
    count_mat=count_mat.transpose()

    # get gene names
    metafeature=GetMetaFeature(seurat_obj,assay)

    # metadata
    metadata=GetMetaData(seurat_obj)
    

    # create anndata
    adata=ad.AnnData(X=count_mat,var=metafeature,obs=metadata)

    return(adata)
