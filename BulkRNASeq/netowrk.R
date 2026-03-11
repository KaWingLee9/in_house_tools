# build regulatory network
BuildNetwork=function(exp_mat,method='WGCNA',
                      power_vector=1:10,power_used=NULL,
                      minModuleSize=100,maxBlockSize=5000,saveTOMFileBase='SpodopteraTOM-blockwise',
                      regulators=NULL){

    if (method=='GENIE3'){
        
        library(GENIE3)
        weighted_mat=GENIE3(exp_mat,regulators=regulators)
        return(weighted_mat)

    }
    
    exp_mat=log2(exp_mat+1)
    exp_mat=t(exp_mat)

    if (method=='pearson'){
        cor_mat=cor(exp_mat,method='pearson')
        return(cor_mat)
    }

    if (method=='spearman'){
        cor_mat=cor(exp_mat,method='spearman')
        return(cor_mat)
    }
    
    if (method=='WGCNA'){
        
        library(WGCNA)
        
        if (is.null(power_used)){
            sft=pickSoftThreshold(exp_mat,powerVector=power_vector)
            power_recommended=sft$powerEstimate
            print(paste0('Recommended power: ',power_recommended))
        } else if (is.numeric(power_used)) {
            power_recommended=power_used
        }
        wgcna_network=blockwiseModules(exp_mat,maxBlockSize=maxBlockSize,minModuleSize=minModuleSize,
                                       power=power_recommended,TOMType='signed',
                                       reassignThreshold=0,mergeCutHeight=0.25,numericLabels=TRUE,
                                       saveTOMs=TRUE,saveTOMFileBase=saveTOMFileBase,verbose=0)
        module_df=data.frame(gene_id=names(wgcna_network$colors),module=wgcna_network$colors)
        
        return(module_df)
    }
}