### Sample QC

SampleQC_cor=function(exp_mat,group1,group2=NA,method='spearman',...){

    library(ComplexHeatmap)

    exp_mat=log(exp_mat+1)
    exp_mat=scale(exp_mat)
    
    col_list1=c('#E64B35FF','#4DBBD5FF','#00A087FF','#3C5488FF','#F39B7FFF',
                '#8491B4FF','#91D1C2FF','#DC0000FF','#7E6148FF','#B09C85FF')
    col_list2=c('#1F77B4FF','#FF7F0EFF','#2CA02CFF','#D62728FF','#9467BDFF',
                '#8C564BFF','#E377C2FF','#7F7F7FFF','#BCBD22FF','#17BECFFF')
    
    cor_mat=cor(exp_mat,method=method)

    if (is.na(group2)){
        ht=Heatmap(cor_mat,show_column_dend=FALSE,
                   name=paste0(method,'\ncorrelation'),
                   left_annotation=rowAnnotation('group1'=group1,
                                                 col=list('group1'= setNames(col_list1[1:length(unique(group1))],unique(group1)) )),
                   ...)        
    } else {
        ht=Heatmap(cor_mat,show_column_dend=FALSE,
                   name=paste0(method,'\ncorrelation'),
                   left_annotation=rowAnnotation('group1'=group1,'group2'=group2,
                                                 col=list('group1'= setNames(col_list1[1:length(unique(group1))],unique(group1)) ,
                                                          'group2'= setNames(col_list2[1:length(unique(group2))],unique(group2)) )),
                   ...)
    }
    return(ht)
    
}

SampleQC_PCA=function(exp_mat,group1,group2=NA,PC=c('PC1','PC2'),plot_label=FALSE){

    library(ggplot2)

    exp_mat=log(exp_mat+1)
    exp_mat=scale(exp_mat)

    pca_result=prcomp(exp_mat)
    pca_result=summary(pca_result)
    pca_coord=pca_result$rotation
    pca_coord=data.frame(pca_coord)
    pca_coord[,'sample']=rownames(pca_coord)
    pca_variance=pca_result$importance['Proportion of Variance',]

    if (is.na(group2)){
        pca_coord[,'group1']=group1
        p=ggplot(pca_coord,aes_string(x=PC[1],y=PC[2],color='group1'))
    } else {
        pca_coord[,'group']=paste0(group1,'_',group2)
        pca_coord[,'group2']=group2
        p=ggplot(pca_coord,aes_string(x=PC[1],y=PC[2],color='group1',shape='group2'))
    }

    p=p+
        geom_point()+
        theme_bw()+
        xlab(paste0(PC[1],' (',pca_variance[PC[1]]*100,'%)'))+
        ylab(paste0(PC[2],' (',pca_variance[PC[2]]*100,'%)'))

    if (plot_label){
        p=p+geom_text(aes_string(label='sample'))
    }

    return(p)
    
}

### Differential expression analysis
# Required packages: DESeq2, edgeR, limma

# exp_mat: expression matrix (gene x sample)
# group: vector for sample groups
# method: one of `DESeq2`, `edgeR`, `limma`,`wilcox`
# mat_type: one of `count`, `tpm`, `fpkm`
# test, fitType, sfType: passed on to `DESeq`
# useQL: whether to use Quasi-likelihood (if FALSE, using v2; if TRUE, using v3/v4)
# NormFactorMethod,: passed on to `edgeR`
# p_adj_method: method for p value adjust, passed on to p.adjust(method=...)

DiffExp=function(exp_mat,group,group_ctrl=NA,method='DESeq2',mat_type='count',
                 test='Wald',fitType='parametric',sfType='ratio',lfc_shrinkage_method=NA,
                 NormFactorMethod='TMM',useQL=TRUE,
                 p_adj_method='fdr',...){
    
    method=match.arg(method,choices=c('DESeq2','edgeR','limma','wilcox'))
    mat_type=match.arg(mat_type,choices=c('count','tpm','fpkm'))

    if (! is.na(group_ctrl)){
        group=relevel(factor(group),ref=group_ctrl)
    }
    
    if (mat_type=='count'){
        count_mat=exp_mat

        if (method=='DESeq2'){
            
            library(DESeq2)
            
            cds=DESeqDataSetFromMatrix(count_mat,colData=DataFrame(group),design=~group)
            cds=DESeq(cds,test=test,fitType=fitType,sfType=sfType)

            if ( ! is.na(lfc_shrinkage_method) ){
                cds=lfcShrink(cds,type=lfc_shrinkage_method)
            }
            
            test_resut=results(cds,pAdjustMethod=p_adj_method)
    
            DESeq2::plotMA(test_resut)
            # plot(test_resut[,'baseMean'],test_resut[,'log2FoldChange'],
            #      col=ifelse(test_resut[,'padj']<=0.05,'red','black'),xlab='logCounts', ylab='logFC',pch=20)
    
            return(test_resut)
            }
        
            if (method=='edgeR'){
        
                library(edgeR)
                
                y=DGEList(count_mat,group=group)
                y=calcNormFactors(y,method=NormFactorMethod)
                design=model.matrix(~group)
                y=estimateDisp(y,design)
                
                if (useQL){
                    fit=glmQLFit(y,design)
                    test_resut=glmQLFTest(fit)
                } else{
                    fit=glmFit(y,design)
                    test_resut=glmLRT(fit)
                }
                test_resut$table['Padj']=p.adjust(test_resut$table[,'PValue'],method='fdr')
                
                test_table=test_resut$table
                plot(test_table[,'logCPM'],test_table[,'logFC'],
                     col=ifelse(test_table[,'Padj']<=0.05,'red','black'),xlab='logCounts', ylab='logFC',pch=20)
                
                return(test_resut)
                
            }
    }

    if (mat_type=='tpm' | mat_type=='fpkm'){

        # if (method=='limma'){}
        # if (method=='wilcox'){}
        
    }

    
}

# DiffExp_pair

# DiffExp_time=function(exp_mat,,method='')

# Volcano plot for single conditions
# Required packages: ggplot2, ggrepel, EnhancedVolcano
DrawVolcano <- function(deg_result,x='log_fc',FCcutoff=1,
                        y='p_val',pCutoff=0.05,
                        col=c('red2','royalblue','grey30'), #
                        selectLab=NULL,alpha=0.5){
    library(ggplot2)
    library(ggrepel)
    library(EnhancedVolcano)
    
    keyvals=rep(col[3],length.out=nrow(deg_result))
    names(keyvals)=rep('Unchanged',nrow(deg_result))

    keyvals[which((deg_result[,x]>=FCcutoff) & (deg_result[,y]<=pCutoff))]=col[1]
    names(keyvals)[which((deg_result[,x]>=FCcutoff) & (deg_result[,y]<=pCutoff))]='Higher expressed'

    keyvals[which((deg_result[,x]<=FCcutoff) & (deg_result[,y]<=pCutoff))]=col[2]
    names(keyvals)[which((deg_result[,x]<=-FCcutoff) & (deg_result[,y]<=pCutoff))]='Lower expressed'

    p=EnhancedVolcano(deg_result,
                lab=rownames(deg_result),
                x=x,FCcutoff=FCcutoff,
                y=y,pCutoff=pCutoff,
                colCustom=keyvals,
                      selectLab=selectLab,
               drawConnectors = TRUE,arrowheads=FALSE,
                title=NULL,subtitle=NULL,colAlpha=alpha
               )

    p=p+scale_color_manual(limits=c('Lower expressed','Unchanged','Higher expressed'),values=c(col[2],col[3],col[1]))+
        xlab(expression(log[2]~FC))+
        ylab(expression(-log[10]~P))

    return(p)
}


### GSEA
