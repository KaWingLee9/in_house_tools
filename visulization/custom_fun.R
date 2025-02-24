# CorPlot - Correlation bubble plot with significance test
# Required packages: ggplot2, Hmisc, dendsort
CorPlot=function(df,cor.method='pearson', # 'pearson', 'spearman'
        size='p.value', # 'p.value', 'p.adj'
        p.adj.method='fdr',
        sig.level=0.05,
        tri='whole', # 'whole', 'lower', 'upper'
        sig.circle=TRUE,reorder.method='reorder.dendrogram'){
    
    library(ggplot2)

    # cor_result and p_matirx could be replaced
    cor_result=Hmisc::rcorr(as.matrix(df),type=cor.method)
    cor_matrix=cor_result$r
    p_matrix=cor_result$P

    # clustering 
    hc=hclust(as.dist(1-cor_matrix),method='ward.D2')

    # reorder hierarchy
    if (reorder.method=='reorder.dendrogram'){
        hc=reorder(as.dendrogram(hc),-rowMeans(cor_matrix))
        hc=as.hclust(hc)
    } else if (reorder.method=='dendsort'){
        hc=dendsort::dendsort(hc)
    } else if (reorder.method=='none'){
        hc=hc
    }
    
    ord.order=hc[['labels']][hc[['order']]]
    cor_matrix=cor_matrix[ord.order,ord.order]

    if (tri=='lower'){
        cor_matrix[upper.tri(cor_result$r,diag=FALSE)]=NA
    }
    if (tri=='upper'){
        cor_matrix[lower.tri(cor_result$r,diag=FALSE)]=NA
    }

    test_result=reshape2::melt(cor_matrix,varnames=c('x','y'),value.name='Correlation')
    test_result$p.value=apply(test_result,1,function(x) p_matrix[x[1],x[2]])
    test_result$p.value[is.na(test_result$Correlation)]=NA
    test_result$p.adj=p.adjust(test_result$p.value,method=p.adj.method)
    test_result$sig=test_result[,size]<=sig.level
    # test_result$sig[is.na(test_result$sig)]=FALSE
    # test_result[test_result[,'x']==test_result[,'y'],'Correlation']=NA
    test_result[test_result[,'x']==test_result[,'y'],'p.value']=0.1
    test_result[test_result[,'x']==test_result[,'y'],'p.adj']=0.1

    if (!sig.circle){
        p=ggplot(test_result,aes_string(x='x',y='y',fill='Correlation',size=size))+
          geom_point(shape=21,color='#FFFFFF00')
    }

    if(sig.circle){
        p=ggplot(test_result,aes_string(x='x',y='y',fill='Correlation',size=size,color='sig'))+
          scale_color_manual(values=c('TRUE'='black','FALSE'='#FFFFFF00'),na.value='#FFFFFF00',limits=c(TRUE,FALSE))+
          geom_point(shape=21,stroke=1)
    }

    p=p+
      scale_fill_gradient2(low='blue',mid='white',high='red')+
      scale_size(trans='reverse')+
      theme_minimal()+ 
      coord_fixed()+
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank())+
      theme(panel.background = element_blank(), 
            panel.grid.major = element_blank(),
            panel.border = element_blank())+
      labs(fill=stringr::str_to_title(paste0(cor.method,'\ncorrelation')))+
      theme(plot.title=element_text(hjust=0.5))+
      labs(color=paste0('Sig\ (',size,'<=',sig.level,')'))+
      scale_size_continuous(range=c(6,0.2))+
      guides(fill=guide_colorbar(order=1),size=guide_legend(order=2),color=guide_legend(order=3))

    if (tri=='lower'){
        p=p+
          scale_y_discrete(position='right')+
          theme(axis.text.x=element_text(angle=45,hjust=1))
    }

    if (tri=='upper'){
        p=p+
          scale_x_discrete(position='top')+
          theme(axis.text.x=element_text(angle=45,hjust=0))
    }

    if (tri=='whole'){
        p=p+
          theme(axis.text.x=element_text(angle=45,hjust=1))
    }

    if (tri=='upper'){
      n=colnames(cor_matrix)
      p=p+geom_text(data=data.frame(x=1.5:(length(n)+0.5),y=1:(length(n)+0),label=n),
                    aes(x=x,y=y,label=label),
                    angle=0, # 0/-45
                    inherit.aes=FALSE,hjust=0)+
      theme(axis.text.x=element_blank(),axis.text.y=element_blank())
    }

    if (tri=='lower'){
      n=colnames(cor_matrix)
      p=p+geom_text(data=data.frame(x=0.5:(length(n)-0.5),y=1:(length(n)+0),label=n),
                    aes(x=x,y=y,label=label),
                    angle=0, # 0/-45
                    inherit.aes=FALSE,hjust=1)+
      theme(axis.text.x=element_blank(),axis.text.y=element_blank())
    }
                              
    return(p)
}

# SHeatmap - Summarized heatmap with significance test
# Required packages: dplyr, coin, ComplexHeatmap
SumHeatmap=function(df,group.col,variable.col,value.col,test.mode='ONEvsVALUE',
                    test.method='t.test',permutated=FALSE,
                    show.significance=TRUE,sig.level=c(0.01,0.05),sig.label=c('**','*'),
                    p.adj=FALSE,p.adj.method='fdr',scale=TRUE,transpose=FALSE,...){
    
    options(warn=-1)
    
    library(dplyr)
    
#     if (scale){
#         df_1=df %>% group_by(!!!syms(variable.col)) %>% mutate(value.col=scale(!!!syms(value.col)))
#         df_1[,value.col]=df_1[,'value.col']
#         df=data.frame(df_1)
#         # heatmap_matrix=scale(heatmap_matrix)
#     }
    
    if (test.method=='wilcox.test'){
        heatmap.aggr.fun=median
    } else {
        heatmap.aggr.fun=mean
    }
    
    heatmap_matrix=reshape2::dcast(df,as.formula(paste0(group.col,'~',variable.col)),value.var=value.col,fun.aggregate=heatmap.aggr.fun) %>% 
        data.frame(row.names=1,check.names=FALSE)
    
    if (scale){
        heatmap_matrix=scale(heatmap_matrix)
    }
    
    if (show.significance){
        p_matrix=lapply(unique(df[,variable.col]),function(x){

            m=unique(df[,group.col])
            df=sapply(m,function(y){

                if (test.mode=='ONEvsVALUE'){
                    x1=df[ (df[,group.col] %in% y) & (df[,variable.col] %in% x) ,value.col]
                    x2=mean(df[ df[,variable.col] %in% x ,value.col])

                    if (test.method=='t.test'){
                        x2=mean(df[ df[,variable.col] %in% x ,value.col])
                        p=t.test(x1,mu=x2)$p.value
                    }

                    if (test.method=='wilcox.test'){
                        x2=median(df[ df[,variable.col] %in% x ,value.col])
                        p=wilcox.test(x1,mu=x2)$p.value
                    }

                    return(p)
                }

                if (test.mode=='ONEvsOTHER'){
                    x1=df[ (df[,group.col] %in% y) & (df[,variable.col] %in% x) ,value.col]
                    x2=df[ (!df[,group.col] %in% y) & ((df[,variable.col] %in% x)) ,value.col ]
                }

                if (test.mode=='ONEvsALL'){
                    x1=df[ (df[,group.col] %in% y) & (df[,variable.col] %in% x) ,value.col]
                    x2=df[ df[,variable.col] %in% x ,value.col]
                }

                if (permutated){

                    df_test=data.frame(value=c(x1,x2),
                    group=as.factor(c(rep('A',length.out=length(x1)),rep('B',length.out=length(x2)))))
                    colnames(df_test)=c('value','group')
                    df_test=data.frame(df_test)

                    if (test.method=='oneway.test'){
                        p=coin::pvalue(coin::oneway_test(value~group,df_test))
                    }

                    if (test.method=='wilcox.test'){
                        p=coin::pvalue(coin::wilcox_test(value~group,df_test))
                    }
                }

                if (!permutated){
                    if (test.method=='t.test'){
                        p=t.test(x1,x2)$p.value
                    }

                    if (test.method=='wilcox.test'){
                        p=wilcox.test(x1,x2)$p.value
                    }
                }

                return(p)

            },USE.NAMES=TRUE) %>% data.frame()

            rownames(df)=m
            colnames(df)=x
            return(df)

        }) %>% dplyr::bind_cols()

        p_matrix=p_matrix[rownames(heatmap_matrix),colnames(heatmap_matrix)]

        if (p.adj){
            p_matrix_adj=unlist(p_matrix) %>% p.adjust(method=p.adj.method) %>% matrix(c(nrow(p_matrix),ncol(p_matrix)))
            rownames(p_matrix_adj)=rownames(p_matrix)
            colnames(p_matrix_adj)=colnames(p_matrix)
            p_matrix=p_matrix_adj
        }
    }
    
    if (transpose){
        heatmap_matrix=t(heatmap_matrix)
        
        if (show.significance){
            p_matrix=t(p_matrix)
        }
    }

    library(ComplexHeatmap)
    ht=Heatmap(heatmap_matrix,cell_fun=function(j,i,x,y,w,h,fill){
        if (show.significance){
            q=min(which(p_matrix[i,j]<=sig.level))
            if (q<=length(sig.level)){
                grid.text(sig.label[q],x,y)
            } else{
                grid.text('',x,y)
            }
        } else {
            grid.text('',x,y)
        }

    },...)
    
    options(warn=1)
    return(ht)

}

                              
# SimilarityClustering - Blocks division in similarity heatmap
# Required packages: NbClust, simplifyEnrichment, ComplexHeatmap, ConsensusClusterPlus
SimilarityClustering=function(data,mode='automatic',select.cutoff=FALSE,
                              similarity.method='pearson',use.fastcluster=FALSE,hc.method='ward.D2',
                              provided_label=NA,min.nc=2,max.nc=15,show_index_result=NA,cluster_num=0,
                              cutoff.seq=seq(0.6,0.98,by=0.01),cutoff=0.85,
                              maxK=10,...){
        
  library(ComplexHeatmap)

  if (!((mode=='ConsensusClusterPlus') & (select.cutoff==FALSE))) {
      
      if (similarity.method=='euclidean') {
          similarity_matrix=dist(data,method='euclidean') %>% as.matrix()
          similarity_matrix=-similarity_matrix
      }
      if (similarity.method %in% c('pearson','spearman')){
          similarity_matrix=cor(t(data),method=similarity.method)
      }
      similarity_matrix[is.na(similarity_matrix)]=0
  }

  col_type=c('#5050FFFF','#CE3D32FF','#749B58FF','#F0E685FF','#466983FF','#BA6338FF','#5DB1DDFF','#802268FF','#6BD76BFF','#D595A7FF','#924822FF',
    '#837B8DFF','#C75127FF','#D58F5CFF','#7A65A5FF','#E4AF69FF','#3B1B53FF','#CDDEB7FF','#612A79FF','#AE1F63FF','#E7C76FFF','#5A655EFF',
    '#CC9900FF','#99CC00FF','#A9A9A9FF','#CC9900FF','#99CC00FF','#33CC00FF','#00CC33FF','#00CC99FF','#0099CCFF','#0A47FFFF','#4775FFFF',
    '#FFC20AFF','#FFD147FF','#990033FF','#991A00FF','#996600FF','#809900FF','#339900FF','#00991AFF','#009966FF','#008099FF','#003399FF',
    '#1A0099FF','#660099FF','#990080FF','#D60047FF','#FF1463FF','#00D68FFF','#14FFB1FF')

  if (use.fastcluster){
    library(fastcluster)
  }

  if (mode=='automatic'){
          
    library(simplifyEnrichment)
          
    if (select.cutoff){
        
        if (hc.method=='kmeans'){
            partition_by_kmeans=function(mat,n_repeats=10) {
                partition_list=lapply(seq_len(n_repeats),function(i) {
                    as.cl_hard_partition(kmeans(mat,2))
                })
            }
        } else {
            partition_fun=function(mat) {
                cutree(hclust(dist(mat),method=hc.method),2)
            }
        }
      return(select_cutoff(similarity_matrix,cutoff=cutoff.seq,verbose=FALSE,partition_fun=partition_by_hclust))
        
    } else {
  
      r=rownames(similarity_matrix)
      c=binary_cut(similarity_matrix,cutoff=cutoff,partition_fun=partition_by_hclust)
      or=r[order(c)]
      c=as.factor(c)
  
      col_type=col_type[1:max(as.numeric(c))]
      names(col_type)=1:max(as.numeric(c))
      
      print(  Heatmap(similarity_matrix,cluster_rows=FALSE,cluster_columns=FALSE,row_order=or,column_order=or,
              show_row_names=FALSE,show_column_names=FALSE,name='Similarity\nmatrix',
              row_split=c,column_split=c,
              left_annotation=rowAnnotation(' '=c,col=list(' '=col_type),show_legend=FALSE),
              top_annotation=HeatmapAnnotation(' '=c,col=list(' '=col_type),show_legend=FALSE),...)  )
      
      names(c)=r
      return(c)
    }
  } 
        
  if (mode=='manual'){
      
    diss_mat=1-similarity_matrix
      
    if (select.cutoff){
        
        if (is.na(show_index_result)){
            
            test_index_1=c("kl","ch","hartigan","ccc","scott","marriot","trcovw","tracew","friedman","rubin",
                            "cindex", "db", "silhouette", 
                            "duda", "pseudot2", "beale", 
                            "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", 
                            "gamma", "gplus", "tau", "dunn","sdindex","sdbw")
            

            test_index_2=c("hubert","dindex")
            
            print(sapply(test_index_1,function(x){
                res=try({
                    y=NbClust::NbClust(data=data,diss=as.dist(diss_mat),distance=NULL,min.nc=min.nc,max.nc=max.nc, 
                                       method=hc.method,alphaBeale=0.1,index=x)                
                },silent=TRUE)
                if (inherits(res,'try-error')) {return(NULL)}
                return(y$Best.nc[['Number_clusters']])
            }))

            sapply(test_index_2,function(x){
                res=try({
                      y=NbClust::NbClust(data=data,diss=as.dist(diss_mat),distance=NULL,min.nc=min.nc,max.nc=max.nc, 
                         method=hc.method,alphaBeale=0.1,index=x)
                },silent=TRUE)
                if (inherits(res,'try-error')) {return(NULL)}
                return(NULL)
            })

          return(NULL)
        } else {
            
            y=NbClust::NbClust(data=data,diss=as.dist(diss_mat),distance=NULL,min.nc=min.nc,max.nc=max.nc, 
                   method=hc.method,alphaBeale=0.1,index=show_index_result)
            plot(min.nc:max.nc,y[['All.index']],xlab='Cluster number',ylab=show_index_result)
            
        }
        
        

    } else {
      
      if (!is.na(provided_label)){
          c=as.factor(provided_label)
          col_type=col_type[1:length(unique(c))]
          names(col_type)=unique(c)
          print(  Heatmap(similarity_matrix,cluster_rows=FALSE,cluster_columns=FALSE,
              show_row_names=FALSE,show_column_names=FALSE,name='Similarity\nmatrix',
              row_split=c,column_split=c,
              left_annotation=rowAnnotation(' '=c,col=list(' '=col_type),show_legend=FALSE),
              top_annotation=HeatmapAnnotation(' '=c,col=list(' '=col_type),show_legend=FALSE),...) )
          print(SumHeatmap2(similarity_matrix,c))
          
      } else if (cluster_num!=0) {
          
          hc=hclust(as.dist(diss_mat),method=hc.method)
          c=cutree(hc,cluster_num)
          hc=dendsort::dendsort(hc)
          or=hc[['labels']][hc[['order']]]
          r=rownames(similarity_matrix)
          c=as.factor(c)
          col_type=col_type[1:max(as.numeric(c))]
          names(col_type)=1:max(as.numeric(c))
          print(  Heatmap(similarity_matrix,cluster_rows=FALSE,cluster_columns=FALSE,row_order=or,column_order=or,
              show_row_names=FALSE,show_column_names=FALSE,name='Similarity\nmatrix',
              row_split=c,column_split=c,
              left_annotation=rowAnnotation(' '=c,col=list(' '=col_type),show_legend=FALSE),
              top_annotation=HeatmapAnnotation(' '=c,col=list(' '=col_type),show_legend=FALSE),...) )
          print(SumHeatmap2(similarity_matrix,c))
          
      }
        
      return(c)
      }    
    }
        
  if (mode=='ConsensusClusterPlus'){
     
     library(ConsensusClusterPlus)
          
     if (select.cutoff){
        ConsensusClustering_result=ConsensusClusterPlus(similarity_matrix,clusterAlg='hc',maxK=maxK,
                       distance='euclidean',innerLinkage=hc.method,finalLinkage=hc.method,title='ConsensusClusteringResult',
                       verbose=FALSE,plot='pdf')
        return(ConsensusClustering_result)
     } else if (cluster_num!=0) {
        ConsensusClustering_result=data
        clustering_result=ConsensusClustering_result[[cluster_num]]
        clustering_matrix=clustering_result[['consensusMatrix']]
        c=clustering_result[['consensusClass']]
             
        col_type=col_type[1:max(as.numeric(c))]
        names(col_type)=1:max(as.numeric(c))
        c=as.factor(c)

        print( Heatmap(clustering_matrix,cluster_rows=FALSE,cluster_columns=FALSE,
              show_row_names=FALSE,show_column_names=FALSE,name='Consensus\nmatrix',
              row_split=c,column_split=c,
              left_annotation=rowAnnotation(' '=c,col=list(' '=col_type),show_legend=FALSE),
              top_annotation=HeatmapAnnotation(' '=c,col=list(' '=col_type),show_legend=FALSE),... ) )

        return (c)   
    }  
  }
        
}
                              
SumHeatmap2=function(similarity_matrix,c){
    
    similarity_matrix_melted=reshape2::melt(similarity_matrix,varnames=c('sample_1','sample_2'),value.name='similarity')
    similarity_matrix_melted=similarity_matrix_melted %>% filter(sample_1!=sample_2)
#     SumHeatmap(similarity_matrix_melted,group.col='group_1',variable.col='group_2',value.col='similarity')
    similarity_matrix_melted$group_1=c[similarity_matrix_melted$sample_1]
    similarity_matrix_melted$group_2=c[similarity_matrix_melted$sample_2]

    df=similarity_matrix_melted
#     p_matrix=lapply(unique(df[,variable.col]),function(x){

#         m=unique(df[,group.col])
#         df=sapply(m,function(y){

#             x1=df[ (df[,group.col] %in% y) & (df[,variable.col] %in% y) ,value.col]
#             x2=df[ (df[,group.col] %in% y) & ((df[,variable.col] %in% x)) ,value.col ]


#             df_test=data.frame(value=c(x1,x2),
#             group=as.factor(c(rep('A',length.out=length(x1)),rep('B',length.out=length(x2)))))
#             colnames(df_test)=c('value','group')
#             df_test=data.frame(df_test)

#             p=coin::pvalue(coin::oneway_test(value~group,df_test))

#             return(p)

#         },USE.NAMES=TRUE) %>% data.frame()

#         rownames(df)=m
#         colnames(df)=x
#         return(df)

#     }) %>% dplyr::bind_cols()

    group.col='group_1'
    variable.col='group_2'
    value.col='similarity'
    
    heatmap_matrix=reshape2::dcast(df,as.formula(paste0(group.col,'~',variable.col)),value.var=value.col,fun.aggregate=mean) %>% 
        data.frame(row.names=1,check.names=FALSE)

#     p_matrix=p_matrix[rownames(heatmap_matrix),colnames(heatmap_matrix)]

#     p.adj.method='bonferroni'
#     p_matrix_adj=unlist(p_matrix) %>% p.adjust(method=p.adj.method) %>% matrix(c(nrow(p_matrix),ncol(p_matrix)))
#     rownames(p_matrix_adj)=rownames(p_matrix)
#     colnames(p_matrix_adj)=colnames(p_matrix)

#     sig.level=c(0.01,0.05)
#     sig.label=c('**','*')

#     sig.level=c(0.01,0.05)
#     sig.label=c('**','*')

    Heatmap(heatmap_matrix,# cluster_rows=FALSE,cluster_columns=FALSE,
#                cell_fun=function(j,i,x,y,w,h,fill){
#         q=min(which(p_matrix_adj[i,j]<=sig.level))
#         if (q<=length(sig.level)){
#             grid.text(sig.label[q],x,y)
#         } else{
#             grid.text('',x,y)
#         }
#     }
           )
}

                              
ClusterCombine=function(c,l,reorder=TRUE){
    
    if ( sum(duplicated( unlist(l) ))!=0 ) {
        stop('Duplicated clusters among combinations!')
    }
    
    for (x in l){
        c[c %in% x]=x[1]
    }
    if (reorder){
        c=factor(c,labels=1:length(unique(c)))
    }
    return(c)
}
                              
# ResetOrder=function(df,by='row'){
    
#     if (by=='row'){
#         o1=apply(df,1,function(x){
#             x=scale(x)
#             which(x==max(x))
#             }) %>% sort()
#         o2=setNames(1:length(o1),names(o1))
#         duplicated_num=o1[duplicated(o1)] %>% unique()
#         for (i in duplicated_num){
#             duplicated_item=o1[i==o1] %>% names()
#             options(warn=-1)
#             duplicated_item=df[duplicated_item,i,drop=FALSE]  %>% t() %>% .[1,] %>% sort(decreasing=TRUE) %>% names()
#             o2[duplicated_item]=sort(o2[duplicated_item])
#         }
#         o2=sort(o2)
#         o2=names(o2)
#         return(df[o2,])
#     }
    
#     if (by=='col'){
#         o1=apply(df,2,function(x){
#             x=scale(x)
#             which(x==max(x))
#         }) %>% sort()
#         o2=setNames(1:length(o1),names(o1))
#         duplicated_num=o1[duplicated(o1)] %>% unique()

#         for (i in duplicated_num){
#             duplicated_item=o1[i==o1] %>% names()
#             options(warn=-1)
#             duplicated_item=df[i,duplicated_item,drop=FALSE]  %>% .[1,] %>% sort(decreasing=TRUE) %>% names()
#             o2[duplicated_item]=sort(o2[duplicated_item])
#         }
#         o2=sort(o2)
#         o2=names(o2)
#         return(df[,o2])
#     }

# }

# Volcano plot for single conditions
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

# Volcano plot for multiple conditions
Draw_MVolcano=function(df,x_col,y_col,gene_col,cluster_col,selected_gene=NULL,yintercept=0,
                       nudge_x=0.8,pnudge_y=0.25,nnudge_y=0,pforce=5,nforce=2.5,nrow=1){
    
    # this function is edited from scRNAtools::markerVolcano, all credits belongs to https://github.com/junjunlab/scRNAtoolVis
    
    library(ggplot2)
    library(ggrepel)
    
    df_selected=df[ df[,gene_col] %in% selected_gene,]
    
    ggplot(df,aes_string(x=x_col,y=y_col))+
        geom_point(color="grey80")+
        geom_hline(yintercept=yintercept,lty="dashed",size=1,color='grey50')+
#         ggrepel::geom_text_repel(data=toppos,aes_string(label=gene_col,color=cluster_col),
#                                  show.legend=FALSE,direction="y",hjust=1,nudge_y=pnudge_y,
#                                  force=pforce,nudge_x=-nudge_x-toppos[,x_col])+
#         ggrepel::geom_text_repel(data=topneg,aes_string(label=gene_col,color=cluster_col),
#                                  show.legend=FALSE,direction="y",hjust=0,nudge_y=nnudge_y,
#                                  force=nforce,nudge_x=nudge_x-topneg[,x_col])+
        ggrepel::geom_text_repel(data=df_selected,
                                 aes_string(label=gene_col,color=cluster_col),
                                 show.legend=FALSE,direction="y",hjust=1,nudge_y=nnudge_y,
                                 force=nforce,# nudge_x=-nudge_x+df[,x_col]
                                 max.overlaps=100
                                )+
        geom_point(data=df_selected,show.legend=FALSE,aes_string(color=cluster_col))+
        facet_grid(eval(parse(text=".~get(cluster_col)")),scale="free")+
        theme_bw(base_size=14)+
        theme(panel.grid=element_blank(),axis.text.x=element_text(angle=45,hjust=1),strip.background=element_rect(color=NA,fill='white')) 
}

                              ForestPlot_single=function(reg_model,point.size=1){
    
    
    coef_int=confint(reg_model)
    colnames(coef_int)=c('lower_limit','higher_limit')
    coef_int=data.frame(coef_int)
    coef_int=na.omit(coef_int)

    coef_int[,'coef']=coef(reg_model)[rownames(coef_int)]
    
    # model-specific settings
    if (inherits(reg_model,'coxph')){
        coef_int=exp(coef_int)
        x_title='Hazard Ratio'
        x_line=1
    } else if (inherits(reg_model,'glm')){
        
        if (reg_model$family$family=='bionomial'){
            
        }
     
        
    }
    
    coef_int[,'variable']=rownames(coef_int)
    coef_int[,'pvalue']=summary(reg_model)$coefficients[rownames(coef_int),'Pr(>|z|)']
    if (0 %in% coef_int[,'pvalue']){
        coef_int[ coef_int[,'pvalue']==0 ,'pvalue']='<2e-16'
    }
    coef_int[,'pvalue']=paste0('P = ',coef_int[,'pvalue'])
    coef_int[,'pvalue']=gsub('= <','< ',coef_int[,'pvalue'])
    
    x_min=min(coef_int[,'lower_limit'],na.rm=TRUE)
    x_max=max(coef_int[,'higher_limit'],na.rm=TRUE)
    coef_int[,'x_pvalue']=x_max+(x_max-x_min)*0.5
    
    # setting coordinates
#     plotrange=c(x_min-(x_max-x_min)*2/3,x_max+(x_max-x_min)*0.15)
    
    p1=ggplot(coef_int)+
        geom_segment(aes(x=lower_limit,xend=higher_limit,y=variable,yend=variable),
                     arrow=arrow(ends='both',angle=90,length=unit(0.1,'cm')) )+
        geom_point(aes(x=coef,y=variable),size=point.size,shape=15)+
        geom_text(aes(x=x_pvalue,y=variable,label=pvalue),hjust=0)+
        xlab(x_title)+
        ylab('')+
        geom_vline(xintercept=x_line,color='red',linetype='dashed')+
        scale_x_log10(breaks=scales::trans_breaks('log10',function(x){10^x}),
                      labels=scales::trans_format('log10',scales::math_format(10^.x))
                      )+
        theme_bw()+
        coord_cartesian(xlim=c(x_min,x_max+ ( (x_max-x_min)*(10^20) )))
    
#     p2=ggplot(coef_int)+
#         geom_label(aes(y=variable,x=0,label=pvalue),)+
#         xlim(c(-0.01,0.01))+
#         theme_classic()+
#         theme()
    
    return(p1)
    
}

# marker expression of ligand and receptor
LinkedPlot=function(df,link_df,x_col,y_col,fill_col,size_col=NULL,
                    color_column=0,
                    widths=c(3,1,3),
                    align='bottom'){
    
    library(ggplot2)
    library(aplot)
    
    colnames(link_df)[1:2]=c('y1_name','y2_name')
    
    df_1=df %>% filter( !!sym(y_col) %in% unique(link_df[,1]) )
    df_2=df %>% filter( !!sym(y_col) %in% unique(link_df[,2]) )

    link_df[,'x1']=2
    link_df[,'x2']=3

    l1=unique(c(link_df[,'y1_name']))
    l2=unique(c(link_df[,'y2_name']))

    ll1=length(l1)-1
    ll2=length(l2)-1

    l=max(length(l1),length(l2))

    coord_1=(0:ll1)/l
    if (align=='bottom'){
        coord_2=(0:ll2)/l
    }
    if (align=='top'){
        coord_2=((0:ll2)-ll2+ll1)/l
    }
    if (align=='center'){
        m=median(0:ll1)
        coord_2=seq(m-ll2/2,m+ll2/2,length.out=ll2+1)/l
    }

    link_df[,'y1_coord']=coord_1[sapply(link_df[,'y1_name'],function(x){which(x==l1)})]
    link_df[,'y2_coord']=coord_2[sapply(link_df[,'y2_name'],function(x){which(x==l2)})]
    link_df[,'col']='1'

    df_1[,'y_coord']=coord_1[sapply(df_1[,y_col],function(x){which(x==l1)})]
    df_2[,'y_coord']=coord_2[sapply(df_2[,y_col],function(x){which(x==l2)})]
    
    if (! is.null(size_col)){
        p1=ggplot(df_1,aes_string(x=x_col,y='y_coord'))+
            geom_point(aes_string(size=size_col,fill=fill_col),shape=21)
        p2=ggplot(df_2,aes_string(x=x_col,y='y_coord'))+
            geom_point(aes_string(size=size_col,fill=fill_col),shape=21)
    } else {
        p1=ggplot(df_1,aes_string(x=x_col,y='y_coord'))+
            geom_tile(aes_string(fill=fill_col))
        p2=ggplot(df_2,aes_string(x=x_col,y='y_coord'))+
            geom_tile(aes_string(fill=fill_col))
    }

    p1=p1+
        scale_y_continuous(position='right',breaks=coord_1,labels=l1,minor_breaks=NULL)+
        # coord_cartesian(ylim=c(0,max(coord_1)))+
        theme_bw()+
        xlab('')+
        ylab('')+
        theme(# legend.position='left',
              axis.title.x=element_blank(),axis.title.y=element_blank())

    p2=p2+
        scale_y_continuous(position='left',breaks=coord_2,labels=l2,minor_breaks=NULL)+
        # coord_cartesian(ylim=c(0,max(coord_2)))+
        theme_bw()+
        xlab('')+
        ylab('')+
        theme(# legend.position='right',
              axis.title.x=element_blank(),axis.title.y=element_blank())
    
    if (color_column==0){
        p3=ggplot(link_df)+
            geom_segment(aes(x=x1,xend=x2,
                             y=y1_coord,yend=y2_coord,color=col))+
            scale_color_manual(values=c('1'='black'))
    }
    
    if (color_column==1) {
        col_1=setNames(rep(ggsci::pal_d3(palette='category10')(10),length.out=length(l1)),
                     l1)
        link_df[,'col']=link_df[,'y1_name']
        
        p1=p1+
            theme(axis.text.y=element_text(color=col_1))
        
        p3=ggplot(link_df)+
            geom_segment(aes(x=x1,xend=x2,
                             y=y1_coord,yend=y2_coord,color=col))+
            scale_color_manual(values=col_1)+
            guides(color='none')
        
    } else if (color_column==2) {
        col_2=setNames(rep(ggsci::pal_d3(palette='category10')(10),length.out=length(l2)),
                     l2)
        link_df[,'col']=link_df[,'y2_name']
        
        p2=p2+
            theme(axis.text.y=element_text(color=col_2))
        
        p3=ggplot(link_df)+
            geom_segment(aes(x=x1,xend=x2,
                             y=y1_coord,yend=y2_coord,color=col))+
            scale_color_manual(values=col_2)+
            guides(color='none')
        
    } else if (! is.null(color_column)){
        link_df[,'col']=link_df[,color_column]
        p3=ggplot(link_df)+
            geom_segment(aes(x=x1,xend=x2,
                             y=y1_coord,yend=y2_coord,color=col))+
            labs(color=color_column)
    }
    
    p3=p3+
        xlab('')+
        ylab('')+
        coord_cartesian(xlim=c(2,3))+
        theme_void()+
        theme(axis.line.x=element_blank(),axis.line.y=element_blank(),
              axis.title.x=element_blank(),axis.title.y=element_blank(),
              axis.text.x=element_blank(),axis.text.y=element_blank(),
              axis.ticks.x=element_blank(),axis.ticks.y=element_blank())
    
    if (ll1>=ll2){
        p=p1 %>% insert_right(p3,width=0.2) %>% insert_right(p2)
    } else if (ll1<ll2) {
        p=p2 %>% insert_left(p3,width=0.2) %>% insert_left(p1)
    }
    
    return(p)
    
}

# AnnotatedPlot - Bar-annotated ggplot
AnnotatedPlot=function(p,
                       df_anno_x=NULL,top_anno_var=NULL,bottom_anno_var=NULL,row_order=NULL,heights=0.1,
                       df_anno_y=NULL,right_anno_var=NULL,left_anno_var=NULL,col_order=NULL,widths=0.1,
                       anno_colors=list()){
    
    library(dplyr)
    library(ggplot2)
    library(RColorBrewer)
    library(aplot)
    
    # reset row order
    if (! is.null(row_order)){
        p_row_order=df_anno_x %>% arrange(!!sym(row_order)) %>% rownames()
        p=p+scale_x_discrete(limits=p_row_order)
    }
    # reset column order
    if (! is.null(col_order)){
        p_col_order=df_anno_y %>% arrange(!!sym(col_order)) %>% rownames()
        p=p+scale_y_discrete(limits=p_col_order)
    }

    # reset order of df_anno_x and df_anno_y
    if (! (is.null(top_anno_var) | is.null(bottom_anno_var))){
        df_anno_x=df_anno_x[ggplot_build(p)$layout$panel_params[[1]]$x$get_labels(),]
    }
    if (! (is.null(left_anno_var) | is.null(right_anno_var))){
        df_anno_y=df_anno_y[ggplot_build(p)$layout$panel_params[[1]]$y$get_labels(),]
    }
    
    j=0
    # annotation on the left
    if (! is.null(left_anno_var)){
        for (i in 1:length(left_anno_var)){
            color_pal=anno_colors[[ left_anno_var[i] ]]
            j=j+1
            assign(paste0('p',j),
                   ggplot(df_anno_y,aes_q(x=0,y=rownames(df_anno_y),fill=df_anno_y[,left_anno_var[i]]))+
                        geom_tile()+
                        theme_classic()+
                        theme(axis.line.x=element_blank(),axis.line.y=element_blank(),
                              axis.title.x=element_blank(),
                              axis.title.y=element_blank(),
                              axis.text.x=element_blank(),axis.text.y=element_blank(),
                              axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+
                        set_color_ggplot(df_anno_y[,left_anno_var[i]],color_pal,legend_name=left_anno_var[i])
                   )

            p=p %>% insert_left(get(paste0('p',j)),width=widths)
        }
    }
    
    # annotation on the right
    if (! is.null(right_anno_var)){
        for (i in 1:length(right_anno_var)){
            color_pal=anno_colors[[ right_anno_var[i] ]]
            j=j+1
            assign(paste0('p',j),
                ggplot(df_anno_y,aes_q(x=0,y=rownames(df_anno_y),fill=df_anno_y[,right_anno_var[i]]))+
                    geom_tile()+
                    theme_classic()+
                    xlab(right_anno_var[i])+
                    theme(axis.line.x=element_blank(),axis.line.y=element_blank(),
                          axis.title.x=element_blank(),
                          axis.title.y=element_blank(),
                          axis.text.x=element_blank(),axis.text.y=element_blank(),
                          axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+
                    set_color_ggplot(df_anno_y[,right_anno_var[i]],color_pal,legend_name=right_anno_var[i])
                )

            p=p %>% insert_right(get(paste0('p',j)),width=widths)
        }
    }
    
    # annotation on the bottom
    if (! is.null(bottom_anno_var)){
        for (i in 1:length(bottom_anno_var)){
            color_pal=anno_colors[[ bottom_anno_var[i] ]]
            
            assign(paste0('p',j),
                   ggplot(df_anno_x,aes_q(y=0,x=rownames(df_anno_x),fill=df_anno_x[,bottom_anno_var[i]]))+
                        geom_tile()+
                        theme_classic()+
                        theme(axis.line.x=element_blank(),axis.line.y=element_blank(),
                              axis.title.x=element_blank(),axis.title.y=element_blank(),
                              axis.text.x=element_blank(),axis.text.y=element_blank(),
                              axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+
                        set_color_ggplot(df_anno_x[,bottom_anno_var[i]],color_pal,legend_name=bottom_anno_var[i])
                  )

            p=p %>% insert_bottom(get(paste0('p',j)),height=widths)
            j=j+1
        }
    }
    
    # annotation on the top
    if (! is.null(top_anno_var)){
        for (i in 1:length(top_anno_var)){
            color_pal=anno_colors[[ top_anno_var[i] ]]
            assign(paste0('p',j),
                   ggplot(df_anno_x,aes_q(y=0,x=rownames(df_anno_x),fill=df_anno_x[,top_anno_var[i]]))+
                        geom_tile()+
                        theme_classic()+
                        theme(axis.line.x=element_blank(),axis.line.y=element_blank(),
                              axis.title.x=element_blank(),axis.title.y=element_blank(),
                              axis.text.x=element_blank(),axis.text.y=element_blank(),
                              axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+
                        set_color_ggplot(df_anno_x[,top_anno_var[i]],color_pal,legend_name=top_anno_var[i])
                   )
            p=p %>% insert_top(get(paste0('p',j)),height=heights)
            j=j+1
        }
    }

    return(p)
}

set_color_ggplot=function(value,color_pal,legend_name){
    if (is.numeric(value)){
        if (! is.null(color_pal)){
            if (! is.colors(color_pal)){
                n=brewer.pal.info[color_pal,][[1]]
                color_pal=brewer.pal(n,color_pal)
            }
        } else {
            color_pal=c("#9a133d", "#b93961", "#d8527c", "#f28aaa", "#f9b4c9", "#f9e0e8", "#ffffff",
            "#eaf3ff", "#c5daf6", "#a1c2ed", "#6996e3", "#4060c8", "#1a318b")
        }
        color_pal=c(color_pal[1],color_pal,color_pal[length(color_pal)])
        p_scale=scale_fill_gradientn(colours=color_pal,
                                     values=scales::rescale(
                                         c(min(value),
                                           seq(quantile(value,0.01),quantile(value,0.99),length.out=length(color_pal)-2),
                                           max(value))
                                     ),name=legend_name)
    } else if (is.character(value)){
        n=length(unique(value))
        if (! is.null(color_pal)){
            if (! is.colors(color_pal)){
                color_pal=brewer.pal(n,color_pal)
            }
        } else {
            color_pal=ggsci::pal_npg(palette='nrc')(n)
        }
        p_scale=scale_fill_manual(values=color_pal,name=legend_name)
    }
    
    return (p_scale)
}

is.colors <- function(x) {
    sapply(x,function(X) {
        tryCatch(is.matrix(col2rgb(X)),error=function(e) FALSE)
    })
}

# OrderedPlot - x-/y- ordered ggplot with or without dendrogram
# Required packages: dplyr, ggplot2, ggtree, aplot
OrderedPlot=function(p,x,y,cluster_value,cluster_var=NULL,
                     cluster_row=FALSE,show_row_dend=FALSE,row_dend_direction='left',row_dend_width=0.1,
                     cluster_column=FALSE,show_column_dend=FALSE,column_dend_direction='top',column_dend_height=0.1,
                     ...){
    
    library(dplyr)
    library(ggplot2)
    library(ggtree)
    library(aplot)
    
    d=p$data
    p1=p
    
    if (cluster_row & (! show_row_dend)){
        if (is.null(cluster_var)){
            d1=reshape2::dcast(d,formula=formula(paste0(y,'~',x)),value.var=cluster_value) %>% data.frame(row.names=1,check.names=FALSE)
        } else {
            d1=reshape2::dcast(d,formula=formula(paste0(y,'~',cluster_var)),value.var=cluster_value) %>% data.frame(row.names=1,check.names=FALSE)
        }
        d1[is.na(d1)]=0
        hc1=hclust(dist(d1),...) %>% dendsort::dendsort()
        hc1.order=hc1[['labels']][hc1[['order']]]
        
        p1=p1+scale_y_discrete(limits=hc1.order)
    }
    
    if (cluster_column & (! show_column_dend)){
        if (is.null(cluster_var)){
            d2=reshape2::dcast(d,formula=formula(paste0(x,'~',y)),value.var=cluster_value) %>% data.frame(row.names=1,check.names=FALSE)
        } else {
            d2=reshape2::dcast(d,formula=formula(paste0(x,'~',cluster_var)),value.var=cluster_value) %>% data.frame(row.names=1,check.names=FALSE)
        }
        d2[is.na(d2)]=0
        hc2=hclust(dist(d2),...) %>% dendsort::dendsort()
        hc2.order=hc2[['labels']][hc2[['order']]]
        
        p1=p1+scale_x_discrete(limits=hc2.order)
    }
    
    if (cluster_row & show_row_dend){
        if (is.null(cluster_var)){
            d1=reshape2::dcast(d,formula=formula(paste0(y,'~',x)),value.var=cluster_value) %>% data.frame(row.names=1,check.names=FALSE)
        } else {
            d1=reshape2::dcast(d,formula=formula(paste0(y,'~',cluster_var)),value.var=cluster_value) %>% data.frame(row.names=1,check.names=FALSE)
        }
        
        d1[is.na(d1)]=0
        hc1=hclust(dist(d1),...) %>% dendsort::dendsort()
        hc1.order=hc1[['labels']][hc1[['order']]]
        
        p2=ggtree(hc1)
        if (row_dend_direction=='left'){
            p1=p1 |> insert_left(p2,width=row_dend_width)
        } else if (row_dend_direction=='right'){
            p1=p1 |> insert_right(p2,width=row_dend_width)
        }
    }
    
    if (cluster_column & show_column_dend){
        if (is.null(cluster_var)){
            d2=reshape2::dcast(d,formula=formula(paste0(x,'~',y)),value.var=cluster_value) %>% data.frame(row.names=1,check.names=FALSE)
        } else {
            d2=reshape2::dcast(d,formula=formula(paste0(x,'~',cluster_var)),value.var=cluster_value) %>% data.frame(row.names=1,check.names=FALSE)
        }
        d2[is.na(d2)]=0
        hc2=hclust(dist(d2),...) %>% dendsort::dendsort()
        hc2.order=hc2[['labels']][hc2[['order']]]
        
        p3=ggtree(hc2)+layout_dendrogram()
        if (column_dend_direction=='top'){
            p1=p1 |> insert_top(p3,height=column_dend_height)
        } else if (column_dend_direction=='right'){
            p1=p1 |> insert_bottom(p3,height=column_dend_height)
        }
    }
    
    return(p1)
    
}

# SunburstPlot - Sunburst 
SunburstPlot=function(df,dims,value=NULL,r=rep(1,length.out=length(dims)),
                      scale=NULL,circular=TRUE,column_keep=NULL,color=NULL){
    
    library(dplyr)
    library(ggplot2)
    library(ggnewscale)

    if (! is.null(value)){
        df=lapply(1:nrow(df),function(x){
            df[x,dims,drop=FALSE] %>% slice(rep(1,df[x,value]))
        }) %>% dplyr::bind_rows()
    }

    # dims=rev(dims)
    df_sta=df[,c(dims,column_keep)]
    dims_ordered=sapply(dims,function(x){df_sta[,x] %>% unique() %>% length()}) %>% 
        sort() %>% names()
    df_sta=df_sta %>% arrange_at(dims_ordered)
    # return(df_sta)
    for (i in dims){
        df_sta[,i]=factor(df_sta[,i],levels=unique(df_sta[,i]))
    }

    r=c(0,r)
    if (is.null(scale)){
        scale=vector('list',length(dims))
    }
    
    if (is.null(color)){
        color=rep(NA,length.out=length(dims))
    }

    p=ggplot()
    for (i in 1:length(dims)){
        df_sta_tmp=df_sta %>% group_by_at(dims[i]) %>% count(name='sunburst_n') %>% data.frame()
        df_sta_tmp[,'sunburst_n']=Reduce(`+`,df_sta_tmp[,'sunburst_n'],accumulate=TRUE)
        df_sta_tmp[,'sunburst_nm']=c(0,df_sta_tmp[-nrow(df_sta_tmp),'sunburst_n'])
        # return(df_tmp)
        assign(paste0('df_sta_',i),df_sta_tmp)
        
        p=p+geom_rect(data=get(paste0('df_sta_',i)),
                      aes_string(xmin=sum(r[1:i]),xmax=sum(r[1:i+1]),ymin='sunburst_nm',ymax='sunburst_n',fill=dims[i]),
                      color=color[i])+
            scale[[i]]+new_scale_fill()
    }
    p=p+theme_void()

    if (circular){
        p=p+coord_polar(theta='y')
    }
    return(p)

}

# layout_circular_community - Network for each community as a circle (like cytoscape)
layout_circular_community=function(community_labels,R=20,k=0.5){
    commnituies_label=unique(community_labels) %>% sort()
    n_communities=length(commnituies_label)
    nodes_per_community=c(table(community_labels))
    
    node_coords=data.frame(
      node=1:length(community_labels),
      community=community_labels,
      x=rep(NA,length.out=length(community_labels)),
      y=rep(NA,length.out=length(community_labels))
    )

    for (i in 1:n_communities) {
        
        theta_center=2*pi*(i-1)/n_communities
        x_center=R*cos(theta_center)
        y_center=R*sin(theta_center)

        m=nodes_per_community[i]
        r=k*m  

        community_nodes=which(community_labels==i)

        for (j in 1:length(community_nodes)) {
            theta_node=2*pi*(j-1)/length(community_nodes)
            x_node=x_center+r*cos(theta_node)
            y_node=y_center+r*sin(theta_node)


            node_coords[ community_nodes[j], "x"]=x_node
            node_coords[ community_nodes[j], "y"]=y_node
            }
        }
    
    return (node_coords[,c('x','y')])
    
}