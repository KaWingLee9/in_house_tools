# CorPlot - Correlation bubble plot with significance test
# Required packages: ggplot2, Hmisc, dendsort
CorPlot=function(df,cor.method='pearson', # 'pearson', 'spearman'
        size='p.value', # 'p.value', 'p.adj'
        p.adj.method='fdr',
        sig.level=0.05,
        tri='whole', # 'whole', 'lower', 'upper'
        sig.circle=TRUE){
    
    library(ggplot2)

    # cor_result and p_matirx could be replaced
    cor_result=Hmisc::rcorr(as.matrix(df),type=cor.method)
    cor_matrix=cor_result$r
    p_matrix=cor_result$P

    # clustering 
    hc=hclust(as.dist(1-cor_matrix),method='ward.D2')
    hc=dendsort::dendsort(hc)
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
                              
ResetOrder=function(df,by='row'){
    
    if (by=='row'){
        o1=apply(df,1,function(x){
            x=scale(x)
            which(x==max(x))
            }) %>% sort()
        o2=setNames(1:length(o1),names(o1))
        duplicated_num=o1[duplicated(o1)] %>% unique()
        for (i in duplicated_num){
            duplicated_item=o1[i==o1] %>% names()
            options(warn=-1)
            duplicated_item=df[duplicated_item,i,drop=FALSE]  %>% t() %>% .[1,] %>% sort(decreasing=TRUE) %>% names()
            o2[duplicated_item]=sort(o2[duplicated_item])
        }
        o2=sort(o2)
        o2=names(o2)
        return(df[o2,])
    }
    
    if (by=='col'){
        o1=apply(df,2,function(x){
            x=scale(x)
            which(x==max(x))
        }) %>% sort()
        o2=setNames(1:length(o1),names(o1))
        duplicated_num=o1[duplicated(o1)] %>% unique()

        for (i in duplicated_num){
            duplicated_item=o1[i==o1] %>% names()
            options(warn=-1)
            duplicated_item=df[i,duplicated_item,drop=FALSE]  %>% .[1,] %>% sort(decreasing=TRUE) %>% names()
            o2[duplicated_item]=sort(o2[duplicated_item])
        }
        o2=sort(o2)
        o2=names(o2)
        return(df[,o2])
    }

}

# Volcano plot for single conditions
DrawVolcano <- function(deg_result,x='log_fc',FCcutoff=1,
                        y='p_val',pCutoff=0.05,
                        col=c('red2','royalblue','grey30'), #
                        selectLab=NULL){
    library(ggplot2)
    library(ggrepel)
    library(EnhancedVolcano)
    
    keyvals=rep(col[3],length.out=nrow(deg_result))
    names(keyvals)=rep('Unchaged',nrow(deg_result))

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
                    title=NULL,subtitle=NULL
               )

    p=p+scale_color_manual(limits=c('Lower expressed','Unchaged','Higher expressed'),values=c(col[2],col[3],col[1]))

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
