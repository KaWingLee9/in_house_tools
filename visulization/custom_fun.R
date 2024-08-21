# CorPlot - Correlation bubble plot with significance test
# Required packages: ggplot2, Hmisc, dendsort
CorPlot=function(data,cor.method='pearson', # 'pearson', 'spearman'
        size='p.value', # 'p.value', 'p.adj'
        p.adj.method='fdr',
        sig.level=0.05,
        tri='whole', # 'whole', 'lower', 'upper'
        sig.circle=TRUE,stroke=1,maxK=20){
    
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
          geom_point(shape=21,stroke=stroke)
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
                    sig.level=c(0.01,0.05),sig.label=c('**','*'),
                    p.adj=FALSE,p.adj.method='fdr',scale=TRUE,transpose=FALSE,...){
    
    options(warn=-1)
    
    library(dplyr)
    library(ComplexHeatmap)
    
    if (scale){
        df_1=df %>% group_by(!!!syms(variable.col)) %>% mutate(value.col=scale(!!!syms(value.col)))
        df_1[,value.col]=df_1[,'value.col']
        df=data.frame(df_1)
    }
    
    if (test.method=='t.test'){
        agg.fun=mean
    } else {
        agg.fun=median
    }
    
    heatmap_matrix=reshape2::dcast(df,as.formula(paste0(group.col,'~',variable.col)),value.var=value.col,fun.aggregate=agg.fun) %>% 
        data.frame(row.names=1,check.names=FALSE)


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
        
    if (transpose){
        heatmap_matrix=t(heatmap_matrix)
        p_matrix=t(p_matrix)
    }
        
    ht=Heatmap(heatmap_matrix,cell_fun=function(j,i,x,y,w,h,fill){
        q=min(which(p_matrix[i,j]<=sig.level))
        if (q<=length(sig.level)){
            grid.text(sig.label[q],x,y)
        } else{
            grid.text('',x,y)
        }
    },...)
    
    options(warn=1)
    return(ht)

}
                              
# SimilarityHeatmap - Blocks division in similarity heatmap
# Required packages: NbClust, simplifyEnrichment, ComplexHeatmap, ConsensusClusterPlus
SimilarityHeatmap=function(data,mode='automatic',select_cutoff=FALSE,
                           min.nc=2,max.nc=15,cluster_num=0,
                           cutoff_seq=seq(0.6,0.98,by=0.01),cutoff=0.85,
                           maxK=10,...){
        
  library(ComplexHeatmap)

  if (!((mode=='ConsensusClusterPlus') & (select_cutoff==FALSE))) {
      similarity_matrix=cor(t(data))
  }

  col_type=c('#5050FFFF','#CE3D32FF','#749B58FF','#F0E685FF','#466983FF','#BA6338FF','#5DB1DDFF','#802268FF','#6BD76BFF','#D595A7FF','#924822FF',
    '#837B8DFF','#C75127FF','#D58F5CFF','#7A65A5FF','#E4AF69FF','#3B1B53FF','#CDDEB7FF','#612A79FF','#AE1F63FF','#E7C76FFF','#5A655EFF',
    '#CC9900FF','#99CC00FF','#A9A9A9FF','#CC9900FF','#99CC00FF','#33CC00FF','#00CC33FF','#00CC99FF','#0099CCFF','#0A47FFFF','#4775FFFF',
    '#FFC20AFF','#FFD147FF','#990033FF','#991A00FF','#996600FF','#809900FF','#339900FF','#00991AFF','#009966FF','#008099FF','#003399FF',
    '#1A0099FF','#660099FF','#990080FF','#D60047FF','#FF1463FF','#00D68FFF','#14FFB1FF')


  if (mode=='automatic'){
          
    library(simplifyEnrichment)
          
    if (select_cutoff){
      return(select_cutoff(similarity_matrix,cutoff=cutoff_seq,verbose=FALSE,partition_fun=partition_by_hclust))
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
        
  if (mode=='maunal'){
    if (select_cutoff){
        test_index_1=c("kl","ch","hartigan","ccc","scott","marriot","trcovw","tracew","friedman","rubin",
                        "cindex", "db", "silhouette", 
                        "duda", "pseudot2", "beale", 
                        "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", 
                        "gamma", "gplus", "tau", "dunn","sdindex","sdbw")
        
        test_index_2=c("hubert","dindex")
            
        sapply(test_index_1,function(x){
            res=try({
                y=NbClust::NbClust(data=df,distance='euclidean',min.nc=min.nc,max.nc=max.nc, 
                                   method='ward.D2',alphaBeale=0.1,index=x)
            },silent=TRUE)
            if (inherits(res,'try-error')) {return(NULL)}
            return(y$Best.nc[['Number_clusters']])
        })

        sapply(test_index_2,function(x){
            res=try({
                y=NbClust::NbClust(data=df,distance='euclidean',min.nc=min.nc,max.nc=max.nc, 
                                   method='ward.D2',alphaBeale=0.1,index=x)
            },silent=TRUE)
            if (inherits(res,'try-error')) {return(NULL)}
            return(NULL)
        })
            
      return(NULL)

    } else if (cluster_num!=0){
      hc=hclust(dist(df),method='ward.D2')
      hc=dendsort::dendsort(hc)
      c=cutree(hc,cluster_num)
  
      r=rownames(similarity_matrix)
      or=hc[['labels']][hc[['order']]]
      c=as.factor(c)
      col_type=col_type[1:max(as.numeric(c))]
      names(col_type)=1:max(as.numeric(c))
      
      print(  Heatmap(similarity_matrix,cluster_rows=FALSE,cluster_columns=FALSE,row_order=or,column_order=or,
              show_row_names=FALSE,show_column_names=FALSE,name='Similarity\nmatrix',
              row_split=c,column_split=c,
              left_annotation=rowAnnotation(' '=c,col=list(' '=col_type),show_legend=FALSE),
              top_annotation=HeatmapAnnotation(' '=c,col=list(' '=col_type),show_legend=FALSE),...) )
      
      return(c)
    }
  }
        
  if (mode=='ConsensusClusterPlus'){
     
     library(ConsensusClusterPlus)
          
     if (select_cutoff){
        ConsensusClustering_result=ConsensusClusterPlus(similarity_matrix,clusterAlg='hc',maxK=maxK,
                       distance='euclidean',innerLinkage="ward.D2",finalLinkage="ward.D2",title='ConsensusClusteringResult',
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
