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

# Required packages: dplyr, coin, ComplexHeatmap
SumHeatmap=function(df,group_col,variable_col,value_col,test.mode='ONEvsVALUE',
                    test.method='t.test',permutated=FALSE,
                    sig.level=c(0.01,0.05)sig.label=c('**','*'),...){

    heatmap_matrix=reshape2::dcast(df,as.formula(paste0(group_col,'~',variable_col)),value.var=value_col,fun.aggregate=mean) %>% 
        data.frame(row.names=1,check.names=FALSE)


    p_matrix=lapply(unique(df[,variable_col]),function(x){
        df=sapply(unique(df[,group_col]),function(y){

            if (test.mode=='ONEvsVALUE'){
                x1=df[ (df[,group_col] %in% y) & (df[,variable_col] %in% x) ,value_col]
                x2=mean(df[ df[,variable_col] %in% x ,value_col])

                if (test.method=='t.test'){
                    x2=mean(df[ df[,variable_col] %in% x ,value_col])
                    p=t.test(x1,mu=x2)$p.value
                }

                if (test.method=='wilcox.test'){
                    x2=median(df[ df[,variable_col] %in% x ,value_col])
                    p=wilcox.test(x1,mu=x2)$p.value
                }
                
                return(p)
            }

            if (test.mode=='ONEvsOTHER'){
                x1=df[ (df[,group_col] %in% y) & (df[,variable_col] %in% x) ,value_col]
                x2=df[ (!df[,group_col] %in% y) & ((df[,variable_col] %in% x)) ,value_col ]
            }

            if (test.mode=='ONEvsALL'){
                x1=df[ (df[,group_col] %in% y) & (df[,variable_col] %in% x) ,value_col]
                x2=df[ df[,variable_col] %in% x ,value_col]
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

        }) %>% data.frame()
        
        colnames(df)=x
        return(df)
        
    }) %>% dplyr::bind_cols()

    p_matrix=p_matrix[rownames(heatmap_matrix),colnames(heatmap_matrix)]

    sig.label=rev(sig.label)
    sig.level=rev(sig.level)
    
    Heatmap(heatmap_matrix,cell_fun=function(j,i,x,y,w,h,fill){
    for (q in length(sig.level)){
        if (p_matrix[i,j]<=sig.level[q]) {
            grid.text(sig.label[q], x, y)
        }}},...)

}
