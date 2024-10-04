# P-value combination
CombnP=function(pval,ES,min_num=5,min_ratio=0.8,method='z.transform'){
    
    # pval: vector with pvals
    # ES: vector with effect size (with sign)
    # min_num: minimal study number with test result
    # min_ratio: % of the results (ignoring NA) have consistent effect would be considered
    # method: "fisher", "z.transform", "logit". Parameter passing to survcomp::combine.test
    
    # if (method=='z.transform'){
    #     pval[pval==1]=0.99
    # }
    
    pval_with_na=pval
    
    pval_without_na=pval[!is.na(pval)]
    p_pos=pval[(ES>0) & (!is.na(pval))]
    p_neg=pval[(ES<0) & (!is.na(pval))]
    
    if (length(p_pos)>=length(p_neg)){
        con_length=length(p_pos)
        p_sign=p_pos
    }else{
        con_length=length(p_neg)
        p_sign=p_neg
    }
    
    if ((length(pval_without_na)<min_num) | (con_length<length(pval_without_na)*min_ratio)){
        combined_p=1
    }else{
        combined_p=survcomp::combine.test(p_sign,method=method)
    }
        
    return(c('total_study'=length(pval_with_na),'available_study'=length(pval_without_na),
             'sig_study'=con_length,
             'percent_same_dir'=con_length/length(pval_without_na),
             'mean_ES'=mean(ES,na.rm=TRUE),'combined_p'=combined_p))
    # return(c(con_length,length(pval),length(p_pos),length(p_neg)))
}

CombnP_DFLs=function(l,p_col='p_val',ES_col='log_fc',min_num=5,
                     min_ratio=0.8,method='z.transform',ncores=30){
    
    # l: list of data frames
    # p_col: column name of p_val
    # ES_col: column name of ES
    # min_num: minimal study number with result
    # min_ratio: % of the results have consistent effect would be considered
    # method: "fisher", "z.transform", "logit". Parameter passing to survcomp::combine.test
    
    gene_names=unique(unlist(lapply(l,rownames)))

    gene_test_result=parallel::mclapply(gene_names,function(x){
        z=sapply(l,function(y){
            y=data.frame(y)
            y[x,c(p_col,ES_col)] %>% unlist
        }) %>% t() %>% data.frame()
        return(z)
    },mc.cores=ncores)

    names(gene_test_result)=gene_names
    
    test_result=sapply(gene_test_result,function(x){
        CombnP(pval=x[,p_col],ES=x[,ES_col],min_num=min_num,min_ratio=min_ratio,method=method)
    },USE.NAMES=TRUE) %>% t() %>% data.frame()
    
    return(test_result)
    
}

# Rank combination
CombRank=function(ES,up_rank,down_rank,min_num=5,min_ratio=0.8,method='RankProd'){
    
    # ES: vector with effect size (with sign)
    # up_rank, down_rank: increasing and decreasing rank of -log p-value
    # min_num: minimal study number with test result
    # min_ratio: % of the results (ignoring NA) have consistent effect would be considered
    # method: "RankSum", "RankProd"
    
    pos_num=length(which(ES>0))
    neg_num=length(which(ES<0))
    total_num=length(ES)
    na_removed_num=length(which(!is.na(ES)))

    if (pos_num>=neg_num){
        con_length=pos_num
        ES=ES[ES>0]
    }else{
        con_length=neg_num
        ES=ES[ES<0]
    }

    # if ((na_removed_num<10) | (con_length<(na_removed_num*min_ratio))){
    #     rank_c=0
    # }else 
    if (pos_num>=neg_num){
        rank_c=up_rank
    }else if (pos_num<neg_num){
        rank_c=down_rank
    }

    rank_c=rank_c[!is.na(rank_c)]

    if (method=='RankProd'){
        rank_combined=prod(rank_c)^(1/length(rank_c))
    }
    if (method=='RankSum'){
        rank_combined=mean(rank_c)
    }
    
    # percent_same_dir: samples percentage with the same direction of the main trend, divided by samples number of detected result
    return(c('total_study'=total_num,'available_study'=na_removed_num,
             'sig_study'=con_length,
             'percent_same_dir'=con_length/na_removed_num,
             'mean_ES'=mean(ES,na.rm=TRUE),'combined_rank'=rank_combined))
    
}

CombRank_DFLs=function(l,p_col='p_val',ES_col='log_fc',
                       min_num=5,min_ratio=0.8,method='RankProd',
                       ncores=30){
    
    # l: list of data frames
    # p_col: column name of p_val
    # ES_col: column name of ES
    # min_num: minimal study number with test result
    # min_ratio: % of the results (ignoring NA) have consistent effect would be considered
    # method: "RankSum", "RankProd"
    
    gene_names=unique(unlist(lapply(l,rownames)))
    
    rank_ls=parallel::mclapply(l,function(x){
        x[,'logP_n']=(-log(x[,p_col]))*sign(x[,ES_col])
        x[,'up_rank']=rank(x[,'logP_n'])/nrow(x)
        x[,'down_rank']=rank(-x[,'logP_n'])/nrow(x)
        return(x[,c('up_rank','down_rank')])
    },mc.cores=ncores)
    
    test_result=sapply(gene_names,function(y){
        ES=sapply(l,function(x){z=x[y,ES_col]})
        up_rank=sapply(rank_ls,function(x){z=x[y,'up_rank']})
        down_rank=sapply(rank_ls,function(x){z=x[y,'down_rank']})
        CombRank(ES,up_rank,down_rank,min_num=min_num,min_ratio=min_ratio,method=method)
    },USE.NAMES=TRUE) %>% t() %>% data.frame()   
}
