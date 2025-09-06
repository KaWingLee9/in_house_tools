# p-value rank combination
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
        rank_combined=median(rank_c)
    }
    
    # percent_same_dir: samples percentage with the same direction of the main trend, divided by samples number of detected result
    return(c('total_study'=total_num,'available_study'=na_removed_num,
             # 'sig_study'=con_length,
             # 'percent_same_dir'=con_length/na_removed_num,
             'median_ES'=median(ES,na.rm=TRUE),'combined_rank'=rank_combined,
             'signed_combined_rank'=sign(median(ES,na.rm=TRUE))*rank_combined))
    
}

CombRank_DFLs=function(l,p_col='p_val',ES_col='log_fc',
                       min_num=5,min_ratio=0.8,method='RankProd',
                       ncores=30# ,quantile_est=1/3
                       ){
    
    # l: list of data frames
    # p_col: column name of p_val
    # ES_col: column name of ES
    # min_num: minimal study number with test result
    # min_ratio: % of the results (ignoring NA) have consistent effect would be considered
    # method: "RankSum", "RankProd"
    # # ncores: number of cores to use
    
    gene_names=unique(unlist(lapply(l,rownames)))
    
    rank_ls=parallel::mclapply(l,function(x){
        x[,'logP_n']=(-log(x[,p_col]))*sign(x[,ES_col])
        x[,'up_rank']=rank(x[,'logP_n'])/nrow(x)
        x[,'down_rank']=rank(-x[,'logP_n'])/nrow(x)
        return(x[,c('up_rank','down_rank')])
    },mc.cores=ncores)
    
    test_result=parallel::mclapply(gene_names,function(y){
        ES=sapply(l,function(x){z=x[y,ES_col]})
        up_rank=sapply(rank_ls,function(x){z=x[y,'up_rank']})
        down_rank=sapply(rank_ls,function(x){z=x[y,'down_rank']})
        CombRank(ES,up_rank,down_rank,min_num=min_num,min_ratio=min_ratio,method=method)
    },mc.cores=ncores) %>% dplyr::bind_rows()

    # test_result[,'quantile_pval']=parallel::mclapply(gene_names,function(y){

    #     s=sapply(l,function(x) { x[y,ES_col] })

    #     if (sum(sign(s),na.rm=TRUE)>=0){
    #         quantile_est=1-quantile_est
    #     }
    #     k=quantile( sapply(l,function(x){ z= -log2(x[y,p_col]) * sign(x[y,ES_col]) }), probs=quantile_est,na.rm=TRUE)
    #     k=2^(-sign(k)*k)*sign(k)
        
    # },mc.cores=ncores)  %>% unlist()

    test_result[,'sig_study_number']=parallel::mclapply(gene_names,function(y){

        s=sapply(l,function(x) { x[y,ES_col] })

        if (sum(sign(s),na.rm=TRUE)>=0){
            k=sum(sapply(l,function(x){ (x[y,p_col]<=0.05) & (x[y,ES_col]>=0) }),na.rm=TRUE)
        } else {
            k=sum(sapply(l,function(x){ (x[y,p_col]<=0.05) & (x[y,ES_col]<=0) }),na.rm=TRUE)
        }

        return(k)
        
    },mc.cores=ncores)  %>% unlist()

    test_result[,'sig_study_ratio']=test_result[,'sig_study_number']/test_result[,'available_study']

    # test_result[,'quantile_ES']=parallel::mclapply(gene_names,function(y){
    #     quantile( sapply(l,function(x){z=x[y,ES_col]}), probs=quantile_est,na.rm=TRUE)
    # },mc.cores=ncores)  %>% unlist()

    test_result=data.frame(test_result)
    rownames(test_result)=gene_names
    return(test_result)
    
}
