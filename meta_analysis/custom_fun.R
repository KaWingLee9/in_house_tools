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
             'mean_ES'=mean(ES,na.rm=TRUE),'combined_rank'=rank_combined,
             'signed_combined_rank'=sign(mean(ES,na.rm=TRUE))*rank_combined))
    
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
    # ncores: number of cores to use
    
    gene_names=unique(unlist(lapply(l,rownames)))
    
    rank_ls=parallel::mclapply(l,function(x){
        x[,'logP_n']=(-log(x[,p_col]))*sign(x[,ES_col])
        x[,'up_rank']=rank(x[,'logP_n'])/nrow(x)
        x[,'down_rank']=rank(-x[,'logP_n'])/nrow(x)
        return(x[,c('up_rank','down_rank')])
    },mc.cores=ncores)
    
    test_result=parallel::mclapply(gene_names,function(y){
        ES=sapply(l_niche,function(x){z=x[y,ES_col]})
        up_rank=sapply(rank_ls,function(x){z=x[y,'up_rank']})
        down_rank=sapply(rank_ls,function(x){z=x[y,'down_rank']})
        CombRank(ES,up_rank,down_rank,min_num=min_num,min_ratio=min_ratio,method=method)
    },mc.cores=ncores) %>% dplyr::bind_rows()
    rownames(test_result)=gene_names
    
    return(test_result)
    
}


# The following codes are from ACAT (https://github.com/yaowuliu/ACAT/)
#'
#' Aggregated Cauchy Assocaition Test
#'
#' A p-value combination method using the Cauchy distribution.
#'
#'
#'
#' @param weights a numeric vector/matrix of non-negative weights for the combined p-values. When it is NULL, the equal weights are used.
#' @param Pvals a numeric vector/matrix of p-values. When it is a matrix, each column of p-values is combined by ACAT.
#' @param is.check logical. Should the validity of \emph{Pvals} (and \emph{weights}) be checked? When the size of \emph{Pvals} is large and one knows \emph{Pvals} is valid, then the checking part can be skipped to save memory.
#' @return The p-value(s) of ACAT.
#' @author Yaowu Liu
#' @examples p.values<-c(2e-02,4e-04,0.2,0.1,0.8);ACAT(Pvals=p.values)
#' @examples ACAT(matrix(runif(1000),ncol=10))
#' @references Liu, Y., & Xie, J. (2019). Cauchy combination test: a powerful test with analytic p-value calculation
#' under arbitrary dependency structures. \emph{Journal of American Statistical Association},115(529), 393-402. (\href{https://amstat.tandfonline.com/doi/abs/10.1080/01621459.2018.1554485}{pub})
#' @export
ACAT<-function(Pvals,weights=NULL,is.check=TRUE){
    Pvals<-as.matrix(Pvals)
    if (is.check){
        #### check if there is NA
        if (sum(is.na(Pvals))>0){
            stop("Cannot have NAs in the p-values!")
        }
        #### check if Pvals are between 0 and 1
        if ((sum(Pvals<0)+sum(Pvals>1))>0){
            stop("P-values must be between 0 and 1!")
        }
        #### check if there are pvals that are either exactly 0 or 1.
        is.zero<-(colSums(Pvals==0)>=1)
        is.one<-(colSums(Pvals==1)>=1)
        if (sum((is.zero+is.one)==2)>0){
            stop("Cannot have both 0 and 1 p-values in the same column!")
        }

        if (sum(is.zero)>0){
            warning("There are p-values that are exactly 0!")
        }
        if (sum(is.one)>0){
            warning("There are p-values that are exactly 1!")
        }

    }
    #### Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
    if (is.null(weights)){
        is.weights.null<-TRUE
    }else{
        is.weights.null<-FALSE
        weights<-as.matrix(weights)
        if (sum(dim(weights)!=dim(Pvals))>0){
            stop("The dimensions of weights and Pvals must be the same!")
        }else if (is.check & (sum(weights<0)>0)){
            stop("All the weights must be nonnegative!")
        }else{
            w.sum<-colSums(weights)
            if (sum(w.sum<=0)>0){
                stop("At least one weight should be positive in each column!")
            }else{
                for (j in 1:ncol(weights)){
                    weights[,j]<-weights[,j]/w.sum[j]
                }
            }
        }

    }

    #### check if there are very small non-zero p values and calcuate the cauchy statistics
    is.small<-(Pvals<1e-15)
    if (is.weights.null){
         Pvals[!is.small]<-tan((0.5-Pvals[!is.small])*pi)
         Pvals[is.small]<-1/Pvals[is.small]/pi
         cct.stat<-colMeans(Pvals)
    }else{
         Pvals[!is.small]<-weights[!is.small]*tan((0.5-Pvals[!is.small])*pi)
         Pvals[is.small]<-(weights[is.small]/Pvals[is.small])/pi
         cct.stat<-colSums(Pvals)
    }
    #### return the ACAT p value(s).
    pval<-pcauchy(cct.stat,lower.tail = F)
    return(pval)
}

#'
#' A set-based test that uses ACAT to combine the variant-level p-values.
#'
#'
#' @param G a numeric matrix or dgCMatrix with each row as a different individual and each column as a separate gene/snp. Each genotype should be coded as 0, 1, 2.
#' @param obj an output object of the \code{\link{NULL_Model}} function.
#' @param weights.beta a numeric vector of parameters for the beta weights for the weighted kernels. If you want to use your own weights, please use the “weights” parameter. It will be ignored if “weights” parameter is not null.
#' @param weights a numeric vector of weights for the SNP p-values. When it is NULL, the beta weight with the “weights.beta” parameter is used.
#' @param mac.thresh a threshold of the minor allele count (MAC). The Burden test will be used to aggregate the SNPs with MAC less than this thrshold.
#' @return The p-value of ACAT-V.
#' @details The Burden test is first used to aggregate very rare variants with Minor Allele Count (MAC) < \emph{mac.thresh} (e.g., 10), and a Burden p-value is obtained. For each of the variants with MAC >= \emph{mac.thresh}, a variant-level p-value is calculated. Then, ACAT is used to combine the variant-level p-values and the Burden test p-value of very rare variants.
#'
#' If \emph{weights.beta} is used, then the weight for the Burden test p-value is demetermined by the average Minor Allele Frequency (MAF) of the variants with MAC < \emph{mac.thresh}; if the user-specified \emph{weights} is used, then the weight for the Burden test p-value is the average of \emph{weights} of the variants with MAC < \emph{mac.thresh}.
#'
#' Note that the \emph{weights} here are for the SNP p-vlaues. In SKAT, the weights are for the SNP score test statistics. To transfrom the SKAT weights to the \emph{weights} here, one can use the formula that \emph{weights} = (skat_weights)^2*MAF*(1-MAF).
#'
#' @author Yaowu Liu
#' @references Liu, Y., et al. (2019). ACAT: A fast and powerful p value combination
#' method for rare-variant analysis in sequencing studies.
#' \emph{American Journal of Humann Genetics 104}(3), 410-421.
#' (\href{https://www.sciencedirect.com/science/article/pii/S0002929719300023}{pub})
#' @examples  library(Matrix)
#' @examples  data(Geno)
#' @examples  G<-Geno[,1:100] # Geno is a dgCMatrix of genotypes
#' @examples  Y<-rnorm(nrow(G)); Z<-matrix(rnorm(nrow(G)*4),ncol=4)
#' @examples  obj<-NULL_Model(Y,Z)
#' @examples  ACAT_V(G,obj)
#' @export
ACAT_V<-function(G,obj,weights.beta=c(1,25),weights=NULL,mac.thresh=10){
    ### check weights
    if (!is.null(weights) && length(weights)!=ncol(G)){
        stop("The length of weights must equal to the number of variants!")
    }

    mac<-Matrix::colSums(G)
    ### remove SNPs with mac=0
    if (sum(mac==0)>0){
        G<-G[,mac>0,drop=FALSE]
        weights<-weights[mac>0]
        mac<-mac[mac>0]
        if (length(mac)==0){
            stop("The genotype matrix do not have non-zero element!")
        }
    }
    ### p and n
    p<-length(mac)
    n<-nrow(G)
    ###

    if (sum(mac>mac.thresh)==0){  ## only Burden
        pval<-Burden(G,obj, weights.beta = weights.beta, weights = weights)
    }else if (sum(mac<=mac.thresh)==0){ ## only cauchy method
        if (is.null(weights)){
            MAF<-mac/(2*n)
            W<-(dbeta(MAF,weights.beta[1],weights.beta[2])/dbeta(MAF,0.5,0.5))^2
        }else{
            W<-weights
        }

        Mpvals<-Get.marginal.pval(G,obj)
        pval<-ACAT(Mpvals,W)
    }else{  ## Burden + Cauchy method
        is.very.rare<-mac<=mac.thresh
        weights.sparse<-weights[is.very.rare]
        weights.dense<-weights[!is.very.rare]
        pval.dense<-Burden(G[,is.very.rare,drop=FALSE],obj, weights.beta = weights.beta, weights = weights.sparse)

        Mpvals<-Get.marginal.pval(G[,!is.very.rare,drop=FALSE],obj)

        Mpvals<-c(Mpvals,pval.dense)
        if (is.null(weights)){
            MAF<-mac/(2*n)
            mafs<-c(MAF[!is.very.rare],mean(MAF[is.very.rare])) ## maf for p-values
            W<-(dbeta(mafs,weights.beta[1],weights.beta[2])/dbeta(mafs,0.5,0.5))^2
        }else{
            W<-c(weights.dense,mean(weights.sparse))
        }


        is.keep<-rep(T,length(Mpvals))
        is.keep[which(Mpvals==1)]<-F  ## remove p-values of 1.
        pval<-ACAT(Mpvals[is.keep],W[is.keep])
    }
    return(pval)
}

#'
#'
#' Get parameters and residuals from the NULL model
#'
#' Compute model parameters and residuals for ACAT-V
#'
#'
#' @param Y a numeric vector of outcome phenotypes.
#' @param Z a numeric matrix of covariates. Z must be full-rank. Do not include intercept in Z. The intercept will be added automatically.
#' @return This function returns an object that has model parameters and residuals of the NULL model of no association between genetic variables and outcome phenotypes. After obtaining it, please use \code{\link{ACAT_V}} function to conduct the association test.
#' @details \emph{Y} could only be continuous or binary. If \emph{Y} is continuous, a linear regression model is fitted. If \emph{Y} is binary, it must be coded as 0,1 and a logistic model is fitted.
#' @author Yaowu Liu
#' @examples  Y<-rnorm(10000)
#' @examples  Z<-matrix(rnorm(10000*4),ncol=4)
#' @examples  obj<-NULL_Model(Y,Z)
#' @export
NULL_Model<-function(Y,Z=NULL){
    n<-length(Y)
    #### check the type of Y
    if ((sum(Y==0)+sum(Y==1))==n){
        out_type<-"D"
    }else{
        out_type<-"C"
    }
    #### Add intercept
    Z.tilde<-cbind(rep(1,length(Y)),Z)
    if (out_type=="C"){
        #### estimate of sigma square
        Z.med<-Z.tilde%*%solve(chol(t(Z.tilde)%*%Z.tilde))   ## Z.med%*%t(Z.med) is the projection matrix of Z.tilde
        Y.res<-as.vector(Y-(Y%*%Z.med)%*%t(Z.med))
        sigma2<-sum(Y.res^2)/(n-ncol(Z.med))
        #### output
        res<-list()
        res[["out_type"]]<-out_type
        res[["Z.med"]]<-Z.med
        res[["Y.res"]]<-Y.res
        res[["sigma2"]]<-sigma2
    }else if (out_type=="D"){
        #### fit null model
        g<-glm(Y~0+Z.tilde,family = "binomial")
        prob.est<-g[["fitted.values"]]
        #### unstandarized residuals
        Y.res<-(Y-prob.est)
        ### Sigma when rho=0
        sigma2.Y<-prob.est*(1-prob.est)  ### variance of each Y_i
        ### output
        res<-list()
        res[["out_type"]]<-out_type
        res[["Z.tilde"]]<-Z.tilde
        res[["Y.res"]]<-Y.res
        res[["sigma2.Y"]]<-sigma2.Y
    }
    return(res)
}




Get.marginal.pval<-function(G,obj){
    ### check obj
    if (names(obj)[1]!="out_type"){
        stop("obj is not calculated from MOAT_NULL_MODEL!")
    }else{
        out_type<-obj[["out_type"]]
        if (out_type=="C"){
            if (!all.equal(names(obj)[2:length(obj)],c("Z.med","Y.res","sigma2"))){
                stop("obj is not calculated from MOAT_NULL_MODEL!")
            }else{
                Z.med<-obj[["Z.med"]]
                Y.res<-obj[["Y.res"]]
                n<-length(Y.res)
                SST<-obj[["sigma2"]]*(n-ncol(Z.med))
            }
        }else if (out_type=="D"){
            if (!all.equal(names(obj)[2:length(obj)],c("Z.tilde","Y.res","sigma2.Y"))){
                stop("obj is not calculated from MOAT_NULL_MODEL!")
            }else{
                Z.tilde<-obj[["Z.tilde"]]
                Y.res<-obj[["Y.res"]]
                sigma2.Y<-obj[["sigma2.Y"]]
                n<-length(Y.res)
            }
        }
    }

    if (class(G)!="matrix" && class(G)!="dgCMatrix"){
        stop("The class of G must be matrix or dgCMatrix!")
    }

    if (out_type=="C"){
        G_tX.med<-as.matrix(Matrix::crossprod(Z.med,G))
        ### Sigma^2 of G
        Sigma2.G<-Matrix::colSums(G^2)-Matrix::colSums(G_tX.med^2)
        SSR<-as.vector((Y.res%*%G)^2/Sigma2.G)
        SSR[Sigma2.G<=0]<-0
        df.2<-n-1-ncol(Z.med)
        t.stat<-suppressWarnings(sqrt(SSR/((SST-SSR)/df.2)))
        marginal.pvals<-2*pt(t.stat,(n-1-ncol(Z.med)),lower.tail = FALSE)
    }else if (out_type=="D"){
        Z.stat0<-as.vector(Y.res%*%G)
        ### Sigma when rho=0
        tG_X.tilde_sigma2<-as.matrix(Matrix::crossprod(G,Z.tilde*sigma2.Y))
        Sigma2.G<-Matrix::colSums(G^2*sigma2.Y)-diag(tG_X.tilde_sigma2%*%solve(t(Z.tilde)%*%(Z.tilde*sigma2.Y))%*%t(tG_X.tilde_sigma2))
        marginal.pvals<-2*pnorm(abs(Z.stat0)/sqrt(Sigma2.G),lower.tail = FALSE)
    }

    return(marginal.pvals)
}


Burden<-function(G,obj,kernel="linear.weighted",weights.beta=c(1,25),weights=NULL){
    ### check obj
    if (names(obj)[1]!="out_type"){
        stop("obj is not calculated from NULL_MODEL!")
    }else{
        out_type<-obj[["out_type"]]
        if (out_type=="C"){
            if (!all.equal(names(obj)[2:length(obj)],c("Z.med","Y.res","sigma2"))){
                stop("obj is not calculated from NULL_MODEL!")
            }else{
                Z.med<-obj[["Z.med"]]
                Y.res<-obj[["Y.res"]]/sqrt(obj[["sigma2"]])  ## rescaled residules
                n<-length(Y.res)
            }
        }else if (out_type=="D"){
            if (!all.equal(names(obj)[2:length(obj)],c("Z.tilde","Y.res","sigma2.Y"))){
                stop("obj is not calculated from NULL_MODEL!")
            }else{
                Z.tilde<-obj[["Z.tilde"]]
                Y.res<-obj[["Y.res"]]
                sigma2.Y<-obj[["sigma2.Y"]]
                n<-length(Y.res)
            }
        }
    }
    ### MAF
    MAF<-Matrix::colSums(G)/(2*dim(G)[1])
    p<-length(MAF)
    #### weights
    if (kernel=="linear.weighted"){
        if (is.null(weights)){
            W<-dbeta(MAF,weights.beta[1],weights.beta[2])
        }else{
            if (length(weights)==p){
                W<-weights
            }else{
                stop("The length of weights must equal to the number of variants!")
            }
        }

    }else if (kernel=="linear"){
        W<-rep(1,p)
    }else{
        stop("The kernel name is not valid!")
    }

    ###### if G is sparse or not
    if (class(G)=="matrix" || class(G)=="dgCMatrix"){
        if (out_type=="C"){
            Z.stat.sum<-as.vector((Y.res%*%G)%*%W)
            Gw<-G%*%W
            sigma.z<-sqrt(sum(Gw^2)-sum((t(Z.med)%*%(Gw))^2))
        }else if (out_type=="D"){
            Z.stat.sum<-as.vector((Y.res%*%G)%*%W)
            Gw<-as.vector(G%*%W)
            sigma.z<-sum(Gw^2*sigma2.Y)-((Gw*sigma2.Y)%*%Z.tilde)%*%solve(t(Z.tilde)%*%(Z.tilde*sigma2.Y))%*%t((Gw*sigma2.Y)%*%Z.tilde)
            sigma.z<-as.vector(sqrt(sigma.z))
        }
    }else{
        stop("The class of G must be matrix or dgCMatrix!")
    }

    V<-Z.stat.sum/sigma.z
    Q<-V^2   ## Q test statistic
    pval<-1-pchisq(Q,df=1)
    return(pval)
}
