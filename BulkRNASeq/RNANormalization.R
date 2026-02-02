KeepProteinGene=function(count_mat,species='human'){
    library(dplyr)
    library(biomaRt)
    
    if (species=='human'){mart=useEnsembl(biomart='genes',dataset='hsapiens_gene_ensembl')}
    if (species=='mouse'){mart=useEnsembl(biomart='genes',dataset='mmusculus_gene_ensembl')}
    protein_coding_genes=getBM(attributes = c('ensembl_gene_id', 'external_gene_name','description', 
                                              'chromosome_name','start_position','end_position','strand'),
                               filters='biotype',values='protein_coding',mart=mart)
    protein_coding_genes=unique(protein_coding_genes[,'external_gene_name'])
    g=intersect(rownames(count_mat),protein_coding_genes)
    count_mat=count_mat[g,]
    return(count_mat)
    
}

NormalizeCount=function(count_mat,species='human',method='tpm',length.type='transcript'){
    
    library(dplyr)
    library(biomaRt)
    library(GenomicFeatures)
    
    if (species=='human'){
        txdb=TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
        symbol_table=toTable(org.Hs.eg.db::org.Hs.egSYMBOL)
    } else if (species=='mouse'){
        txdb=TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
        symbol_table=toTable(org.Mm.eg.db::org.Mm.egSYMBOL)
    }
    
    if (length.type=='transcript'){
        transcript.lengths=transcriptLengths(txdb)
        transcript.lengths=merge(transcript.lengths,symbol_table,by='gene_id')
        transcript.lengths=transcript.lengths %>% group_by(symbol) %>% summarise(mean_length=max(tx_len))
        lengths=setNames(transcript.lengths$mean_length,transcript.lengths$symbol)
    }
    
    if (length.type=='exon'){
        exons.list.per.gene=exonsBy(txdb,by="gene")
        exons.lengths=sum(width(GenomicRanges::reduce(exons.list.per.gene)))
        exons.lengths=data.frame(gene_id=names(exons.lengths),length=exons.lengths)
        exons.lengths=merge(exons.lengths,symbol_table,by='gene_id')
        lengths=setNames(exons.lengths$length,exons.lengths$symbol)
    }
    
    lengths=lengths[rownames(count_mat)]
    lengths=lengths[!is.na(lengths)]
    
    count_mat=count_mat[names(lengths),]
    
    if (method=='tpm'){
        tpm_mat=apply(count_mat,2,function(x){
            tpm(x,lengths)
        })
        rownames(tpm_mat)=rownames(count_mat)
        return(tpm_mat)
    }
    
    if (method=='rpkm'){
        rpkm_mat=apply(count_mat,2,function(x){
            rpkm(x,lengths)
        })
        rownames(rpkm_mat)=rownames(count_mat)
        return(rpkm_mat)
    }
    
    
}

# functions to calculate RPKM and TPM were completely copied from https://gist.github.com/slowkow/6e34ccb4d1311b8fe62e
rpkm=function(counts,lengths){
  rate=counts/lengths 
  rate=rate/sum(counts)*1e6
}

tpm=function(counts,lengths){
  rate=counts/lengths
  rate=rate/sum(rate)*1e6
}