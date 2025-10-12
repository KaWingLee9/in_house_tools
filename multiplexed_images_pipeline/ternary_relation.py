import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.spatial import distance

# Spatial distribution pattern of three cell types (compare the distance between CT1 and CT2, and that between CT1 and CT3)
# Reference: Immune cell topography predicts response to PD-1 blockade in cutaneous T cell lymphoma
def CalSpatialScore(adata,sample_key,cluster_key,CT1,CT2,CT3):
    
    assert all([item in adata.obs[cluster_key].unique() for item in [CT1,CT2,CT3]]) , 'Cell types do not exisit.'
    
    coord_mat=pd.DataFrame(adata.obsm['spatial'],index=adata.obs_names,columns=['x','y'])
    
    sample_ls=adata.obs[sample_key].unique()
    
    SpatialScore_sample={}
    
    d={}
    d['CT1']=CT1
    d['CT2']=CT3
    d['CT3']=CT3
    
    for sample_id in sample_ls:
        df_1=coord_mat[np.array(adata.obs[sample_key]==sample_id) & np.array(adata.obs[cluster_key]==CT1)]
        df_2=coord_mat[np.array(adata.obs[sample_key]==sample_id) & np.array(adata.obs[cluster_key]==CT2)]
        df_3=coord_mat[np.array(adata.obs[sample_key]==sample_id) & np.array(adata.obs[cluster_key]==CT3)]
        
        if any([df_1.shape[0]==0,df_2.shape[0]==0,df_3.shape[0]==0]):
            continue
        
        dist_12=np.sqrt(np.sum(( np.array(df_1)[:,np.newaxis,:]- np.array(df_2))**2,axis=2))
        dist_12_min=np.apply_along_axis(np.min,axis=1,arr=dist_12)
        dist_13=np.sqrt(np.sum(( np.array(df_1)[:,np.newaxis,:]- np.array(df_3))**2,axis=2))
        dist_13_min=np.apply_along_axis(np.min,axis=1,arr=dist_13)
        
        df=pd.DataFrame({'CT2':dist_12_min,'CT3':dist_13_min,'ratio':dist_12_min/dist_13_min})
        df.index=df_1.index
        
        d[sample_id]=df
        SpatialScore_sample[sample_id]=np.mean(df['ratio'])
    
    d['SpatialScore_sample']=SpatialScore_sample
    adata.uns['Spatialscore']=d
    return(adata)

# Relative distance and relative number considering permutations
from joblib import Parallel,delayed
from scipy.spatial import distance
def _CalRelativeDistance(coord_mat_sample,ct1_id,ct2_id,ct3_id):
    
    df_1=coord_mat_sample.loc[ct1_id,:]
    df_2=coord_mat_sample.loc[ct2_id,:]
    df_3=coord_mat_sample.loc[ct3_id,:]
    
    if ct2_id.shape[0]!=0:
        # dist_12=np.sqrt(np.sum(( np.array(df_1)[:,np.newaxis,:]- np.array(df_2))**2,axis=2))
        dist_12=distance.cdist(df_1,df_2,'euclidean')
        dist_12_min=np.apply_along_axis(np.min,axis=1,arr=dist_12)
    else:
        dist_12_min=np.array([float('inf') for i in range(0,len(ct1_id))])
    
    if ct3_id.shape[0]!=0:
        # dist_13=np.sqrt(np.sum(( np.array(df_1)[:,np.newaxis,:]- np.array(df_3))**2,axis=2))
        dist_13=distance.cdist(df_1,df_3,'euclidean')
        dist_13_min=np.apply_along_axis(np.min,axis=1,arr=dist_13)
    else:
        dist_13_min=np.array([float('inf') for i in range(0,len(ct1_id))])
    
    df=pd.DataFrame({'CT2':dist_12_min,'CT3':dist_13_min,'ratio':dist_12_min/dist_13_min})
    return(np.mean(df['ratio']))

# def _CalRelativeNumber(coord_mat_sample,ct1_id,ct2_id,ct3_id):
    
#     df_1=coord_mat_sample.loc[ct1_id,:]
#     df_2=coord_mat_sample.loc[ct2_id,:]
#     df_3=coord_mat_sample.loc[ct3_id,:]
    
#     dist_12=np.sqrt(np.sum(( np.array(df_1)[:,np.newaxis,:]- np.array(df_2))**2,axis=2))
#     dist_12_min=np.apply_along_axis(np.min,axis=1,arr=dist_12)
#     dist_13=np.sqrt(np.sum(( np.array(df_1)[:,np.newaxis,:]- np.array(df_3))**2,axis=2))
#     dist_13_min=np.apply_along_axis(np.min,axis=1,arr=dist_13)
    
#     df=pd.DataFrame({'CT2':dist_12_min,'CT3':dist_13_min,'ratio':dist_12_min/dist_13_min})
#     return( (np.sum(df['ratio']>1)+1)/(np.sum(df['ratio']<1)+1) )

def _CalRelativeNumber(coord_mat_sample,ct1_id,ct2_id,ct3_id):
    
    df_1=coord_mat_sample.loc[ct1_id,:]
    df_2=coord_mat_sample.loc[ct2_id,:]
    df_3=coord_mat_sample.loc[ct3_id,:]
    
    if ct2_id.shape[0]!=0:
        dist_12=np.sqrt(np.sum(( np.array(df_1)[:,np.newaxis,:]- np.array(df_2) )**2,axis=2))
        dist_12_min=np.apply_along_axis(np.min,axis=1,arr=dist_12)
    else:
        dist_12_min=np.array([0 for i in range(0,len(ct1_id))])
    
    if ct3_id.shape[0]!=0:
        dist_13=np.sqrt(np.sum(( np.array(df_1)[:,np.newaxis,:]- np.array(df_3))**2,axis=2))
        dist_13_min=np.apply_along_axis(np.min,axis=1,arr=dist_13)
    else:
        dist_13_min=np.array([0 for i in range(0,len(ct1_id))])
    
    df=pd.DataFrame({'CT2':dist_12_min,'CT3':dist_13_min,'ratio':dist_12_min/dist_13_min})
    return( [(np.sum(df['ratio']>1)+1)/(np.sum(df['ratio']<1)+1),np.sum(df['ratio']>1),np.sum(df['ratio']<1)] )

def _Calpermutation(shuffle_type,ct1_id,ct2_id,ct3_id,coord_mat_sample,quantatative_method):
    
    if shuffle_type=='shuffle_23':
        x=np.array(list(ct2_id)+list(ct3_id))
        y=np.random.choice(x,size=len(ct2_id)+len(ct3_id),replace=False)
        ct1_perm_id=ct1_id
        ct2_perm_id=y[0:len(ct2_id)]
        ct3_perm_id=y[len(ct2_id):len(y)]
        
    if shuffle_type=='shuffle_123':
        x=np.array(list(ct1_id)+list(ct2_id)+list(ct3_id))
        y=np.random.choice(x,size=len(ct1_id)+len(ct2_id)+len(ct3_id),replace=False)
        ct1_perm_id=y[0:len(ct1_id)]
        ct2_perm_id=y[len(ct1_id):(len(ct1_id)+len(ct2_id))]
        ct3_perm_id=y[(len(ct1_id)+len(ct2_id)):len(y)]     
    
    if shuffle_type=='shuffle_except1':
        x=np.setdiff1d(np.array(coord_mat_sample.index),ct1_id)
        y=np.random.choice(x,size=len(ct2_id)+len(ct3_id),replace=False)
        ct1_perm_id=ct1_id
        ct2_perm_id=y[0:len(ct2_id)]
        ct3_perm_id=y[len(ct2_id):len(y)]
        # ct2_perm_id=np.random.choice(x,size=len(ct2_id),replace=False)
        # ct3_perm_id=np.random.choice(np.setdiff1d(x,ct2_perm_id),size=len(ct3_id),replace=False)
        
    if shuffle_type=='shuffle_all':
        x=np.array(coord_mat_sample.index)
        y=np.random.choice(x,size=len(x),replace=False)
        ct1_perm_id=y[0:len(ct1_id)]
        ct2_perm_id=y[len(ct1_id):(len(ct1_id)+len(ct2_id))]
        ct3_perm_id=y[(len(ct1_id)+len(ct2_id)):len(y)]
        
    if quantatative_method=='RD':
        RD_perm=_CalRelativeDistance(coord_mat_sample,ct1_perm_id,ct2_perm_id,ct3_perm_id)
        return(RD_perm)
    
    if quantatative_method=='RN':
        RN_perm=_CalRelativeNumber(coord_mat_sample,ct1_perm_id,ct2_perm_id,ct3_perm_id)
        return(RN_perm[0])

def relative_distance_analysis(adata,sample_key,cluster_key,CT1,CT2,CT3,
                               shuffle_type='shuffle_except1',quantatative_method='RD',shuffle_times=100,n_jobs=20):
    
    assert all([item in adata.obs[cluster_key].unique() for item in [CT1,CT2,CT3]]) , 'Cell types do not exisit.'
    
    coord_mat=pd.DataFrame(adata.obsm['spatial'],index=adata.obs_names,columns=['x','y'])
    
    sample_ls=adata.obs[sample_key].unique()
    
    RD_true_dict={}
    RD_zscore_dict={}
    
    RN_true_dict={}
    RN_zscore_dict={}
    
    cell_num={}
    CT1_num={}
    CT2_num={}
    CT3_num={}
    
    RN_CT1_2_dict={}
    RN_CT1_3_dict={}
    
    for sample_id in sample_ls:
        
        coord_mat_sample=coord_mat[adata.obs[sample_key]==sample_id]
        cell_type_tmp=adata.obs[adata.obs[sample_key]==sample_id][cluster_key]
        
        ct1_id=coord_mat_sample.index[cell_type_tmp==CT1]
        ct2_id=coord_mat_sample.index[cell_type_tmp==CT2]
        ct3_id=coord_mat_sample.index[cell_type_tmp==CT3]
        
        if any([ct1_id.shape[0]==0,ct2_id.shape[0]==0,ct3_id.shape[0]==0]):
            continue
        
        if quantatative_method=='RD':
            RD_true=_CalRelativeDistance(coord_mat_sample,ct1_id,ct2_id,ct3_id)
            RD_perms=Parallel(n_jobs=n_jobs)(delayed(_Calpermutation)(shuffle_type=shuffle_type,
                                                                      ct1_id=ct1_id,ct2_id=ct2_id,ct3_id=ct3_id,
                                                                      coord_mat_sample=coord_mat_sample,quantatative_method=quantatative_method) 
                                             for i in range(0,shuffle_times))
            RD_zscore=(RD_true-np.mean(RD_perms))/np.std(RD_perms)
            RD_true_dict[sample_id]=RD_true
            RD_zscore_dict[sample_id]=RD_zscore
            cell_num[sample_id]=coord_mat_sample.shape[0]
            CT1_num[sample_id]=ct1_id.shape[0]
            CT2_num[sample_id]=ct2_id.shape[0]
            CT3_num[sample_id]=ct3_id.shape[0]
            
        if quantatative_method=='RN':
            RN_true=_CalRelativeNumber(coord_mat_sample,ct1_id,ct2_id,ct3_id)
            
            RN_true_dict[sample_id]=RN_true[0]
            RN_CT1_2_dict[sample_id]=RN_true[2]
            RN_CT1_3_dict[sample_id]=RN_true[1]
            
            RN_perms=Parallel(n_jobs=n_jobs)(delayed(_Calpermutation)(shuffle_type=shuffle_type,
                                                                      ct1_id=ct1_id,ct2_id=ct2_id,ct3_id=ct3_id,
                                                                      coord_mat_sample=coord_mat_sample,quantatative_method=quantatative_method) 
                                             for i in range(0,shuffle_times))
            # RN_perms=[RN_perms[i][0] for i in range(0,len(RN_perms))]
            RN_zscore=(RN_true[0]-np.mean(RN_perms))/np.std(RN_perms)
            RN_zscore_dict[sample_id]=RN_zscore
            cell_num[sample_id]=coord_mat_sample.shape[0]
            CT1_num[sample_id]=ct1_id.shape[0]
            CT2_num[sample_id]=ct2_id.shape[0]
            CT3_num[sample_id]=ct3_id.shape[0]
            
    # summarize result
    if quantatative_method=='RD':
        df=pd.DataFrame({'Cell_num':cell_num,'CT1_num':CT1_num,'CT2_num':CT2_num,'CT3_num':CT3_num,
                         'RD':RD_true_dict,'RD_zscore':RD_zscore_dict})
        adata.uns['RelativeDistance']={'CT1':CT1,'CT2':CT2,'CT3':CT3,
                                       'shuffle_type':shuffle_type,'quantatative_method':quantatative_method,'shuffle_times':shuffle_times,
                                       'RD_result':df}
    if quantatative_method=='RN':
        df=pd.DataFrame({'Cell_num':cell_num,'CT1_num':CT1_num,'CT2_num':CT2_num,'CT3_num':CT3_num,
                         'CT1_CT2':RN_CT1_2_dict,'CT1_CT3':RN_CT1_3_dict,
                         'RN':RN_true_dict,'RN_zscore':RN_zscore_dict})
        adata.uns['RelativeDistance']={'CT1':CT1,'CT2':CT2,'CT3':CT3,
                                       'shuffle_type':shuffle_type,'quantatative_method':quantatative_method,'shuffle_times':shuffle_times,
                                       'RN_result':df}
        
    return(adata)

# cell annotation of CT1
def subtype_assignment_by_RelativeDistance(adata_x,cluster_key,CT1,CT2,CT3):
    
    coord_mat=pd.DataFrame(adata_x.obsm['spatial'],index=adata_x.obs_names,columns=['x','y'])

    cell_type_tmp=adata_x.obs[cluster_key]

    ct1_id=coord_mat.index[cell_type_tmp==CT1]
    ct2_id=coord_mat.index[cell_type_tmp==CT2]
    ct3_id=coord_mat.index[cell_type_tmp==CT3]

    df_1=coord_mat.loc[ct1_id,:]
    df_2=coord_mat.loc[ct2_id,:]
    df_3=coord_mat.loc[ct3_id,:]


    if ct2_id.shape[0]!=0:
        # dist_12=np.sqrt(np.sum(( np.array(df_1)[:,np.newaxis,:]- np.array(df_2))**2,axis=2))
        dist_12=distance.cdist(df_1,df_2,'euclidean')
        dist_12_min=np.apply_along_axis(np.min,axis=1,arr=dist_12)
    else:
        dist_12_min=np.array([float('inf') for i in range(0,len(ct1_id))])

    if ct3_id.shape[0]!=0:
        # dist_13=np.sqrt(np.sum(( np.array(df_1)[:,np.newaxis,:]- np.array(df_3))**2,axis=2))
        dist_13=distance.cdist(df_1,df_3,'euclidean')
        dist_13_min=np.apply_along_axis(np.min,axis=1,arr=dist_13)
    else:
        dist_13_min=np.array([float('inf') for i in range(0,len(ct1_id))])

    df=pd.DataFrame({'CT2':dist_12_min,'CT3':dist_13_min})
    x=df.apply(lambda x: x[0]<x[1],axis=1)
    df['relative_subtype']=[ 'CT1_CT2' if x[i] else 'CT1_CT3' for i in range(0,len(x)) ]

    df.index=ct1_id

    adata_x.obs['relative_analysis_subtype']=None
    adata_x.obs.loc[ct1_id,'relative_analysis_subtype']=df.loc[ct1_id,'relative_subtype']
    adata_x.obs.loc[ct2_id,'relative_analysis_subtype']=CT2
    adata_x.obs.loc[ct3_id,'relative_analysis_subtype']=CT3
    
    return(adata_x)