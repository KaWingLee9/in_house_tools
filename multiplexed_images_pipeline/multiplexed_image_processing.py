# Author: Jiarong Li (Kelvin Lee), Zhixuan Tang, Jiao Yuan, Heqi Wang

# Tutorial
# colors could be used for `pseudo_color`, `plot_pixel`: white, red, green, blue, yellow, magenta, cyan
# `show_legend` in `plot_pixel`: whether to show legend of the markers; if False, the output image has the same size with the raw one
# `show_boundary` in `plot_pixel`: whether to show cell boundary 
# adata=pct2adata(img,mask,channel_names=channel_names,exp_removed=[0,1])
# plot_pixel(adata,color_panel={'red':'HepPar1','green':'Vimentin','blue':'DNA1','cyan':'CD45'})
# plot_cell(adata,tag='cluster')
# pseudo_color(img,color_panel={'white':(30,'CD68'),'red':(17,'CD4'),'green':(19,'CD68'),'blue':(22,'CD8A')})


import scanpy as sc
from gray2color import gray2color
from matplotlib.patches import Patch
import tifffile
import anndata as ad
import cv2
import scipy
from skimage.measure import label, regionprops, regionprops_table
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from skimage.segmentation import find_boundaries

def pct2adata(img,mask,channel_names='',
              exp_removed='',max_quantile=0.98):

    # boundary=find_boundaries(mask,connectivity=1,mode='inner')
    # mask_erode=mask.copy()
    # mask_erode[boundary]=0
    
    # x=list(np.unique(mask_erode))
    # replacement={x[i]:i for i in range(len(x))}
    # replacement.update({i:0 for i in range(mask.max()) if i not in x})
    
    # mask=np.array([[replacement[j] for j in list(i)] for i in mask])
    # mask_erode=np.array([[replacement[j] for j in list(i)] for i in mask_erode])
    
    # mask_erode[mask_erode!=0]=1
    # mask_erode=mask_erode.astype(np.uint8)
    # retval,mask_new,stats,centroids=cv2.connectedComponentsWithStats(mask_erode,connectivity=8)

    # centroids=[cal_centroid([tuple(x) for x in np.swapaxes(np.array(np.where(mask==i)),1,0)]) for i in range(1,mask.max()+1)]
    # centroids=np.array(centroids)

    props_table=regionprops_table(mask,intensity_image=np.moveaxis(img, 0, -1),
                                  properties=["label","area","centroid","axis_major_length",
                                              "axis_minor_length","eccentricity","solidity"])

    centroids=np.array((props_table['centroid-0'],props_table['centroid-1']))
    centroids=centroids.transpose()
    
    # exp_mat=np.zeros((mask.max(),img.shape[0]))
    
    # for g in range(img.shape[0]):
    #     img_g=img[g,:,:]
    #     if img_g.sum()!=0:
    #         img_g[img_g>=np.quantile(img_g[img_g!=0],max_quantile)]=np.quantile(img_g[img_g!=0],max_quantile)
    #     for i in range(1,mask.max()+1):
    #         l= mask==i
    #         exp_mat[i-1,g]=np.mean(img_g[l])

    # for g in range(img.shape[0]):#
        
    #     img_g = img[g, :, :]
    #     print('正在计算的通道数：',g)
    #     if img_g.sum() != 0:
    #         img_g[img_g >= np.quantile(img_g[img_g != 0], max_quantile)] = np.quantile(img_g[img_g != 0],max_quantile)
    #     img_norm[g, :, :] = img_g                                                                            
    # start_time = time.time()
    # with tqdm(total=mask.max(), desc="Processing Cells") as pbar:
    # for i in range(1, mask.max() + 1):       
    #     l = mask == i
    #     true_rows, true_cols = np.where(l)
    #     exp_mat[i - 1, :] = np.mean(img_norm[:, true_rows, true_cols], axis=1) # 均值表达
    #     pbar.update(1)
    # pbar.close()

    object_ids=np.unique(mask[mask != 0])
    exp_mat=np.array([scipy.ndimage.mean(img[i],labels=mask,index=object_ids) for i in range(img.shape[0])]).T
    
    if exp_removed:
        for g in exp_removed:
            exp_mat[:,g]=0
    
    adata=ad.AnnData(exp_mat)
    
    if channel_names:
        adata.var_names=channel_names
    adata.obs_names=range(1,len(adata.obs_names)+1)
    
    adata.uns['img']=img
    adata.uns['mask']=mask
    # adata.uns['size']=img.shape
    adata.uns['var_for_analysis']=[channel_names[i] for i in range(len(channel_names)) if exp_removed.count(i)==0]

    adata.obsm['spatial']=centroids
    adata.obs['area']=props_table['area']
    adata.obs['axis_major_length']=props_table['axis_major_length']
    adata.obs['axis_minor_length']=props_table['axis_minor_length']
    adata.obs['eccentricity']=props_table['eccentricity']
    adata.obs['solidity']=props_table['solidity']
    return(adata)
    
# def cal_centroid(points):
#     x_coords = [p[0] for p in points]
#     y_coords = [p[1] for p in points]
#     _len = len(points)
#     centroid_x = sum(x_coords)/len(points)
#     centroid_y = sum(y_coords)/len(points)
#     return [centroid_x, centroid_y]

# def plot_pixel(adata,color_panel,max_quantile=0.98,show_legend=True,return_mat=False):
    
#     idx=list(adata.var_names)
    
#     if sum([1 for v in color_panel.values() if v in idx])==len(color_panel):
#         channel_num=[idx.index(v) for v in color_panel.values()  if v in idx]
#         color_panel_values=[(i,idx[i]) for i in channel_num]
#         k=list(color_panel.keys())
#         color_panel={k[i]:color_panel_values[i] for i in range(len(color_panel_values))}
    
#     if return_mat:
#         return(pseudo_color(adata,color_panel,max_quantile=max_quantile,return_mat=return_mat))
    
#     pseudo_color(adata,color_panel,max_quantile=max_quantile,return_mat=return_mat)


def _stain(img_gray,color,r,g,b):

    color_map = {
        'blue': (0, 0, img_gray),
        'red': (img_gray, 0, 0),
        'green': (0, img_gray, 0),
        'yellow': (img_gray, img_gray, 0),
        'cyan': (0, img_gray, img_gray),
        'magenta': (img_gray, 0, img_gray),
        'white': (img_gray, img_gray, img_gray)
    }
    r[:, :, 0] += color_map[color][0]
    g[:, :, 0] += color_map[color][1]
    b[:, :, 0] += color_map[color][2]
    return (r, g, b)



def hex_to_rgb(hex_col):
    return tuple(int(hex_col[i:i+2], 16) for i in (1, 3, 5))

# def plot_cell(adata,color='cluster',remove_boundary=True,col=None):
    
#     pic=adata.uns['mask'].copy()
#     mask=adata.uns['mask'].copy()
#     label=adata.obs[color]
    
#     if col==None:
#         color_name=label.cat.categories
#         label=label.cat.rename_categories(dict(zip(color_name,range(len(label)+1,len(label)+len(np.unique(label))+1))))
#         col_default=['#5050FFFF','#CE3D32FF','#749B58FF','#F0E685FF','#466983FF','#BA6338FF','#5DB1DDFF','#802268FF','#6BD76BFF','#D595A7FF','#924822FF',
#                      '#837B8DFF','#C75127FF','#D58F5CFF','#7A65A5FF','#E4AF69FF','#3B1B53FF','#CDDEB7FF','#612A79FF','#AE1F63FF','#E7C76FFF','#5A655EFF',
#                      '#CC9900FF','#99CC00FF','#A9A9A9FF','#CC9900FF','#99CC00FF','#33CC00FF','#00CC33FF','#00CC99FF','#0099CCFF','#0A47FFFF','#4775FFFF',
#                      '#FFC20AFF','#FFD147FF','#990033FF','#991A00FF','#996600FF','#809900FF','#339900FF','#00991AFF','#009966FF','#008099FF','#003399FF',
#                      '#1A0099FF','#660099FF','#990080FF','#D60047FF','#FF1463FF','#00D68FFF','#14FFB1FF']
#         col=col_default[0:(len(np.unique(label)))]
#         col_0=['#FFFFFF']
#         col_0=col_0+col
        
#     else:
#         color_name=list(col.keys())
#         label=label.cat.rename_categories(dict(zip(color_name,range(len(label)+1,len(label)+len(np.unique(label))+1))))
#         col=list(col.values())
#         col_0=['#FFFFFF']+col
    
#     for i in range(0,label.shape[0]):
#         pic[pic==i+1]=label[i]
    
#     pic[pic!=0]=pic[pic!=0]-len(label)
#     colormap=np.array([[hex_to_rgb(i) for i in col_0[0:len(col_0)]]],np.uint8)/255
#     rgb=gray2color(pic.astype(np.uint8), use_pallet=None, custom_pallet=colormap)
    
#     if remove_boundary:
#         boundary=find_boundaries(mask,connectivity=1,mode='inner')
#         pic[boundary>0,:]=255
    
#     legend_elements=[Patch(color=col[i],label=color_name[i]) for i in range(0,len(color_name))]
#     plt.imshow(rgb)
#     plt.axis('off')
#     # plt.invert_yaxis()
#     plt.legend(handles=legend_elements,loc=(1.001,0.5),frameon=False)

def plot_cell_cluster(adata,tag='cluster',col=None):
    
    pic=adata.uns['mask'].copy()
    label=adata.obs[tag]
    label=label.astype('category')
    
    if col==None:
        tag_name=label.cat.categories
        label=label.cat.rename_categories(dict(zip(tag_name,range(len(label)+1,len(label)+len(np.unique(label))+1))))
        col_default=['#5050FFFF','#CE3D32FF','#749B58FF','#F0E685FF','#466983FF','#BA6338FF','#5DB1DDFF','#802268FF','#6BD76BFF','#D595A7FF','#924822FF',
                     '#837B8DFF','#C75127FF','#D58F5CFF','#7A65A5FF','#E4AF69FF','#3B1B53FF','#CDDEB7FF','#612A79FF','#AE1F63FF','#E7C76FFF','#5A655EFF',
                     '#CC9900FF','#99CC00FF','#A9A9A9FF','#CC9900FF','#99CC00FF','#33CC00FF','#00CC33FF','#00CC99FF','#0099CCFF','#0A47FFFF','#4775FFFF',
                     '#FFC20AFF','#FFD147FF','#990033FF','#991A00FF','#996600FF','#809900FF','#339900FF','#00991AFF','#009966FF','#008099FF','#003399FF',
                     '#1A0099FF','#660099FF','#990080FF','#D60047FF','#FF1463FF','#00D68FFF','#14FFB1FF']
        col=col_default[0:(len(np.unique(label)))]
        col_0=['#FFFFFF']
        col_0=col_0+col
        
    else:
        tag_name=list(col.keys())
        label=label.cat.rename_categories(dict(zip(tag_name,range(len(label)+1,len(label)+len(np.unique(label))+1))))
        col=list(col.values())
        col_0=['#FFFFFF']+col
    
    for i in range(0,label.shape[0]):
        pic[pic==i+1]=label[i]
    
    pic[pic!=0]=pic[pic!=0]-len(label)
    colormap=np.array([[hex_to_rgb(i) for i in col_0[0:len(col_0)]]],np.uint8)/255
    rgb=gray2color(pic.astype(np.uint8), use_pallet=None, custom_pallet=colormap)
    legend_elements=[Patch(color=col[i],label=tag_name[i]) for i in range(0,len(tag_name))]
    plt.imshow(rgb)
    plt.axis('off')
    # plt.invert_yaxis()
    plt.legend(handles=legend_elements,loc=(1.001,0.5),frameon=False)


def plot_cell_exp(adata,feature,show_boundary=True,
                  cmap=None,show_colorbar=False,color='red',
                  max_quantile=0.99,min_quantile=0):
    
    mask=adata.uns['mask']
    
    object_exp=np.array(adata[:,feature].copy().X)
    object_exp=np.array([i[0] for i in object_exp])
    
    object_exp[object_exp>=np.quantile(object_exp[object_exp!=0],max_quantile)]=np.quantile(object_exp[object_exp!=0],max_quantile)
    object_exp[ np.array( object_exp<=np.quantile(object_exp[object_exp!=0],min_quantile) ) ]=0
    
    object_exp={(i+1): val for i, val in enumerate(object_exp)}
    
    output_mat=np.zeros_like(mask,dtype=float)
    for obj_id, value in object_exp.items():
        output_mat[mask==obj_id]=value
    
    if show_boundary:  
        boundary=find_boundaries(mask,connectivity=1,mode='inner')
        boundary[boundary>0]=255
        boundary_mat=np.copy(mask)
        m=np.max(boundary_mat)+1
        boundary_mat[boundary]=m
        boundary_mat[boundary_mat!=m]=0
        boundary_mat[boundary_mat!=0]=1
        

    colors=['#000000',color]
    cmap=LinearSegmentedColormap.from_list('custom_cmap',colors,N=500)
        
    plt.imshow(output_mat,cmap=cmap)
    if show_colorbar:
        plt.colorbar()
    plt.imshow(boundary_mat,cmap=LinearSegmentedColormap.from_list('custom_cmap',['#00000000','#FFFFFFFF'],N=100))
    plt.axis('off')
    plt.show()

def plot_pixel(adata,color_panel,
               max_quantile=1,min_quantile=0.01,return_mat=False,show_legend=True,
               show_boundary=False):
    
    idx=list(adata.var_names)
    
    if all(v in idx for v in color_panel.values()):
        channel_num=[idx.index(v) for v in color_panel.values()  if v in idx]
        color_panel_values=[(i,idx[i]) for i in channel_num]
        k=list(color_panel.keys())
        color_panel={k[i]:color_panel_values[i] for i in range(len(color_panel_values))}
        
    img=adata.uns['img']
    # pure black image
    r=np.zeros((img.shape[1],img.shape[2],1))
    g=np.zeros((img.shape[1],img.shape[2],1))
    b=np.zeros((img.shape[1],img.shape[2],1))
    
    color_name=list(color_panel.keys())
    
    if type(list(color_panel.values())[0])==type((0,1)):
        channel_name=[i[1] for i in color_panel.values()]
    if type(list(color_panel.values())[0])==type(0):
        channel_name=[i for i in color_panel.values()]
        
    for color,channel in color_panel.items():
        if type(channel)==type((0,1)):
            channel_num=channel[0]
        if type(channel)==type(0):
            channel_num=channel
        img_gray=img[channel_num,:,:].copy()
        img_gray=np.log2(img_gray+1)
        img_gray[img_gray>=np.quantile(img_gray[img_gray!=0],max_quantile) ]=np.quantile(img_gray[img_gray!=0],max_quantile)
        # img_gray[img_gray>=np.quantile(img_gray,max_quantile) ]=np.quantile(img_gray,max_quantile)
        
        img_gray[ np.array( img_gray<=np.quantile(img_gray[img_gray!=0],min_quantile) ) ]=0
        # img_gray[ np.array( img_gray<=np.quantile(img_gray,min_quantile) ) ]=0
        
        img_gray=img_gray/img_gray.max()*255
        (r,g,b)=_stain(img_gray,color,r,g,b)
        
    rgb_img=np.concatenate([r,g,b],axis=2)
    rgb_img=rgb_img.astype('uint16')
    if return_mat:
        return (rgb_img)

    if show_boundary:  
        mask=adata.uns['mask']
        boundary=find_boundaries(mask,connectivity=1,mode='inner')
        rgb_img=np.copy(rgb_img)
        rgb_img[boundary>0,:]=255
    plt.imshow(rgb_img)
    plt.axis('off')

    if show_legend:
    
        if 'white' in color_name:
            color_name[color_name.index('white')]='black'
        
        legend_elements=[Line2D([0], [0], marker='o',color='black', markerfacecolor='w', label=channel_name[i], markersize=0,linestyle='None')
                for i in range(len(channel_name))]
        plt.legend(handles=legend_elements,loc='lower center',bbox_to_anchor=(0.45,0.98),frameon=False,
                labelcolor=color_name,ncol=min(len(color_name),4),columnspacing=0.5)


def pseudo_color(img,color_panel,max_quantile=1,min_quantile=0.01,
                 show_legend=True,return_mat=False):
    
    import numpy as np
    import pandas as pd
    from matplotlib.lines import Line2D
    import matplotlib.pyplot as plt
    
    # pure black image
    r=np.zeros((img.shape[1],img.shape[2],1))
    g=np.zeros((img.shape[1],img.shape[2],1))
    b=np.zeros((img.shape[1],img.shape[2],1))
    
    def stain(img_gray,color,image_color=[r,g,b]):
        if color=='blue':
            b[:,:,0]=b[:,:,0]+img_gray
        if color=='red':
            r[:,:,0]=r[:,:,0]+img_gray
        if color=='green':
            g[:,:,0]=g[:,:,0]+img_gray
        if color=='yellow':
            r[:,:,0]=r[:,:,0]+img_gray
            g[:,:,0]=g[:,:,0]+img_gray
        if color=='cyan':
            b[:,:,0]=b[:,:,0]+img_gray
            g[:,:,0]=g[:,:,0]+img_gray
        if color=='magenta':
            r[:,:,0]=r[:,:,0]+img_gray
            b[:,:,0]=b[:,:,0]+img_gray
        if color=='white':
            r[:,:,0]=r[:,:,0]+img_gray
            b[:,:,0]=b[:,:,0]+img_gray
            g[:,:,0]=g[:,:,0]+img_gray
        return ((r,g,b))

    color_name=list(color_panel.keys())
    
    if type(list(color_panel.values())[0])==type((0,1)):
        channel_name=[i[1] for i in color_panel.values()]
    if type(list(color_panel.values())[0])==type(0):
        channel_name=[i for i in color_panel.values()]
        
    for color,channel in color_panel.items():
        if type(channel)==type((0,1)):
            channel_num=channel[0]
        if type(channel)==type(0):
            channel_num=channel
        img_gray=img[channel_num,:,:].copy()
        img_gray=np.log2(img_gray+1)

        img_gray[img_gray>=np.quantile(img_gray[img_gray!=0],max_quantile) ]=np.quantile(img_gray[img_gray!=0],max_quantile)
        # img_gray[img_gray>=np.quantile(img_gray,max_quantile) ]=np.quantile(img_gray,max_quantile)
        img_gray[ np.array( img_gray<=np.quantile(img_gray[img_gray!=0],min_quantile) ) ]=0
        # img_gray[ np.array( img_gray<=np.quantile(img_gray,min_quantile) ) ]=0
        
        img_gray=img_gray/img_gray.max()*255
        (r,g,b)=stain(img_gray,color,image_color=(r,g,b))
        
    rgb_img=np.concatenate([r,g,b],axis=2)
    rgb_img=rgb_img.astype('uint16')
    if return_mat:
        return (rgb_img)
        
    img_color=np.concatenate([r,g,b],axis=2)
    img_color=img_color.astype('uint16')
    plt.imshow(img_color)
    plt.axis('off')

    if show_legend:

        if 'white' in color_name:
            color_name[color_name.index('white')]='black'
            
        legend_elements=[Line2D([0], [0], marker='o',color='black', markerfacecolor='w', label=channel_name[i], markersize=0,linestyle='None')
                for i in range(len(channel_name))]
        plt.legend(handles=legend_elements,loc='lower center',bbox_to_anchor=(0.45,0.98),frameon=False,
                labelcolor=color_name,ncol=min(len(color_name),4),columnspacing=0.5)

# detect protein expression in background, cytoplasm and nuclear
def subcellular_exp_qc(img,cell_mask,nuclei_mask,qthres=0.99):

    scale_fun=lambda x: (x-np.mean(x))/np.std(x)
    
    def clip_fun(x,qthres=qthres):
        q=np.quantile(x,qthres)
        x[x>=q]=q
        return(x)

    cell_mask[cell_mask!=0]=1
    nuclei_mask[nuclei_mask!=0]=1

    bkg_mask=cell_mask+nuclei_mask
    bkg_mask[bkg_mask!=0]=3
    bkg_mask[bkg_mask==0]=1
    bkg_mask[bkg_mask==3]=0

    cytoplasm_mask=cell_mask
    cytoplasm_mask[nuclei_mask==1]=0

    img=np.apply_along_axis(clip_fun,0,img)
    img=np.apply_along_axis(scale_fun,0,img)

    b=np.apply_along_axis(np.mean,1,img[:,bkg_mask==1])
    c=np.apply_along_axis(np.mean,1,img[:,cytoplasm_mask==1])
    n=np.apply_along_axis(np.mean,1,img[:,nuclei_mask==1])
    
    df=pd.DataFrame({'Background':b,'Cytoplasm':c,'Nuclei':n})
    df=df.apply(lambda x: (x-np.mean(x))/np.std(x),axis=1)
    return(df)


# scale within samples
# this function is totally copied from https://github.com/scverse/scanpy/issues/2142
def scale_by_batch(adata,batch_key,
                   max_value=10) -> ad.AnnData:
    return ad.concat(
        {
            k: sc.pp.scale(adata[idx],copy=True,max_value=max_value)
            for k, idx in adata.obs.groupby(batch_key).indices.items()
        },
        merge="first"
    )


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


# considering permutations
from joblib import Parallel,delayed
def _CalRelativeDistance(coord_mat_sample,ct1_id,ct2_id,ct3_id):
    
    df_1=coord_mat_sample.loc[ct1_id,:]
    df_2=coord_mat_sample.loc[ct2_id,:]
    df_3=coord_mat_sample.loc[ct3_id,:]
    
    dist_12=np.sqrt(np.sum(( np.array(df_1)[:,np.newaxis,:]- np.array(df_2))**2,axis=2))
    dist_12_min=np.apply_along_axis(np.min,axis=1,arr=dist_12)
    dist_13=np.sqrt(np.sum(( np.array(df_1)[:,np.newaxis,:]- np.array(df_3))**2,axis=2))
    dist_13_min=np.apply_along_axis(np.min,axis=1,arr=dist_13)
    
    df=pd.DataFrame({'CT2':dist_12_min,'CT3':dist_13_min,'ratio':dist_12_min/dist_13_min})
    return(np.mean(df['ratio']))

def _CalRDpermutation(shuffle_type,ct2_id,ct3_id,coord_mat_sample):
    
    if shuffle_type=='shuffle_23':
        x=np.array(list(ct2_id)+list(ct3_id))
        ct2_perm_id=np.random.choice(x,size=len(ct2_id),replace=False)
        ct3_perm_id=np.setdiff1d(x,ct2_perm_id)
    
    if shuffle_type=='shuffle_except1':
        x=np.setdiff1d(np.array(coord_mat_sample.index),ct1_id)
        ct2_perm_id=np.random.choice(x,size=len(ct2_id),replace=False)
        ct3_perm_id=np.random.choice(np.setdiff1d(x,ct2_perm_id),size=len(ct3_id),replace=False)
        
    RD_perm=_CalRelativeDistance(coord_mat_sample,ct2_perm_id,ct2_id,ct3_perm_id)
        
    return(RD_perm)

def relative_distance_analysis(adata,sample_key,cluster_key,CT1,CT2,CT3,
                               shuffle_type='shuffle_23',shuffle_times=100,n_jobs=20):
    
    assert all([item in adata.obs[cluster_key].unique() for item in [CT1,CT2,CT3]]) , 'Cell types do not exisit.'
    
    coord_mat=pd.DataFrame(adata.obsm['spatial'],index=adata.obs_names,columns=['x','y'])
    
    sample_ls=adata.obs[sample_key].unique()
    
    RD_true_dict={}
    RD_zscore_dict={}
    
    for sample_id in sample_ls:
        
        coord_mat_sample=coord_mat[adata.obs[sample_key]==sample_id]
        cell_type_tmp=adata.obs[adata.obs[sample_key]==sample_id][cluster_key]
        
        ct1_id=coord_mat_sample.index[cell_type_tmp==CT1]
        ct2_id=coord_mat_sample.index[cell_type_tmp==CT2]
        ct3_id=coord_mat_sample.index[cell_type_tmp==CT3]
        
        if any([ct1_id.shape[0]==0,ct2_id.shape[0]==0,ct3_id.shape[0]==0]):
            continue
        
        RD_true=_CalRelativeDistance(coord_mat_sample,ct1_id,ct2_id,ct3_id)
        
        RD_perms=Parallel(n_jobs=n_jobs)(delayed(_CalRDpermutation)(shuffle_type=shuffle_type,
                                                                    ct2_id=ct2_id,ct3_id=ct3_id,coord_mat_sample=coord_mat_sample) for i in range(0,shuffle_times))
        RD_zscore=(RD_true-np.mean(RD_perms))/np.std(RD_perms)
        
        RD_true_dict[sample_id]=RD_true
        RD_zscore_dict[sample_id]=RD_zscore
        
    df=pd.DataFrame({'RD':RD_true_dict,'RD_zscore':RD_zscore_dict})
    adata.uns['RelativeDistance']={'CT1':CT1,'CT2':CT2,'CT3':CT3,
                                   'shuffle_type':shuffle_type,'shuffle_times':shuffle_times,
                                   'RD_result':df}
    
    return(adata)