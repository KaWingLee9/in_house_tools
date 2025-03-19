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

def plot_cell(adata,tag='cluster',col=None):
    
    pic=adata.uns['mask'].copy()
    label=adata.obs[tag]
    
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
    

def plot_pixel(adata,color_panel,
                  max_quantile=0.98,return_mat=False,show_boundary = False):
    
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
        img_gray[img_gray>=np.quantile(img_gray[img_gray!=0],max_quantile)]=np.quantile(img_gray[img_gray!=0],max_quantile)
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
        rgb_img[boundary>0,:] = 255
    plt.imshow(rgb_img)
    plt.axis('off')
    
    if 'white' in color_name:
        color_name[color_name.index('white')]='black'
    
    legend_elements=[Line2D([0], [0], marker='o',color='black', markerfacecolor='w', label=channel_name[i], markersize=0,linestyle='None')
            for i in range(len(channel_name))]
    plt.legend(handles=legend_elements,loc='lower center',bbox_to_anchor=(0.45,0.98),frameon=False,
               labelcolor=color_name,ncol=min(len(color_name),4),columnspacing=0.5)


def pseudo_color(img,color_panel,max_quantile=0.98):
    
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
        img_gray[img_gray>=np.quantile(img_gray[img_gray!=0],max_quantile)]=np.quantile(img_gray[img_gray!=0],max_quantile)
        img_gray=img_gray/img_gray.max()*255
        (r,g,b)=stain(img_gray,color,image_color=(r,g,b))
        
    img_color=np.concatenate([r,g,b],axis=2)
    img_color=img_color.astype('uint16')
    plt.imshow(img_color)
    plt.axis('off')

    if 'white' in color_name:
    color_name[color_name.index('white')]='black'
    
    legend_elements=[Line2D([0], [0], marker='o',color='black', markerfacecolor='w', label=channel_name[i], markersize=0,linestyle='None')
            for i in range(len(channel_name))]
    plt.legend(handles=legend_elements,loc='lower center',bbox_to_anchor=(0.45,0.98),frameon=False,
               labelcolor=color_name,ncol=min(len(color_name),4),columnspacing=0.5)


def subcellular_exp_qc(img,cell_mask,nuclei_mask):

    cell_mask[cell_mask!=0]=1
    nuclei_mask[nuclei_mask!=0]=1

    bkg_mask=cell_mask+nuclei_mask
    bkg_mask[bkg_mask!=0]=3
    bkg_mask[bkg_mask==0]=1
    bkg_mask[bkg_mask==3]=0

    cytoplasm_mask=cell_mask
    cytoplasm_mask[nuclei_mask==1]=0

    b=np.apply_along_axis(np.mean,1,img[:,bkg_mask==1])
    c=np.apply_along_axis(np.mean,1,img[:,cytoplasm_mask==1])
    n=np.apply_along_axis(np.mean,1,img[:,nuclei_mask==1])
    
    df=pd.DataFrame({'Background':b,'Cytoplasm':c,'Nuclei':n})
    return(df)
