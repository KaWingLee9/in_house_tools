# Data preprocessing of multiplexed imaging data (Python)
This pipeline is used for data preprocessing and visualization for multiplexed imaging data, which could be perfectly adapted for scanpy and SOAPy.  
## Cell segmentation using DeepCell - Mesmer
To run Mesmer in docker, please refer to https://deepcell.readthedocs.io/en/master/.  
Required packages: numpy, pandas, matplotlib, tifffile, deepcell
``` python
import os
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import tifffile
from PIL import Image
from deepcell.utils.plot_utils import create_rgb_image, make_outline_overlay
from deepcell_toolbox import erode_edges
# from skimage.exposure import rescale_intensity
# from skimage.measure import label, regionprops, regionprops_table

from mesmer_function import *
from deepcell.applications import Mesmer
app = Mesmer()

# The order of the image and the marker_list must be the same
marker_list = ['Au','Background','Beta_catenin','Ca','CD11b','CD11c','CD138','CD16','CD20','CD209','CD3',
 'CD31','CD4','CD45','CD45RO','CD56','CD63','CD68','CD8','dsDNA','EGFR','Fe','FoxP3','H3K27me3',
 'H3K9ac','HLA-DR','HLA-I','IDO','CK17','CK6','Ki67','Lag3','MPO','Na','P','p53','PanCK','PD-L1',
 'PD-1','pS6','Si','SMA','Ta','Vimentin']
# the keys of the dict can be whatever you like
in_nuc_marker_kwargs = {'nuc_channel': ['DNA1']}
in_mem_marker_kwargs = {'mem_channel_1': ['PanCK','CK6','CK17','EGFR','Beta_catenin'], 
                        'mem_channel_3':['SMA','Vimentin','CD31'],
                        'mem_channel_2': ['CD45','CD3', 'CD68',  'CD8A', 'MPO']}

file='./Point4.tiff'
img_multipledxed=tifffile.imread(file)
# remove the highest and lowest 2% of the pixels
img_norm=Denoise_img(img_multipledxed,98)
img_mesmer=Generate_mesmer_input(img_norm,marker_list,in_nuc_marker_kwargs,in_mem_marker_kwargs)
segmentation_predictions=app.predict(img_mesmer,preprocess_kwargs={'percentile':99.9},image_mpp=0.5,compartment='whole-cell')
# save mask file
tifffile.imwrite(file.split('.')[0]+'_mask.tiff', segmentation_predictions[0,:,:,0])
```

Note: Parameter `iamge_mpp` in `app.predict` specifies the resolution of the image (aka $\mu m$ of each pixel).  
IMC: 1, MIBI_TOF: 0.5
## Expression quantification, cell type classification and in-situ visualization
Raw multiplexed `.tiff` file and segmentation result `_mask.tiff` are used as input. The codes are in `.ipynb`.  
SOAPy could be used for further spatial-related analysis.

## Data generation from different plaforms
### IMC (Imaging mass cytometry)
Convert IMC `.txt` files to multiplexed `.tiff` files.  
Reaquired packages: readimc
``` python
import os
import numpy as np
import tifffile
from readimc import TXTFile

for file in file_ls:
    with TXTFile(file) as f:
        img=f.read_acquisition()
    tifffile.imwrite(file.split('.')[0]+'.tiff',data=img)
```

Reference: Windhager, J., Zanotelli, V.R.T., Schulz, D. et al. An end-to-end workflow for multiplexed image processing and analysis. Nat Protoc 18, 3565–3613 (2023). https://doi.org/10.1038/s41596-023-00881-0

## Citation
Wang H, Li J, Jing S, et al. SOAPy: a Python package to dissect spatial architecture, dynamics and communication[J]. bioRxiv, 2023: 2023.12. 21.572725.
