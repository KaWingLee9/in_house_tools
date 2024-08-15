# Data preprocessing of multiplexed images data
This pipeline is used for data preprocessing and visualization for multiplexed data, which could be perfectly adapted for scanpy and SOAPy.  
## Cell segmentation using DeepCell - Mesmer
To run Mesmer in docker, please refer to https://deepcell.readthedocs.io/en/master/.  


## Expression quantification, cell type classification and in-situ visualization
Raw multiplexed `.tiff` file and segmentation result `_mask.tiff` are used as input. The analysis pipeline is in `` file.

## Citation
SOAPy: a Python package to dissect spatial architecture, dynamics and communication
