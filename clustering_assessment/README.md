# Clustering result assessment

## Sensitivity analysis

__Author__: Zhixuan Tang, Jiarong Li  
__Reference__: Petralia F, Ma W, Yaron TM, et al. Pan-cancer proteogenomics characterization of tumor immunity. Cell. 2024;187(5):1255-1277.e27. doi:10.1016/j.cell.2024.01.027
```
Sensitivity analysis to assess the impact of estimation errors in the decomposition results on immune subtype clustering

To evaluate the impact of estimation errors in the decomposition results on the subtypes clustering, we have run additional
computational experiments to assess the robustness of the clustering results. Specifically, we perturbed the cell-type fraction
estimates and proteomic-based pathway scores by adding independent Gaussian noises with varying standard deviations (5%, 10%
and 20% of the original standard deviation). We then evaluated how the clustering results might change based on the perturbed
data matrices. This experiment was repeated 100 times. The Rand index between the original immune subtypes and the clusters
derived from perturbed data matrices were above 0.87 with a median above 0.9 for all the SD levels, which indicates that the
detected immune subtypes are rather robust to the variability in the cell-type fraction estimates.
```
