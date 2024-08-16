# Some custom visualization methods in R

```
library(ggplot2)
library(patchwork)
source('https://github.com/KaWingLee9/in_house_tools/blob/main/visulization/custom_fun.R')
```

### Correlation bubble plot with significance test
Required packages: Hmisc, dendsort  
Parameters of `CorPlot`:
+ `df`: data frame (observations x variables);
+ `cor.method`: methods to calculate correlation, could be `pearson` or `spearman`;
+ `tri`: show the `whole` matrix or `upper`/`lower` triangular ones;
+ `size`: `p.value` or `p.adj` reflected by bubble size;
+ `p.adj.method`: method for p value adjustment (e.g. `bonferroni` or `fdr`), passed on to `p.adjust(method=...)`;
+ `sig.circle`: whether to show outlines of the bubbles if significant;
+ `sig.level`: significance level for hypothesis test, reflected in sig.circle;

```
library(patchwork)
```
