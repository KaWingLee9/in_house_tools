# Some custom visualization in R

```
source('https://github.com/KaWingLee9/in_house_tools/blob/main/visulization/custom_fun.R')
```

### Correlation bubble plot with significance test
Required packages: ggplot2, Hmisc, dendsort  
Parameters of `CorPlot`:
+ df: data frame (observations x variables)
+ cor.method: methods to calculate correlation, could be `pearson` or `spearman`
+ size: size of variable
```
library(patchwork)
```
