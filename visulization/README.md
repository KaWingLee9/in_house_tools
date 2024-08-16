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
+ `tri`: show the `whole` matrix or `upper`/`lower` triangular matrix;
+ `size`: `p.value` or `p.adj` reflected by bubble size;
+ `p.adj.method`: method for p value adjustment (e.g. `bonferroni` or `fdr`), passed on to `p.adjust(method=...)`;
+ `sig.circle`: whether to show outlines of the bubbles if significant;
+ `sig.level`: significance level for hypothesis test, reflected in sig.circle;

```
library(ggplot2)
library(patchwork)

# load dataset
data(mtcars)

p1=CorPlot(mtcars,cor.method='pearson',tri='whole',size='p.adj',p.adj.method='fdr',sig.level=0.01,sig.circle=TRUE)
p2=CorPlot(mtcars,cor.method='pearson',tri='upper',size='p.adj',p.adj.method='fdr',sig.level=0.01,sig.circle=FALSE)+expand_limits(x=c(0,12))
p3=CorPlot(mtcars,cor.method='pearson',tri='lower',size='p.value',sig.level=0.01,sig.circle=TRUE)+expand_limits(x=c(-1,11))

options(repr.plot.height=5,repr.plot.width=5*3)
p1+p2+p3
```

<p align="center">
  <img height="80" src="pct/CorPlot.png">
</p>
