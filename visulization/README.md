# Custom R visualization functions

## Table of Contents
+ [CorPlot - Correlation bubble plot with significance test](#correlation)
+ [SHeatmap - Summarized heatmap of grouped data with significance test](#heatmap1)
+ [SimilarityHeatmap - Blocks division in similarity heatmap](#heatmap2)

## <a name="correlation">CorPlot - Correlation bubble plot with significance test</a>
Required packages: Hmisc, dendsort  
``` r
library(ggplot2)
library(patchwork)
source('https://github.com/KaWingLee9/in_house_tools/blob/main/visulization/custom_fun.R')

# load dataset
data(mtcars)

p1=CorPlot(mtcars,cor.method='pearson',tri='whole',size='p.adj',p.adj.method='fdr',sig.level=0.01,sig.circle=TRUE)
p2=CorPlot(mtcars,cor.method='pearson',tri='upper',size='p.adj',p.adj.method='fdr',sig.level=0.01,sig.circle=FALSE)+expand_limits(x=c(0,12))
p3=CorPlot(mtcars,cor.method='pearson',tri='lower',size='p.value',sig.level=0.01,sig.circle=TRUE)+expand_limits(x=c(-1,11))

options(repr.plot.height=5,repr.plot.width=5*3)
p1+p2+p3
```
<p align="center">
  <img height="400" src="pct/CorPlot.png">
</p>

Parameters of `CorPlot`:
+ `df`: data frame (observations x variables)
+ `cor.method`: methods to calculate correlation, could be `pearson` or `spearman`, passed on to `Hmisc::rcorr(type=...)`
+ `tri`: show the `whole` matrix or `upper`/`lower` triangular matrix
+ `size`: `p.value` or `p.adj` reflected by bubble size
+ `p.adj.method`: method for p value adjustment (e.g. `bonferroni` or `fdr`), passed on to `p.adjust(method=...)`
+ `sig.circle`: whether to show outlines of the bubbles if significant
+ `sig.level`: significance level for hypothesis test, reflected in sig.circle
+ `stroke`: outline thickness of the circle

## <a name="heatmap1">SHeatmap - Summarized heatmap of grouped data with significance test</a>
Required packages: dplyr, ComplexHeatmap  
``` r
# load dataset
data(iris)
head(iris)
  Sepal.Length Sepal.Width Petal.Length Petal.Width Species
1          5.1         3.5          1.4         0.2  setosa
2          4.9         3.0          1.4         0.2  setosa
3          4.7         3.2          1.3         0.2  setosa
4          4.6         3.1          1.5         0.2  setosa
5          5.0         3.6          1.4         0.2  setosa
6          5.4         3.9          1.7         0.4  setosa

# generate a long data frame as input data
df=reshape2::melt(iris,id.vars='Species',measure.vars=c('Sepal.Length','Sepal.Width','Petal.Length','Petal.Width'))
head(df)
  Species     variable value
1  setosa Sepal.Length   5.1
2  setosa Sepal.Length   4.9
3  setosa Sepal.Length   4.7
4  setosa Sepal.Length   4.6
5  setosa Sepal.Length   5.0
6  setosa Sepal.Length   5.4

ht_1=SumHeatmap(df,group.col='Species',variable.col='variable',value.col='value',test.mode='ONEvsVALUE',
                name='Raw\nmatrix',
                test.method='t.test',permutated=FALSE,
                sig.level=c(0.00001,0.01),sig.label=c('**','*'),p.adj=TRUE,scale=FALSE)

ht_2=SumHeatmap(df,group.col='Species',variable.col='variable',value.col='value',test.mode='ONEvsVALUE',
                name='Scaled\nmatrix',col=circlize::colorRamp2(c(-2,0,4),c('#1a318b','#ffffff','#9a133d')),
                test.method='t.test',permutated=FALSE,
                sig.level=c(0.00001,0.01),sig.label=c('**','*'),p.adj=TRUE,scale=TRUE)

options(repr.plot.height=5,repr.plot.width=5*2)
draw(ht_1+ht_2,auto_adjust=FALSE)
```
<p align="center">
  <img height="400" src="pct/SHeatmap.png">
</p>

Parameters of `SHeatmap`:
+ `df`: a long data frame 
+ `group.col`: column name represents sample group
+ `variable.col`: column name represents variable
+ `value.col`: column name represents values
+ `test.mode`: one of  `ONEvsVALUE`, `ONEvsOTHER` and `ONEvsALL`
+ `permutated`: whether to use permutation-based method. Permutation-based methods are from `coin` package
+ `test.method`: one of `t.test`, `wilcox.test` and `oneway.test`
+ `sig.level`, `sig.label`: significance level and corresponding labels, should be __increasing__ ordered
+ `p.adj`, `p.adj.method`: whether to adjust p value and method for p value adjustment (passed on to `p.adjust(method=...)`)
+ `scale`: whether to scale data within each variable. Default: `TRUE`
+ `transpose`: Default:`FALSE`, show heatmap in group x variable. If `TRUE`, show heatmap in variable x group
+ `...`: other arguments passed on to `ComplexHeatmap::Heatmap`

<p align="center">
  <img height="400" src="pct/SHeatmap_description.jpg">
</p>

## <a name="heatmap2">SimilarityHeatmap - Blocks division in similarity heatmap</a>
__Required packages__: simplifyEnrichment, ComplexHeatmap
``` r
library(simplifyEnrichment)
library(ComplexHeatmap)

library(leukemiasEset)
data(leukemiasEset)

df=t(leukemiasEset@assayData$exprs)
df[1:5,1:5]
	ENSG00000000003	ENSG00000000005	ENSG00000000419	ENSG00000000457	ENSG00000000460
GSM330151.CEL	3.386743	3.539030	9.822758	4.747283	3.307188
GSM330153.CEL	3.687029	3.836208	7.969170	4.866344	4.046402
GSM330154.CEL	3.360517	3.246327	9.457491	4.981642	5.529369
GSM330157.CEL	3.459388	3.063286	9.591018	5.982854	4.619444
GSM330171.CEL	3.598589	3.307543	9.863687	5.779449	3.352696
```
The following codes are for samples clustering. You can input transposed expression matrix for genes clustering.  
There are three `mode`s for blocks identification with similarity matrix: `mode`=c(`manual`, `automatic`, `ConsunsusClusterPlus`). For each mode, `select_cutoff=TRUE` presents the effects under different cluster assignment parameters.

__The first mode__: `mode='manual'`; Manually determine cluster number using `ward.D2` hierarchical clustering method. Clstering effects are assessed by `NbClust`.
+ `cluster_num`: cluster number
```r
c1=SimilarityHeatmap(df,mode='manual',cluster_num=4)
```
<p align="center">
  <img height="400" src="pct/SimilarityHeatmap_c1.png">
</p>

__The second mode__: (Default) `mode='automatic'`; Automatically blocks division. Clstering effects are assessed by `simplifyEnrichment`.  
Clustering performance with different cutoff. Lower cutoff -> More clusters.
``` r
SimilarityHeatmap(df,automatic_clustering=TRUE,select_cutoff=TRUE,cutoff_seq=seq(0.5,0.8,by=0.01))
```
<p align="center">
  <img height="400" src="pct/SimilarityHeatmap_select_cutoff.png">
</p>

``` r
c2=SimilarityHeatmap(df,automatic_clustering=TRUE,select_cutoff=FALSE,cutoff=0.52)
```
<p align="center">
  <img height="400" src="pct/SimilarityHeatmap_c2.png">
</p>

__The third mode__: `mode='ConsunsusClusterPlus'`; Blocks identification using `ConsensusClusterPlus`.  

``` r
ConsensusClustering_result=SimilarityHeatmap(df,mode='ConsunsusClusterPlus,select_cutoff=TRUE,cutoff_seq=seq(0.5,0.8,by=0.01))
```

``` r
ConsensusClustering_result=SimilarityHeatmap(df,mode='ConsunsusClusterPlus,select_cutoff=FALSE,cluster_num=11)
```

The function shows the similarity heatmap and returns a vector of the cluster assignment of each sample.
