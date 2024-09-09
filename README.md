# iDAS

Interpretable Differential Abundance signature

Single-cell technologies have revolutionized our understanding of cellular dynamics by allowing researchers to investigate individual cell responses under various conditions, such as comparing diseased versus healthy states. Many differential abundance methods have been developed in this field, however, the understanding of the gene signatures obtained from those methods is often incomplete, requiring the integration of cell type information and other biological factors to yield interpretable and meaningful results. To better interpret the gene signatures generated in the differential abundance analysis, we developed iDAS to classify the gene signatures into multiple categories.  

## Installation

```{r}
devtools::install_github("SydneyBioX/iDAS")
```


## Example

Run a two-way analysis to classify genes based on two relevant factors.

```{r}
res=iDAS_2F(Z= X,
 f1=pcelltype,f2=pcell_stats,random=NULL,test_func="lm",
 Sig_cutoff=0.02, Sig = 0.1,Int= 0.01, F1 = 0.01,F2 = 0.01,
 adj_method="BH",f1name=NULL,f2name=NULL,randomname=NULL)

```

Run a three-way analysis to classify genes based on three relevant factors.


```{r}
res=iDAS_3F(Z = X,
f1 = all.timepoint, f2 = all.pcellstats, f3 = all.pcelltype, random = all.pid,
test_func = "lmer",
Sig_cutoff = 0.02,Int = 0.01,F1 = 0.01,F2=0.02,F3=0.01,adj_method = "BH")
```


The iDAS package is still under development to meet Bioconductor standards. If you have any questions, please don't hesitate to open an issue.
