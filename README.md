# iDAS

Interpretable Differential Abundance Signature

Single-cell technologies have revolutionized our understanding of cellular dynamics by allowing researchers to investigate individual cell responses under various conditions, such as comparing diseased versus healthy states. Many differential abundance methods have been developed in this field, however, the understanding of the gene signatures obtained from those methods is often incomplete, requiring the integration of cell type information and other biological factors to yield interpretable and meaningful results. To better interpret the gene signatures generated in the differential abundance analysis, we developed iDAS to classify the gene signatures into multiple categories.  

## Installation

```{r}
devtools::install_github("SydneyBioX/iDAS")
```



## Example


```{r}
  # Example using two factors
  result_twoway <- iDAS(Z = my_feature_matrix, factor1 = group1, factor2 = group2)

  # Example using three factors
  result_threeway <- iDAS(Z = my_feature_matrix, factor1 = group1, factor2 = group2, factor3 = timepoint)
```




```{r}
## Not run: 
# Generate sample data
set.seed(123)
Z <- matrix(rnorm(1000), ncol = 10)
colnames(Z)=paste0("gene",1:10)
factor1 <- as.factor(rep(1:2, each = 5))
factor2 <- as.factor(rep(1:2, times = 5))
factor3 <- as.factor(rep(1:2, length.out = 10))

# Run the differential analysis using iDAS
result <- threefactors(
  Z, factor1, factor2, factor3,
  model_fit_function = "lm",
  test_function = "anova_test",
  pval_quantile_cutoff = 0.02,
  pval_cutoff_full = 0.05,
  pval_cutoff_interaction = 0.01,
  pval_cutoff_factor1 = 0.01,
  pval_cutoff_factor2 = 0.01,
  pval_cutoff_factor3 = 0.01,
  pval_cutoff_int12 = 0.01,
  pval_cutoff_int13 = 0.01,
  pval_cutoff_int23 = 0.01,
  pval_cutoff_int123 = 0.01,
  p_adjust_method = "BH"
)

# Inspect results
head(result$pval_matrix)
head(result$stat_matrix)
head(result$class_df)

## End(Not run)
```

The iDAS package is still under development to meet Bioconductor standards. If you have any questions, please don't hesitate to open an issue.
