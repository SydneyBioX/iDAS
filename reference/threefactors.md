# iDAS: Interpretable Differential Analysis of Genes with Three Factors

This function implements the iDAS (Interpretable Differential Analysis
Signature) framework to identify features associated with three
experimental factors (`factor1`, `factor2`, and `factor3`), as well as
their interactions. The analysis involves an overall model test,
interaction tests (two-way and three-way), and main effects tests.
Results include adjusted p-values and test statistics for each feature.

## Usage

``` r
threefactors(
  Z,
  factor1,
  factor2,
  factor3,
  random_effect = NULL,
  model_fit_function = "lm",
  test_function = "parametric",
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
  p_adjust_method = "BH",
  factor1_name = NULL,
  factor2_name = NULL,
  factor3_name = NULL,
  random_effect_name = NULL,
  ...
)
```

## Arguments

- Z:

  A numeric matrix or data frame where each column represents a feature
  and each row represents an observation.

- factor1:

  A factor or vector representing the primary experimental factor.

- factor2:

  A factor or vector representing the secondary experimental factor.

- factor3:

  A factor or vector representing the tertiary experimental factor.

- random_effect:

  An optional factor or vector for random effects (e.g., subject ID).
  Use `NULL` if not applicable.

- model_fit_function:

  A character string specifying the model-fitting function (e.g., `"lm"`
  or `"lmer"`). Defaults to `"lm"`.

- test_function:

  A character string specifying the testing function to use (e.g.,
  `"parametric"` or `"permutation"`). Defaults to `"parametric"`.

- pval_quantile_cutoff:

  A numeric threshold for the quantile-based filtering of overall
  p-values. Defaults to `0.02`.

- pval_cutoff_full:

  A numeric p-value cutoff for the overall model test. Defaults to
  `0.05`.

- pval_cutoff_interaction:

  A numeric p-value cutoff for the omnibus interaction test. Defaults to
  `0.01`.

- pval_cutoff_factor1:

  A numeric p-value cutoff for testing the main effect of `factor1`.
  Defaults to `0.01`.

- pval_cutoff_factor2:

  A numeric p-value cutoff for testing the main effect of `factor2`.
  Defaults to `0.01`.

- pval_cutoff_factor3:

  A numeric p-value cutoff for testing the main effect of `factor3`.
  Defaults to `0.01`.

- pval_cutoff_int12:

  Numeric p-value cutoff for the two-way interaction between factor1 and
  factor2. Defaults to `0.01`.

- pval_cutoff_int13:

  Numeric p-value cutoff for the two-way interaction between factor1 and
  factor3. Defaults to `0.01`.

- pval_cutoff_int23:

  Numeric p-value cutoff for the two-way interaction between factor2 and
  factor3. Defaults to `0.01`.

- pval_cutoff_int123:

  Numeric p-value cutoff for the three-way interaction among all
  factors. Defaults to `0.01`.

- p_adjust_method:

  A character string specifying the method used to adjust p-values
  (e.g., `"BH"`). Defaults to `"BH"`.

- factor1_name:

  Optional label for `factor1`.

- factor2_name:

  Optional label for `factor2`.

- factor3_name:

  Optional label for `factor3`.

- random_effect_name:

  Optional label for the random effect.

- ...:

  Additional arguments passed to internal functions, model-fitting
  routines, or test functions.

## Value

A list containing:

- `pval_matrix`:

  A matrix of adjusted p-values for each gene, including main effects
  and interactions.

- `stat_matrix`:

  A matrix of corresponding test statistics.

- `class_df`:

  A data frame classifying each gene based on the significance of main
  effects and interactions.

## Details

Internally, the function:

1.  Builds the appropriate model formulas for each gene, depending on
    `model_fit_function` (e.g., `lm` vs. `lmer`) and whether
    `random_effect` is provided.

2.  Performs an overall significance test for each gene (the
    `pval_cutoff_full` threshold).

3.  For those genes passing the overall test, conducts an omnibus
    interaction test and further specific tests (main effects or
    two-way/three-way interactions) controlled by the respective p-value
    cutoffs.

Multiple testing corrections are applied based on `p_adjust_method`.

## Examples

``` r
if (FALSE) { # \dontrun{
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
  test_function = "parametric",
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
} # }
```
