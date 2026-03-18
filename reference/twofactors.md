# iDAS: Interpretable Differential Analysis of Genes with Two Factors

This function implements the iDAS (Interpretable Differential Analysis
Signature) framework to identify features associated with two
experimental factors (`factor1` and `factor2`), as well as their
interaction. The analysis involves an overall model test, followed by
specific tests for interactions and main effects. Results include
adjusted p-values and test statistics for each feature.

## Usage

``` r
twofactors(
  Z,
  factor1,
  factor2,
  random_effect = NULL,
  model_fit_function = "lm",
  p_adjust_method_for_factors_and_interation = FALSE,
  pval_quantile_cutoff = 0.02,
  pval_cutoff_full = 0.05,
  test_function = "parametric",
  pval_cutoff_interaction = 0.01,
  pval_cutoff_factor1 = 0.01,
  pval_cutoff_factor2 = 0.01,
  p_adjust_method = "BH",
  factor1_name = NULL,
  factor2_name = NULL,
  random_effect_name = NULL,
  ...
)
```

## Arguments

- Z:

  A numeric matrix or data frame where each column represents a feature
  (e.g., gene expression) and each row represents an observation.

- factor1:

  A factor or vector representing the first categorical variable.

- factor2:

  A factor or vector representing the second categorical variable.

- random_effect:

  An optional factor or vector representing a random effect (e.g.,
  subject ID). Use `NULL` if not applicable.

- model_fit_function:

  A character string specifying the model fitting function to use.
  Acceptable values are `"lm"` for linear models or `"lmer"` for
  mixed-effects models.

- p_adjust_method_for_factors_and_interation:

  Logical indicating whether to apply p-value adjustment for follow-up
  tests on factors and interaction. Defaults to `FALSE`.

- pval_quantile_cutoff:

  Numeric; a quantile cutoff used for an overall significance filter in
  the analysis. Defaults to `0.02`.

- pval_cutoff_full:

  Numeric; the p-value threshold for the overall model test. Defaults to
  `0.05`.

- test_function:

  A character string specifying which test function to use for
  statistical comparison. Defaults to `"parametric"`. Other valid
  functions (if implemented) might include `"permutation"`, etc.

- pval_cutoff_interaction:

  Numeric; the p-value threshold for testing the interaction effect.
  Defaults to `0.01`.

- pval_cutoff_factor1:

  Numeric; the p-value threshold for testing the main effect of
  `factor1`. Defaults to `0.01`.

- pval_cutoff_factor2:

  Numeric; the p-value threshold for testing the main effect of
  `factor2`. Defaults to `0.01`.

- p_adjust_method:

  A character string specifying the method for p-value adjustment (e.g.,
  `"BH"`). Defaults to `"BH"`.

- factor1_name:

  An optional character string to label `factor1` in model formulas; if
  `NULL`, a default is used.

- factor2_name:

  An optional character string to label `factor2` in model formulas; if
  `NULL`, a default is used.

- random_effect_name:

  An optional character string to label the random effect in model
  formulas; if `NULL`, a default is used.

- ...:

  Additional arguments passed to internal functions (e.g., model-fitting
  functions or test functions).

## Value

A list containing:

- `pval_matrix`:

  A matrix of adjusted p-values for the overall test, interaction test,
  and main effects.

- `stat_matrix`:

  A matrix of test statistics corresponding to the computed p-values.

- `class_df`:

  A data frame summarizing the significance classification for each
  feature based on the analysis.

## Examples

``` r
if (FALSE) { # \dontrun{
results <- twofactors(
  Z = my_data_matrix,
  factor1 = my_factor1,
  factor2 = my_factor2,
  model_fit_function = "lm",
  p_adjust_method_for_factors_and_interation = FALSE,
  pval_quantile_cutoff = 0.02,
  pval_cutoff_full = 0.05,
  test_function = "parametric",
  pval_cutoff_interaction = 0.01,
  pval_cutoff_factor1 = 0.01,
  pval_cutoff_factor2 = 0.01,
  p_adjust_method = "BH",
  factor1_name = "Group",
  factor2_name = "Treatment"
)
# Inspect the output
print(results$pval_matrix)
print(results$stat_matrix)
print(results$class_df)
} # }
```
