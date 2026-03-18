# iDAS: Interpretable Differential Abundance Analysis

This function implements the iDAS (Interpretable Differential Abundance
Analysis Signature) framework for analyzing differential abundance gene
signatures. It serves as a comprehensive wrapper supporting both
two-factor and three-factor experimental designs. When a third factor is
provided, a three-factor analysis is performed via the `threefactors`
function; otherwise, a two-factor analysis is executed via the
`twofactors` function.

## Usage

``` r
iDAS(
  Z,
  factor1,
  factor2,
  factor3 = NULL,
  random_effect = NULL,
  model_fit_function = "lm",
  p_adjust_method_for_factors_and_interation = FALSE,
  pval_quantile_cutoff = 0.02,
  pval_cutoff_full = 0.05,
  pval_cutoff_interaction = 0.01,
  pval_cutoff_factor1 = 0.01,
  pval_cutoff_factor2 = 0.01,
  pval_cutoff_factor3 = NULL,
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
  (e.g., microbial taxa, metabolites) to be analyzed.

- factor1:

  A vector or factor representing the first experimental factor.

- factor2:

  A vector or factor representing the second experimental factor.

- factor3:

  An optional vector or factor representing the third experimental
  factor. If provided, a three-factor analysis is performed. Default is
  `NULL`.

- random_effect:

  An optional vector or factor representing a random effect (e.g.,
  subject ID). Default is `NULL`.

- model_fit_function:

  A character string indicating the model fitting function to use (e.g.,
  `"lm"` for linear models or `"lmer"` for mixed-effects models).
  Default is `"lm"`.

- p_adjust_method_for_factors_and_interation:

  Logical or character, specifying whether p-values for factors and
  interactions should be adjusted. Default is `FALSE`.

- pval_quantile_cutoff:

  Numeric value representing the quantile cutoff for overall
  significance testing. Default is `0.02`.

- pval_cutoff_full:

  Numeric p-value cutoff for the overall (full) model test. Default is
  `0.05`.

- pval_cutoff_interaction:

  Numeric p-value cutoff for the interaction test. Default is `0.01`.

- pval_cutoff_factor1:

  Numeric p-value cutoff for testing the main effect of the first
  factor. Default is `0.01`.

- pval_cutoff_factor2:

  Numeric p-value cutoff for testing the main effect of the second
  factor. Default is `0.01`.

- pval_cutoff_factor3:

  Numeric p-value cutoff for testing the main effect of the third factor
  (if `factor3` is provided). Default is `NULL`.

- p_adjust_method:

  Character string specifying the method for p-value adjustment (e.g.,
  `"BH"`). Default is `"BH"`.

- factor1_name:

  Optional character string for naming the first factor. Default is
  `NULL`.

- factor2_name:

  Optional character string for naming the second factor. Default is
  `NULL`.

- random_effect_name:

  Optional character string for naming the random effect. Default is
  `NULL`.

- ...:

  Additional arguments passed to the `twofactors` or `threefactors`
  functions.

## Value

A list containing the results from the differential abundance analysis.
The output typically includes matrices of p-values and test statistics,
and a data frame classifying features based on significance.

## Details

The function distinguishes between a two-factor design and a
three-factor design based on whether `factor3` is provided. It prints
`"Running three-factor model"` when executing a three-factor analysis
and `"Running two-factor model"` for a two-factor analysis (for
debugging purposes).

## Examples

``` r
if (FALSE) { # \dontrun{
  # Example using two factors
  result_two <- iDAS(Z = my_feature_matrix, factor1 = group1, factor2 = group2)

  # Example using three factors
  result_three <- iDAS(Z = my_feature_matrix, factor1 = group1, factor2 = group2,
  factor3 = timepoint)
} # }
```
