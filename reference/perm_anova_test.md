# Permutation Test for ANOVA Model Comparison

This function performs a permutation test for comparing two nested
models using an ANOVA F-test. It fits two models (a null model and an
alternative model) to the original data, extracts the observed F
statistic, and then generates a null distribution of F statistics by
permuting the response variable. Parallel execution is supported via the
BiocParallel framework.

## Usage

``` r
perm_anova_test(
  data,
  formula0,
  formula1,
  model_fit_function,
  n_perm = 50,
  BPPARAM = MulticoreParam()
)
```

## Arguments

- data:

  A data frame containing the variables used in the models. The response
  variable must be named `Y`.

- formula0:

  A formula or character string specifying the null (reduced) model.

- formula1:

  A formula or character string specifying the alternative (full) model.

- model_fit_function:

  A character string naming the function used to fit the models (e.g.,
  `"lm"`, `"glm"`, `"lmer"`). The function is retrieved via
  [`get`](https://rdrr.io/r/base/get.html) and must be available in the
  current environment.

- n_perm:

  An integer specifying the number of permutations to perform. Default
  is `50`.

- BPPARAM:

  A
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  object specifying the parallel backend used for permutations. Defaults
  to
  [`MulticoreParam()`](https://rdrr.io/pkg/BiocParallel/man/MulticoreParam-class.html).

## Value

A numeric vector of length two:

- The first element is the permutation p-value (`perm_p_value`).

- The second is the observed F statistic (`observed_F`).

## Details

In each permutation iteration, the response variable `Y` in `data` is
randomly shuffled, while the predictors remain unchanged. This generates
a null distribution for the F statistic under the hypothesis of no
association between the response and the predictors.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example using lm for two nested models
data <- data.frame(
  Y = rnorm(100),
  factor1 = gl(2, 50),
  factor2 = gl(2, 25, length = 100)
)

formula0 <- "Y ~ factor1 + factor2"
formula1 <- "Y ~ factor1 * factor2"

result <- perm_anova_test(
  data = data,
  formula0 = formula0,
  formula1 = formula1,
  model_fit_function = "lm",
  n_perm = 100
)

# The first element is the permutation p-value
print(result[1])
# The second element is the observed F statistic
print(result[2])
} # }
```
