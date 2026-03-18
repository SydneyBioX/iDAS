# Run Either an ANOVA or Permutation-Based ANOVA Test

This function calls either
[`anova_test`](https://sydneybiox.github.io/iDAS/reference/anova_test.md)
or
[`perm_anova_test`](https://sydneybiox.github.io/iDAS/reference/perm_anova_test.md)
based on the specified `test_function`. If `test_function` is
`"parametric"`, a standard ANOVA is performed. If `test_function` is
`"permutation"`, a permutation-based ANOVA is performed.

## Usage

``` r
run_test(
  data,
  formula0,
  formula1,
  model_fit_function,
  test_function,
  n_perm = 1000,
  BPPARAM = MulticoreParam()
)
```

## Arguments

- data:

  A data frame containing the variables referenced in `formula0` and
  `formula1`.

- formula0:

  A formula specifying the reduced model.

- formula1:

  A formula specifying the full model.

- model_fit_function:

  A function used to fit the model (e.g. `lm`, `glm`, etc.).

- test_function:

  A character string specifying which test to run. Must be either
  `"parametric"` or `"permutation"`.

- n_perm:

  An integer specifying the number of permutations to use when
  `test_function = "permutation"`. Default is `1000`.

- BPPARAM:

  An optional
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  object for parallel computation. Defaults to
  [`MulticoreParam()`](https://rdrr.io/pkg/BiocParallel/man/MulticoreParam-class.html).

## Value

The result object returned by either
[`anova_test`](https://sydneybiox.github.io/iDAS/reference/anova_test.md)
or
[`perm_anova_test`](https://sydneybiox.github.io/iDAS/reference/perm_anova_test.md),
depending on `test_function`.

## See also

[`anova_test`](https://sydneybiox.github.io/iDAS/reference/anova_test.md),
[`perm_anova_test`](https://sydneybiox.github.io/iDAS/reference/perm_anova_test.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  # Example using a linear model and standard ANOVA:
  run_test(
    data = mtcars,
    formula0 = mpg ~ 1,
    formula1 = mpg ~ cyl,
    model_fit_function = "lm",
    test_function = "parametric"
  )

  # Example using a linear model and permutation-based ANOVA:
  run_test(
    data = mtcars,
    formula0 = mpg ~ 1,
    formula1 = mpg ~ cyl,
    model_fit_function = "lm",
    test_function = "permutation",
    n_perm = 500
  )
} # }
```
