# Run Model Test and Return p-value and Statistic

Fits two models using either `lm` or `lmer` (if a random effect is
provided) and performs an ANOVA to compare them. Returns the p-value and
a test statistic.

## Usage

``` r
anova_test(data, formula0, formula1, model_fit_function)
```

## Arguments

- data:

  A data frame containing the data for the models.

- formula0:

  A character string representing the formula for the null model.

- formula1:

  A character string representing the formula for the alternative model.

- model_fit_function:

  A character string specifying the model fitting function to use ("lm"
  or "lmer").

## Value

A numeric vector of length 2, where the first element is the p-value and
the second element is the test statistic.

## Examples

``` r
if (FALSE) { # \dontrun{
# Using linear models:
run_test(data = my_data, "Y ~ x", "Y ~ x + z", model_fit_function = "lm")
} # }
```
