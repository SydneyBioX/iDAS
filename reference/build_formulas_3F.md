# Build Model Formulas for Three-Factor Differential Analysis

This function constructs various model formulas for differential
analysis involving three factors. It generates formulas for the full
model with all interactions, alternative two-way interaction models,
models excluding specific interactions, a null interaction model, and
individual main effects models. The formulas vary depending on whether a
random effect is included and the specified model fitting function.

## Usage

``` r
build_formulas_3F(formatted_factors, test_func, random_effect)
```

## Arguments

- formatted_factors:

  A list containing formatted factor names for the analysis. Expected
  elements include `factor1_tmp`, `factor2_tmp`, `factor3_tmp`, and
  optionally `random_effect_tmp`.

- test_func:

  A character string specifying the model fitting function to use (e.g.,
  `"lm"`).

- random_effect:

  An optional argument representing the random effect. If `NULL`,
  formulas for fixed effects are generated.

## Value

A list of character strings containing the following formulas:

- lm_full:

  The full model formula including all three-way interactions.

- lm_int_alt2way:

  An alternative model formula including only two-way interactions.

- lm_int_altf1f2:

  A model formula excluding the interaction between factor1 and factor2.

- lm_int_altf1f3:

  A model formula excluding the interaction between factor1 and factor3.

- lm_int_altf2f3:

  A model formula excluding the interaction between factor2 and factor3.

- lm_int_null:

  A null model formula with only main effects (no interactions).

- lm_f1:

  A model formula for the main effect of factor1 only.

- lm_f2:

  A model formula for the main effect of factor2 only.

- lm_f3:

  A model formula for the main effect of factor3 only.

## Details

Depending on whether a random effect is provided, the function creates
formulas for either a fixed effects model or a mixed-effects model
(including random effects using `(1|random_effect)`).

## Examples

``` r
if (FALSE) { # \dontrun{
  # Example for a fixed effects model:
  formatted_factors <- list(factor1_tmp = "Group", factor2_tmp = "Treatment",
  factor3_tmp = "Time")
  formulas <- build_formulas_3F(formatted_factors, test_func = "lm", random_effect = NULL)
  print(formulas$lm_full)

  # Example for a random effects model:
  formatted_factors$random_effect_tmp <- "Subject"
  formulas_re <- build_formulas_3F(formatted_factors, test_func = "lm",
  random_effect = formatted_factors$random_effect_tmp)
  print(formulas_re$lm_full)
} # }
```
