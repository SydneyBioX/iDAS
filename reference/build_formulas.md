# Construct Model Formulas for Two- or Three-Factor Analysis

This function builds model formula strings for either two or three
factors, with optional random effects. When `test_func` is `"lm"`, the
function assumes no random effect (`random_effect = NULL`). When
`test_func` is `"lmer"`, a random effect must be provided. Supported
factors:

- Two-factor scenario: `factor1_tmp` and `factor2_tmp`.

- Three-factor scenario: `factor1_tmp`, `factor2_tmp`, and
  `factor3_tmp`.

## Usage

``` r
build_formulas(formatted_factors, test_func, random_effect)
```

## Arguments

- formatted_factors:

  A named list containing the formatted factor names (e.g.,
  `factor1_tmp`, `factor2_tmp`, `factor3_tmp`) and optionally
  `random_effect_tmp` if using a mixed-effects model. Typically produced
  by an internal or auxiliary formatting function.

- test_func:

  A character string indicating which model-fitting approach to use.
  Currently supports `"lm"` (no random effect) and `"lmer"`
  (mixed-effects).

- random_effect:

  A factor or `NULL`. If `test_func = "lmer"`, this should be a valid
  factor; if `test_func = "lm"`, it must be `NULL`.

## Value

A named list of formula strings. The exact list structure depends on the
number of non-`NULL` factors in `formatted_factors`:

- **Two-factor case**:
  `list( full = <full_model_formula>, int_null = <interaction_null_formula>, int_alt = <interaction_alt_formula>, f1 = <factor1_formula>, f2 = <factor2_formula> )`

      \item **Three-factor case**:
        \code{list(
              lm_full = <full_model_formula>,
              lm_int_alt2way = <two_way_interaction_formula>,
              lm_int_altf1f2 = <omitting_f1f2_interaction_formula>,
              lm_int_altf1f3 = <omitting_f1f3_interaction_formula>,
              lm_int_altf2f3 = <omitting_f2f3_interaction_formula>,
              lm_int_null = <no_interaction_formula>,
              lm_f1 = <factor1_formula>,
              lm_f2 = <factor2_formula>,
              lm_f3 = <factor3_formula>
            )}

## Details

The function calculates `num_factors` internally by counting how many
non-`NULL` items exist in `formatted_factors` (apart from the random
effect) and dividing by 2. This must yield either 2 or 3 (representing
two- or three-factor designs). Otherwise, it raises an error.

1.  For `test_func = "lm"`, `random_effect` must be `NULL`.

2.  For `test_func = "lmer"`, `random_effect` must be a valid factor.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example for a two-factor linear model:
my_factors <- list(
  factor1_tmp = "Group",
  factor2_tmp = "Treatment",
  random_effect_tmp = NULL  # Not used for 'lm'
)
formulas_2f <- build_formulas(my_factors, test_func = "lm", random_effect = NULL)
print(formulas_2f)

# Example for a three-factor mixed-effects model:
my_factors_3 <- list(
  factor1_tmp = "Group",
  factor2_tmp = "Treatment",
  factor3_tmp = "Time",
  random_effect_tmp = "Subject"  # random effect for 'lmer'
)
formulas_3f <- build_formulas(my_factors_3, test_func = "lmer", random_effect = factor("Subject"))
print(formulas_3f)
} # }
```
