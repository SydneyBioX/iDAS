# Check Factor Name

Processes input factors and a random effect by converting them into data
frames with specified column names. It also collapses the column names
into single strings if multiple columns are present. Supports an
optional third factor.

## Usage

``` r
check_factor_name(
  factor1_name,
  factor2_name,
  factor3_name = NULL,
  random_effect_name,
  factor1,
  factor2,
  factor3 = NULL,
  random_effect = NULL
)
```

## Arguments

- factor1_name:

  A character string representing the name for factor1. Defaults to
  "factor1" if NULL.

- factor2_name:

  A character string representing the name for factor2. Defaults to
  "factor2" if NULL.

- factor3_name:

  A character string representing the name for factor3. Defaults to
  "factor3" if `factor3` is provided and factor3_name is NULL.

- random_effect_name:

  A character string representing the name for the random effect.
  Defaults to "random_effect" if NULL.

- factor1:

  A factor or vector representing the first factor.

- factor2:

  A factor or vector representing the second factor.

- factor3:

  A factor or vector representing the third factor. Optional; defaults
  to `NULL`.

- random_effect:

  A factor or vector representing the random effect. May be `NULL`.

## Value

A list containing:

- factor1_tmp:

  Collapsed name for factor1.

- factor2_tmp:

  Collapsed name for factor2.

- factor3_tmp:

  Collapsed name for factor3 (or `NULL` if `factor3` is not provided).

- random_effect_tmp:

  Collapsed name for the random effect (or `NULL` if `random_effect` is
  `NULL`).

- factor1:

  Data frame for factor1.

- factor2:

  Data frame for factor2.

- factor3:

  Data frame for factor3 (or `NULL`).

- random_effect:

  Data frame for random_effect (or `NULL`).
