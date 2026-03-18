# Process Input

Converts a factor or vector to a data frame and assigns a specified
column name.

## Usage

``` r
process_input(x, name)
```

## Arguments

- x:

  A factor or vector that should be converted.

- name:

  A character value to assign as the column name.

## Value

A data frame with one column named `name` containing `x`, or the
original `x` if it is not a factor or vector.
