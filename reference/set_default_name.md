# Set Default Name

Returns the default name if the provided name is NULL.

## Usage

``` r
set_default_name(name, default)
```

## Arguments

- name:

  A character value representing the current name.

- default:

  A character value representing the default name to use if `name` is
  NULL.

## Value

The original `name` if not NULL, otherwise `default`.
