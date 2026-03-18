# Interpretable differential abundance analysis (three-way analysis)

This function will be no longer be used.

## Usage

``` r
iDAS_3F(
  Z,
  f1,
  f2,
  f3,
  random = NULL,
  test_func = "lm",
  Sig_cutoff = 0.02,
  Sig = 0.05,
  Int = 0.01,
  F1 = 0.01,
  F2 = 0.01,
  F3 = 0.01,
  F1F2 = 0.01,
  F1F3 = 0.01,
  F2F3 = 0.01,
  threeways = 0.01,
  adj_method = "BH",
  f1name = NULL,
  f2name = NULL,
  f3name = NULL,
  randomname = NULL
)
```

## Arguments

- Z:

  A matrix/dataframe of omics or gene expression data, row as sample.

- f1:

  A vector of factor 1 variables.

- f2:

  A vector of factor 2 variables.

- f3:

  A vector of factor 3 variables.

- random:

  A vector of random effect term of ANOVA analysis, by default is NULL,
  which means the model doesn't include random effect term.

- test_func:

  Testing function used, either stats::lm or lme4::lmer. By default is
  "lm".

- Sig_cutoff:

  No effect test significance level is defined by a fraction value to
  indicate when ordering the p-values and defining the top X% as the
  significance level. If both Sig_cutoff and Sig are set, the algorithm
  will, by default, use Sig_cutoff to determine the significance level.

- Sig:

  No effect test significance level is defined directly by a fraction
  value, by default is 0.05

- Int:

  Interaction effect test significance level,by default is 0.01

- F1:

  Factor 1 effect test significance level, by default is 0.01

- F2:

  Factor 2 effect test significance level, by default is 0.01

- F3:

  Factor 3 effect test significance level, by default is 0.01

- F1F2:

  F1 and F2 interaction effect test significance level, by default is
  0.01

- F1F3:

  F1 and F3 interaction effect test significance level, by default is
  0.01

- F2F3:

  F2 and F3 interaction effect test significance level, by default is
  0.01

- threeways:

  three-way effect test significance level, by default is 0.01

- adj_method:

  Pvalue adjust method. See p.adjust. By default is "BH".

- f1name:

  The column name of factor 1, by default is F1.

- f2name:

  The column name of factor 2, by default is F2.

- f3name:

  The column name of factor 3, by default is F3.

- randomname:

  The column name of random effect term, by default is Random.

## Value

A list of hypothesis test outcome, P_mat is the pvalue matrix of all
tests, S_mat is the statistics matrix of all test, cls_df is the
Classification data frame of all tests.

## Examples

``` r
#res=iDAS_3F(Z = X,
#f1 = all.timepoint, f2 = all.pcellstats, f3 = all.pcelltype, random = all.pid, test_func = "lmer",
#Sig_cutoff = 0.02,Int = 0.01,F1 = 0.01,F2=0.02,F3=0.01,adj_method = "BH")
```
