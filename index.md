# iDAS [![](https://i.imgur.com/O25crxx.png "iDAS hex sticker")](https://github.com/SydneyBioX/iDAS)

**Interpretable Differential Abundance Signature**

iDAS is an R package for interpreting gene signatures obtained after a
differential abundance (DA) analysis. It uses a two-stage, ANOVA-based
framework to classify features according to the experimental factors and
interactions that explain their variation. This helps distinguish, for
example, signatures shared across cell states from cell-state-specific
treatment-response signatures.

iDAS supports:

- two- and three-factor experimental designs;
- fixed-effects models (`lm`) and mixed-effects models (`lmer`);
- parametric ANOVA and permutation-based tests; and
- input derived from any upstream DA method, including
  [miloR](https://github.com/MarioniLab/miloR),
  [DA-seq](https://github.com/KlugerLab/DAseq), and
  [Scissor](https://github.com/sunduanchen/Scissor).

The method was developed for pseudobulk single-cell data and can also be
used with other feature matrices and spatial summary measures when the
model assumptions are appropriate. See the [iDAS
paper](https://doi.org/10.1002/smtd.202500572) for the complete
framework, case studies, and discussion of limitations.

![Overview of the iDAS workflow](./reference/figures/iDAS_workflow.png)

Overview of the iDAS workflow

## Installation

iDAS is under active development toward Bioconductor standards. Install
the development version from GitHub:

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("SydneyBioX/iDAS")
```

Then load the package:

``` r

library(iDAS)
```

Please report bugs or request features through the [GitHub issue
tracker](https://github.com/SydneyBioX/iDAS/issues).

## Input format

The main [`iDAS()`](https://sydneybiox.github.io/iDAS/reference/iDAS.md)
function expects:

- `Z`: a numeric matrix or data frame with **observations in rows** and
  **features (for example, genes) in columns**;
- `factor1` and `factor2`: vectors with one value per row of `Z`;
- `factor3`: an optional third factor; and
- `random_effect`: an optional grouping vector, such as patient ID, for
  repeated-measures designs.

For single-cell RNA-seq, we recommend aggregating expression into
pseudobulk observations (for example, one row per patient and cell
state), then applying an appropriate transformation before iDAS. The
paper uses log2-transformed pseudobulk expression. Avoid treating
individual cells from the same biological sample as independent
replicates.

``` r

stopifnot(
  nrow(Z) == length(factor1),
  nrow(Z) == length(factor2),
  !is.null(colnames(Z))
)
```

## Quick start: two-factor analysis

This reproducible example simulates 40 pseudobulk observations and 20
genes. The first five genes have a treatment-by-cell-state interaction.

``` r

set.seed(123)

n <- 40
cell_state <- factor(rep(c("State_A", "State_B"), each = n / 2))
response <- factor(rep(rep(c("Responder", "Non_responder"), each = n / 4), 2))

Z <- matrix(rnorm(n * 20), nrow = n, ncol = 20)
colnames(Z) <- paste0("gene", seq_len(ncol(Z)))

# Add an interaction signal to genes 1-5.
interaction_group <- cell_state == "State_B" & response == "Responder"
Z[interaction_group, 1:5] <- Z[interaction_group, 1:5] + 3

result <- iDAS(
  Z = Z,
  factor1 = cell_state,
  factor2 = response,
  pval_quantile_cutoff = NULL, # use pval_cutoff_full instead
  pval_cutoff_full = 0.05
)

head(result$class_df)
head(result$pval_matrix)
head(result$stat_matrix)
```

By default, `pval_quantile_cutoff = 0.02` retains the top 2% of features
ranked by adjusted full-model p-value. Set it to `NULL`, as above, to
use the fixed `pval_cutoff_full` threshold instead.

## Repeated measures and three-factor designs

Use a mixed-effects model when observations are repeated within
subjects:

``` r

result_mixed <- iDAS(
  Z = Z,
  factor1 = cell_state,
  factor2 = response,
  random_effect = patient_id,
  model_fit_function = "lmer"
)
```

Supplying `factor3` automatically selects the three-factor analysis. For
example, a longitudinal study might model cell state, response, and
timepoint, with patient ID as a random effect:

``` r

result_three_factor <- iDAS(
  Z = Z_longitudinal,
  factor1 = cell_state,
  factor2 = response,
  factor3 = timepoint,
  random_effect = patient_id,
  model_fit_function = "lmer",
  pval_cutoff_factor3 = 0.01
)
```

Use `model_fit_function = "lm"` only when `random_effect = NULL`; use
`model_fit_function = "lmer"` when a random effect is supplied.

## Understanding the results

[`iDAS()`](https://sydneybiox.github.io/iDAS/reference/iDAS.md) returns
a list with three elements:

- `pval_matrix`: p-values (adjusted where requested) for the nested
  model comparisons;
- `stat_matrix`: the corresponding test statistics; and
- `class_df`: the feature-level signature classifications.

For a two-factor analysis, `class_df$Sig0` reports whether a feature
passed the full-model screen. For retained features, `class_df$Sig1` is
one of:

- `Interaction`: the association depends on the combination of factors;
- `Factor1`: primarily associated with factor 1;
- `Factor2`: primarily associated with factor 2; or
- `Additive`: associated with the additive model but not assigned
  exclusively to one factor.

The two-factor p-value and statistic matrices contain `Full`,
`Interaction`, `Factor1`, and `Factor2` columns. Three-factor results
additionally separate main effects (`F1`, `F2`, and `F3`), two-way
interactions (`F1F2`, `F1F3`, and `F2F3`), and three-way interactions.
Follow-up entries are `NA` when a feature did not reach the preceding
stage of the hierarchical testing procedure.

For full argument documentation, run:

``` r

?iDAS
vignette("iDAS", package = "iDAS")
```

## Statistical considerations

iDAS is intended to interpret candidate signatures, not to replace the
upstream DA analysis. Results can depend on the chosen DA method and on
the quality of the labels passed as factors. The parametric test assumes
that the linear-model assumptions are reasonable. A permutation test is
available for settings where these assumptions are questionable:

``` r

result_permutation <- iDAS(
  Z = Z,
  factor1 = cell_state,
  factor2 = response,
  test_function = "permutation",
  n_perm = 1000
)
```

Permutation testing is substantially slower, especially for large
feature sets. iDAS currently supports at most three factors; interaction
complexity and power become limiting as factors are added.

## Citation

If you use iDAS in your research, please cite:

> Yu, L., Lin, Y., Xu, X., Yang, P., and Yang, J. Y. H. (2026).
> Interpretable Differential Abundance Signature (iDAS). *Small
> Methods*, **10**(2), e2500572.
> <https://doi.org/10.1002/smtd.202500572>

BibTeX:

``` bibtex
@article{Yu2026iDAS,
  author  = {Yu, Lijia and Lin, Yingxin and Xu, Xiangnan and Yang, Pengyi and Yang, Jean Y. H.},
  title   = {Interpretable Differential Abundance Signature ({iDAS})},
  journal = {Small Methods},
  year    = {2026},
  volume  = {10},
  number  = {2},
  pages   = {e2500572},
  doi     = {10.1002/smtd.202500572}
}
```

## License

iDAS is released under the GPL license.
