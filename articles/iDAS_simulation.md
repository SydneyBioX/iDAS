# Simulate Pseudobulk and apply iDAS to find factor associated gene list

## Summary

We start with simulate two-factor RNA-seq Pseudobulk and apply iDAS
functions on the simulated data.

Our models are:

Model1:
\\y\_{ijm}=\mu+\alpha_i+\beta_j+(\alpha\beta)\_{ij}+\epsilon\_{ijm}\\

Model2: \\y\_{ijm}=\mu+\alpha_i+\beta_j+\epsilon\_{ijm}\\

Model3: \\y\_{ijm}=\mu+\alpha_i+\epsilon\_{ijm}\\

Model4: \\y\_{ijm}=\mu+\beta_j+\epsilon\_{ijm}\\

Model5: \\y\_{ijm}=\mu+\epsilon\_{ijm}\\

Where \\\alpha\\ and \\\beta\\ are two factors, and
\\(\alpha\beta)\_{ij}\\ represents the interaction effect,which I will
denote as \\\gamma\\ in my simulation code.

Thus, my task is to simulate different groups of genes using Models 1,
2, 3, 4, and 5.

And also modify the \\\epsilon\_{ijm}\\ term to match the three
conditions.

### Simulation strategy

- Data that meets the ANOVA assumptions:

  1.  First, we generated prior parameters for \\\mu_0\\, \\\alpha_0\\,
      \\\beta_0\\, \\(\alpha\beta)\_0\\
  2.  We then generated parameters of \\\mu_i\\,, \\\alpha_i\\,
      \\\beta_i\\,\\(\alpha\beta)\_i\\ follow normal distribution
      \\N(\mu_0,1)\\,\\N(\alpha_0,1)\\,\\N(\beta_0,1)\\,\\N((\alpha\beta)\_0,1)\\
  3.  Error terms was generated from \\N(0,\Sigma)\\, \\\Sigma\\ is the
      covariance matrix, the diagonal element is 1 and the off-diagonal
      element is ρ, represents the correlation between the genes.
  4.  Generate samples from \\Model1,…Model5\\ based on the steps above.

- Data with heteroscedastic residuals:

  1.  Error terms were generated from a normal distribution, with
      \\\Sigma\\ based on group variance. For each group, variance
      values were simulated from a uniform distribution.

- Data with non-normal residuals:

  1.  Error terms were generated from a gamma distribution, and all
      values were shifted by subtracting the mean to ensure that the
      residuals included both positive and negative values.

#### Functions for simulating Pseudobulk and meta information

``` r
generate_sample_metadata <- function(n_samples,
                                     cell_type_props,
                                     treatment_props) {
  # set.seed(42)  # For reproducibility

  n_cell_types <- length(cell_type_props)
  n_treatment_levels <- length(treatment_props)

  if (sum(cell_type_props) != 1) stop("Cell type proportions must sum to 1.")
  if (sum(treatment_props) != 1) stop("Treatment proportions must sum to 1.")

  # Compute number of samples per cell type
  cell_type_sizes <- round(n_samples * cell_type_props)
  cell_type_labels <- rep(1:n_cell_types, times = cell_type_sizes)

  # Adjust rounding errors
  diff_samples <- n_samples - sum(cell_type_sizes)
  if (diff_samples != 0) {
    cell_type_sizes[which.max(cell_type_sizes)] <-
      cell_type_sizes[which.max(cell_type_sizes)] + diff_samples
  }

  # Assign treatment effects based on user-defined proportions
  treatment_labels <- numeric(n_samples)
  sample_index <- 1

  for (cell in 1:n_cell_types) {
    n_cell_samples <- cell_type_sizes[cell]
    treatment_sizes <- round(n_cell_samples * treatment_props)

    # Adjust rounding errors
    diff_treatment <- n_cell_samples - sum(treatment_sizes)
    if (diff_treatment != 0) {
      treatment_sizes[which.max(treatment_sizes)] <-
        treatment_sizes[which.max(treatment_sizes)] + diff_treatment
    }

    # Assign treatment labels to this cell type
    treatment_dist <- rep(1:length(treatment_props), times = treatment_sizes)
    treatment_labels[sample_index:(sample_index + n_cell_samples - 1)] <-
      sample(treatment_dist)
    sample_index <- sample_index + n_cell_samples
  }

  return(list(
    cell_type_labels = cell_type_labels,
    treatment_labels = treatment_labels
  ))
}



simulate_pseudobulk_rnaseq <- function(n_samples, n_genes,
                                       cell_type_labels, treatment_labels,
                                       model_type = 1, SD = 3) {
  # set.seed(42)  # Ensure reproducibility

  n_cell_types <- length(unique(cell_type_labels))
  n_treatment_levels <- length(unique(treatment_labels))

  # **Simulating priors using rnorm()**
  # mu_prior <- rnorm(1, mean = 20, sd = SD)
  # Treatment effects
  alpha_prior <- c(0, rnorm(n_treatment_levels - 1, mean = 0, sd = SD^2))
  # Cell type effects
  beta_prior <- stats::rnorm(n_cell_types, mean = 0, sd = SD^2)
  # Interaction terms
  gamma_prior <- stats::rnorm(n_treatment_levels * n_cell_types, mean = 0, sd = SD^2)

  # **Simulating each effect using mvrnorm() for correlation structure**
  # mu_effect <-  mvrnorm(1, mu = rep(mu_prior, n_genes), Sigma = diag(n_genes))
  treatment_effect <- MASS::mvrnorm(1,
    mu = alpha_prior[treatment_labels],
    Sigma = diag(n_samples)
  )
  cell_type_effect <- MASS::mvrnorm(1,
    mu = beta_prior[cell_type_labels],
    Sigma = diag(n_samples)
  )
  interaction_effect <- MASS::mvrnorm(1,
    mu = gamma_prior[(cell_type_labels - 1) *
      n_treatment_levels +
      treatment_labels],
    Sigma = diag(n_samples)
  )

  # **Modify Expression Matrix Based on Model Type**
  if (model_type == 1) { # Full Model (alpha + beta + gamma)
    expression_matrix <- # mu_effect +
      matrix(treatment_effect, nrow = n_genes, ncol = n_samples, byrow = TRUE) +
      matrix(cell_type_effect, nrow = n_genes, ncol = n_samples, byrow = TRUE) +
      matrix(interaction_effect, nrow = n_genes, ncol = n_samples, byrow = TRUE) #+
    # error_terms
  } else if (model_type == 2) { # Additive Model (alpha + beta)
    expression_matrix <- # mu_effect +
      matrix(treatment_effect, nrow = n_genes, ncol = n_samples, byrow = TRUE) +
      matrix(cell_type_effect, nrow = n_genes, ncol = n_samples, byrow = TRUE) #+
    # error_terms
  } else if (model_type == 3) { # Only Treatment Effect (alpha)
    expression_matrix <- # mu_effect +
      matrix(treatment_effect, nrow = n_genes, ncol = n_samples, byrow = TRUE) #+
    # error_terms
  } else if (model_type == 4) { # Only Cell Type Effect (beta)
    expression_matrix <- # mu_effect +
      matrix(cell_type_effect, nrow = n_genes, ncol = n_samples, byrow = TRUE) #+
    # error_terms
  } else if (model_type == 5) { # Null Model (No Effects)
    expression_matrix <- # mu_effect+
      matrix(0, nrow = n_genes, ncol = n_samples) # + error_terms
  } else {
    stop("Invalid model selection. Choose 1, 2, 3, 4, or 5.")
  }

  return(expression_matrix)
}
```

#### Functions for simulating residuals

We not only simulate the ideal situation but also simulate situations
where the assumptions of parametric-based tests may not hold, such as
heteroscedastic residuals and non-normal residuals.

``` r
generate_heteroscedastic_residuals <- function(n_samples, n_genes,
                                               cell_type_labels,
                                               treatment_labels,
                                               min_sd = 1, max_sd = 5) {
  # set.seed(42)  # Ensure reproducibility

  # Create unique group labels (T1_C1, T2_C2, etc.)
  group_labels <- paste0("T", treatment_labels, "_C", cell_type_labels)
  unique_groups <- unique(group_labels)

  # Assign different standard deviations per group
  group_variances <- stats::runif(length(unique_groups),
    min = min_sd, max = max_sd
  )
  names(group_variances) <- unique_groups

  # Generate heteroscedastic residuals
  error_terms <- matrix(0, nrow = n_genes, ncol = n_samples)

  for (i in 1:n_samples) {
    group <- group_labels[i]
    sigma <- group_variances[group] # Get group-specific standard deviation
    # Sample residuals
    error_terms[, i] <- stats::rnorm(n_genes, mean = 0, sd = sigma)
  }

  return(error_terms)
}

generate_nonnormal_residuals <- function(n_samples,
                                         n_genes,
                                         cell_type_labels,
                                         treatment_labels,
                                         dist_type = "gamma") {
  # set.seed(42)  # Ensure reproducibility

  # Create unique group labels (T1_C1, T2_C2, etc.)
  group_labels <- paste0("T", treatment_labels, "_C", cell_type_labels)
  unique_groups <- unique(group_labels)

  # Assign different standard deviations per group
  group_variances <- stats::runif(length(unique_groups), min = 0, max = 1.5)
  names(group_variances) <- unique_groups

  # Generate non-normal residuals
  error_terms <- matrix(0, nrow = n_genes, ncol = n_samples)

  for (i in 1:n_samples) {
    group <- group_labels[i]
    sigma <- group_variances[group] # Get group-specific standard deviation

    # Generate residuals based on the specified distribution
    if (dist_type == "gamma") {
      temp <- stats::rgamma(n_genes, shape = 0.5, scale = sigma) #-2 * sigma
      error_terms[, i] <- temp - mean(temp) # Right-skewed
    } else if (dist_type == "lognormal") {
      # Long-tailed
      error_terms[, i] <- stats::rlnorm(n_genes,
        meanlog = log(sigma),
        sdlog = 0.5
      ) - sigma
    } else if (dist_type == "laplace") {
      # Heavy-tailed
      error_terms[, i] <- stats::rnorm(n_genes,
        mean = 0,
        sd = sigma
      ) * sign(runif(n_genes) - 0.5)
    } else {
      stop("Invalid distribution type. Choose 'gamma',
           'lognormal', or 'laplace'.")
    }
  }

  return(error_terms)
}
```

#### Simulate data

We want to simulate 60 patient samples. There are two cell types within
these samples (cell type 1 and cell type 2) and two treatment
phenotypes. The pseudobulk data is generated at the cell type level,
where each sample represents a cell type with its treatment phenotype
(either responding or non-responding to the therapy).

``` r
n_samples <- 60 # Total number of samples
n_genes_list <- c(50, 50, 50, 50, 50) # Number of genes per model type
names(n_genes_list) <- c("Interaction", "Additive", "Treatment", "Celltype", "Null")
cell_type_props <- c(0.5, 0.5) # 50% cell type 1, 50% cell type 2
treatment_props <- c(0.5, 0.5) # 50% Treatment 1, 50% Treatment 2

metadata <- generate_sample_metadata(n_samples, cell_type_props, treatment_props)

table(metadata$cell_type_labels, metadata$treatment_labels)
#>    
#>      1  2
#>   1 15 15
#>   2 15 15
```

##### 1.Simulate expression data for each model (Models 1 to 5)

``` r
expr_list <- lapply(1:5, function(i) {
  simulate_pseudobulk_rnaseq(n_samples,
    n_genes_list[i],
    metadata$cell_type_labels,
    metadata$treatment_labels,
    model_type = i
  )
})

# Combine the simulated expression matrices by row
expr_combined <- do.call(rbind, expr_list)
total_genes <- sum(n_genes_list)
```

##### 2.Generate error/residual matrices

``` r
rho <- 0
Sigma <- (1 - rho) * diag(total_genes) + rho * matrix(1, nrow = total_genes, ncol = total_genes)

error_terms <- t(MASS::mvrnorm(n_samples, mu = rep(0, total_genes), Sigma = Sigma))
error_terms_het <- generate_heteroscedastic_residuals(n_samples, total_genes,
  metadata$cell_type_labels,
  metadata$treatment_labels,
  min_sd = 1, max_sd = sqrt(10)
)
error_terms_nonnormal <- generate_nonnormal_residuals(n_samples, total_genes,
  metadata$cell_type_labels,
  metadata$treatment_labels,
  dist_type = "gamma"
)
```

##### 3.Create three expression matrices by adding different error terms

`expr_combined_null` represents the ideal simulation, where the data has
two factor effects and normal residuals.

`expr_combined_het` represents the situation where the assumption of
homoscedastic residuals does not hold.

`expr_combined_nonnormal` represents the situation where the assumption
of normal residuals does not hold.

``` r
expr_combined_null <- expr_combined + error_terms
expr_combined_het <- expr_combined + error_terms_het
expr_combined_nonnormal <- expr_combined + error_terms_nonnormal
```

##### 4.Generate gene names

True labels are embedded as the first part, e.g., “Model1”.

``` r
gene_names_combined <- unlist(lapply(1:5, function(i) {
  offset <- if (i == 1) 0 else sum(n_genes_list[1:(i - 1)])
  paste0("Model", i, "_Gene", seq_len(n_genes_list[i]) + offset)
}))
```

##### 5.Generate sample names from metadata

``` r
sample_names <- paste0("Sample", 1:n_samples, "_T", 
                       metadata$treatment_labels, "_C", 
                       metadata$cell_type_labels)

# Assign row and column names to all three matrices
rownames(expr_combined_null) <- gene_names_combined
rownames(expr_combined_het) <- gene_names_combined
rownames(expr_combined_nonnormal) <- gene_names_combined

colnames(expr_combined_null) <- sample_names
colnames(expr_combined_het) <- sample_names
colnames(expr_combined_nonnormal) <- sample_names
```

##### 6.Create sample annotation (same for all simulations)

``` r
sample_annotation <- data.frame(
  Treatment = factor(metadata$treatment_labels,
    labels = paste0("T", 1:length(unique(metadata$treatment_labels)))
  ),
  Cell_Type = factor(metadata$cell_type_labels,
    labels = paste0("C", 1:length(unique(metadata$cell_type_labels)))
  )
)
rownames(sample_annotation) <- sample_names
```

### Run iDAS analysis

Please note that permutation tests may take a very long time to produce
results.

``` r
oc_null <- iDAS::iDAS(
  Z = data.frame(t(expr_combined_null), check.names = FALSE),
  factor1 = sample_annotation$Cell_Type,
  factor2 = sample_annotation$Treatment,
  model_fit_function = "lm", pval_quantile_cutoff = NULL,
  pval_cutoff_full = 0.05, pval_cutoff_interaction = 0.05,
  pval_cutoff_factor1 = 0.05, pval_cutoff_factor2 = 0.05,
  p_adjust_method = "BH"
)
#> [1] "Running two-factor model"
#> [1] "full model"
#> [1] "interaction,f1,f2"

# oc_null_permu=iDAS::iDAS(Z = data.frame(t(expr_combined_null), check.names = FALSE),
#                      factor1= sample_annotation$Cell_Type,
#                      factor2= sample_annotation$Treatment,
#                      model_fit_function = "lm",
#                      pval_cutoff_full = 0.05, pval_cutoff_interaction = 0.05,
#                      pval_cutoff_factor1 = 0.05, pval_cutoff_factor2 = 0.05,
#                      p_adjust_method = "BH",test_function = "perm_anova_test",
#                      n_perm = 100,BPPARAM = SerialParam())


oc_het <- iDAS::iDAS(
  Z = data.frame(t(expr_combined_het), check.names = FALSE),
  factor1 = sample_annotation$Cell_Type,
  factor2 = sample_annotation$Treatment,
  model_fit_function = "lm",
  pval_cutoff_full = 0.05, pval_cutoff_interaction = 0.05,
  pval_cutoff_factor1 = 0.05, pval_cutoff_factor2 = 0.05,
  p_adjust_method = "BH"
)
#> [1] "Running two-factor model"
#> [1] "full model"
#> [1] "interaction,f1,f2"

# oc_het_permu=iDAS::iDAS(Z = data.frame(t(expr_combined_het),check.names = FALSE),
#                      factor1= sample_annotation$Cell_Type,
#                      factor2= sample_annotation$Treatment,
#                      model_fit_function = "lm",
#                      pval_cutoff_full = 0.05, pval_cutoff_interaction = 0.05,
#                      pval_cutoff_factor1 = 0.05, pval_cutoff_factor2 = 0.05,
#                      p_adjust_method = "BH",test_function = "perm_anova_test",
#                      n_perm = 100,BPPARAM = MulticoreParam(10))

oc_nonnormal <- iDAS::iDAS(
  Z = data.frame(t(expr_combined_nonnormal), check.names = FALSE),
  factor1 = sample_annotation$Cell_Type,
  factor2 = sample_annotation$Treatment,
  model_fit_function = "lm",
  pval_cutoff_full = 0.05, pval_cutoff_interaction = 0.05,
  pval_cutoff_factor1 = 0.05, pval_cutoff_factor2 = 0.05,
  p_adjust_method = "BH"
)
#> [1] "Running two-factor model"
#> [1] "full model"
#> [1] "interaction,f1,f2"

# oc_nonnormal_permu=iDAS::iDAS(Z = data.frame(t(expr_combined_nonnormal),check.names = FALSE),
#                      factor1= sample_annotation$Cell_Type,
#                      factor2= sample_annotation$Treatment,
#                      model_fit_function = "lm",
#                      pval_cutoff_full = 0.05, pval_cutoff_interaction = 0.05,
#                      pval_cutoff_factor1 = 0.05, pval_cutoff_factor2 = 0.05,
#                      p_adjust_method = "BH",test_function = "perm_anova_test",
#                      n_perm = 100,BPPARAM = MulticoreParam(10))
```

#### Compare Parametric-based and Permuation-based test results

``` r
# Extract cluster results and replace NAs in Sig1 with "non-sig"
clusterres_null <- data.frame(oc_null$class_df)
clusterres_null$Sig1[is.na(clusterres_null$Sig1)] <- "non-sig"

clusterres_het <- data.frame(oc_het$class_df)
clusterres_het$Sig1[is.na(clusterres_het$Sig1)] <- "non-sig"

clusterres_nonnormal <- data.frame(oc_nonnormal$class_df)
clusterres_nonnormal$Sig1[is.na(clusterres_nonnormal$Sig1)] <- "non-sig"

# clusterres_null_permu <- data.frame(oc_null_permu$class_df)
# clusterres_null_permu$Sig1[is.na(clusterres_null_permu$Sig1)] <- "non-sig"
  
# clusterres_het_permu <- data.frame(oc_het_permu$class_df)
# clusterres_het_permu$Sig1[is.na(clusterres_het_permu$Sig1)] <- "non-sig"
#
# clusterres_nonnormal_permu <- data.frame(oc_nonnormal_permu$class_df)
# clusterres_nonnormal_permu$Sig1[is.na(clusterres_nonnormal_permu$Sig1)] <- "non-sig"


# Extract the true labels from gene names (e.g., "Model1" from "Model1_GeneX")
true_labels <- do.call(rbind, strsplit(clusterres_null$varname, "_"))[, 1]

# Compute ARI between true labels and predicted labels
ari_true_null <- aricode::ARI(true_labels, clusterres_null$Sig1)
ari_true_het <- aricode::ARI(true_labels, clusterres_het$Sig1)
ari_true_nonnormal <- aricode::ARI(true_labels, clusterres_nonnormal$Sig1)



# ari_true_null_permu      <- aricode::ARI(true_labels, clusterres_null_permu$Sig1)
# ari_true_het_permu       <- ARI(true_labels, clusterres_het_permu$Sig1)
# ari_true_nonnormal_permu <- ARI(true_labels, clusterres_nonnormal_permu$Sig1)


print(paste0("ARI of parametric-based test: ", ari_true_null))
#> [1] "ARI of parametric-based test: 0.750682420499178"
# print(paste0("concordance result from permutation test: ",ari_true_null_permu))


print("ARI of parametric-based test with heterscedastic residuals: ")
#> [1] "ARI of parametric-based test with heterscedastic residuals: "
print(ari_true_het)
#> [1] 0.001059704
print("ARI of parametric-based test with non-normal residuals does not hold.: ")
#> [1] "ARI of parametric-based test with non-normal residuals does not hold.: "
print(ari_true_nonnormal)
#> [1] 0.001839015
```

**Session Information**

``` r
sessionInfo()
#> R version 4.5.3 (2026-03-11)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] BiocStyle_2.38.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] Matrix_1.7-4        jsonlite_2.0.0      compiler_4.5.3     
#>  [4] BiocManager_1.30.27 Rcpp_1.1.1          iDAS_0.1.1         
#>  [7] parallel_4.5.3      jquerylib_0.1.4     splines_4.5.3      
#> [10] systemfonts_1.3.2   textshaping_1.0.5   boot_1.3-32        
#> [13] BiocParallel_1.44.0 yaml_2.3.12         fastmap_1.2.0      
#> [16] aricode_1.0.3       lattice_0.22-9      R6_2.6.1           
#> [19] knitr_1.51          rbibutils_2.4.1     MASS_7.3-65        
#> [22] nloptr_2.2.1        bookdown_0.46       desc_1.4.3         
#> [25] minqa_1.2.8         bslib_0.10.0        rlang_1.1.7        
#> [28] cachem_1.1.0        xfun_0.56           fs_1.6.7           
#> [31] sass_0.4.10         cli_3.6.5           pkgdown_2.2.0      
#> [34] Rdpack_2.6.6        digest_0.6.39       grid_4.5.3         
#> [37] lme4_2.0-1          lifecycle_1.0.5     nlme_3.1-168       
#> [40] reformulas_0.4.4    evaluate_1.0.5      codetools_0.2-20   
#> [43] ragg_1.5.1          rmarkdown_2.30      tools_4.5.3        
#> [46] htmltools_0.5.9
```
