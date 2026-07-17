#' Interpretable differential abundance analysis (two-way analysis)
#'
#' This function will be no longer be used.
#' @param Z A matrix/dataframe of omics or gene expression data, row as sample.
#' @param factor1 A vector of the first factor variable.
#' @param factor2 A vector of the second factor variable.
#' @param random_effect A vector of the random effect term of ANOVA analysis,
#' by default is NULL, which means the model doesn't include a random effect term.
#' @param model_fit_function Model fitting function used, either stats::lm or lme4::lmer. By default is "lm".
#' @param pval_quantile_cutoff A fraction used to determine the significance threshold
#' for the overall (full) model p-values.
#' @param pval_cutoff_full The p-value threshold for the overall (full) model, by default is 0.05.
#' @param pval_cutoff_interaction The p-value threshold for the interaction effect, by default is 0.01.
#' @param pval_cutoff_factor1 The p-value threshold for the main effect of factor1, by default is 0.01.
#' @param pval_cutoff_factor2 The p-value threshold for the main effect of factor2, by default is 0.01.
#' @param p_adjust_method P-value adjustment method. See p.adjust. By default is "BH".
#' @param factor1_name The column name of the first factor variable, by default is "factor1".
#' @param factor2_name The column name of the second factor variable, by default is "factor2".
#' @param random_effect_name The column name of the random effect term, by default is "random_effect".
#'
#' @return A list of hypothesis test outcomes. pval_matrix is the matrix of p-values for all tests,
#' stat_matrix is the matrix of test statistics, and class_df is the data frame of class results.
# Legacy implementation retained for development reference only.
#' @importFrom stats anova formula p.adjust
#' @importFrom lme4 lmer
#' @keywords internal
#' @examples
#' # res = iDAS_2F(Z = X,
#' #               factor1 = pcelltype, factor2 = pcell_stats, random_effect = NULL,
#' #               model_fit_function = "lm",
#' #               pval_quantile_cutoff = 0.02, pval_cutoff_full = 0.05,
#' #               pval_cutoff_interaction = 0.01, pval_cutoff_factor1 = 0.01,
#' #               pval_cutoff_factor2 = 0.01,
#' #               p_adjust_method = "BH", factor1_name = NULL, factor2_name = NULL,
#' #               random_effect_name = NULL)
iDAS_2F = function(Z, factor1, factor2, random_effect = NULL,
                   model_fit_function = "lm",
                   pval_quantile_cutoff = 0.02, pval_cutoff_full = 0.05,
                   pval_cutoff_interaction = 0.01, pval_cutoff_factor1 = 0.01, pval_cutoff_factor2 = 0.01,
                   p_adjust_method = "BH", factor1_name = NULL, factor2_name = NULL, random_effect_name = NULL) {

  pval_matrix = stat_matrix <- matrix(NA, nrow = ncol(Z), ncol = 4)

  # Format factor variables and names
  formatted_factors = check_factor_name2(factor1_name, factor2_name, random_effect_name, factor1, factor2, random_effect)

  # Build model formulas based on the fitting function and random effects
  if(model_fit_function == "lm" & is.null(random_effect)) {
    formula_full <- paste("Y ~ (", formatted_factors$factor1_tmp, ")*(", formatted_factors$factor2_tmp, ")", sep = "")
    formula_interaction_alt <- paste("Y ~ (", formatted_factors$factor1_tmp, "):(", formatted_factors$factor2_tmp, ")", sep = "")
    formula_interaction_null <- paste("Y ~ (", formatted_factors$factor1_tmp, ")+(", formatted_factors$factor2_tmp, ")", sep = "")
    formula_factor1 <- paste("Y ~ (", formatted_factors$factor1_tmp, ")", sep = "")
    formula_factor2 <- paste("Y ~ (", formatted_factors$factor2_tmp, ")", sep = "")
  } else if(model_fit_function == "lmer" & (!is.null(random_effect))) {
    formula_full <- paste("Y ~ (", formatted_factors$factor1_tmp, ")*(", formatted_factors$factor2_tmp,
                          ") + (1|", formatted_factors$random_effect_tmp, ")", sep = "")
    formula_interaction_alt <- paste("Y ~ (", formatted_factors$factor1_tmp, "):(", formatted_factors$factor2_tmp,
                                     ") + (1|", formatted_factors$random_effect_tmp, ")", sep = "")
    formula_interaction_null <- paste("Y ~ (", formatted_factors$factor1_tmp, ")+(", formatted_factors$factor2_tmp,
                                      ") + (1|", formatted_factors$random_effect_tmp, ")", sep = "")
    formula_factor1 <- paste("Y ~ (", formatted_factors$factor1_tmp, ") + (1|", formatted_factors$random_effect_tmp, ")", sep = "")
    formula_factor2 <- paste("Y ~ (", formatted_factors$factor2_tmp, ") + (1|", formatted_factors$random_effect_tmp, ")", sep = "")
  } else {
    stop("Invalid model specification: For 'lm', random_effect must be NULL; for 'lmer', random_effect must be provided.")
  }

  lm_formula = c("full" = formula_full, "int_null" = formula_interaction_null,
                 "int_alt" = formula_interaction_alt,
                 "f1" = formula_factor1, "f2" = formula_factor2)

  pval_full = stat_full = c()  # p-values and stats from the full (overall) test

  if(model_fit_function == "lm" & is.null(random_effect)) {
    null_model_formula <- "Y ~ 1"
    fit_model_fn = get("lm")
    alt_model_formula <- lm_formula[1]
  } else {
    null_model_formula <- "Y ~ 1 + (1|random_effect)"
    fit_model_fn = get("lmer")
    alt_model_formula <- lm_formula[1]
  }

  # Loop over features in Z
  for (i in 1:ncol(Z)) {
    if(model_fit_function == "lm" & is.null(random_effect)) {
      dat <- data.frame(Y = Z[, i], formatted_factors$factor1, formatted_factors$factor2)
      M0 <- fit_model_fn(formula(null_model_formula), data = dat)
      M1 <- fit_model_fn(formula(alt_model_formula), data = dat)
      anova_res <- anova(M0, M1, test = "F")
      p <- anova_res[2, ncol(anova_res)]
      S <- anova_res[2, ncol(anova_res) - 1]
    } else {
      dat <- data.frame(Y = Z[, i], formatted_factors$factor1, formatted_factors$factor2,
                        random_effect = formatted_factors$random_effect)
      M0 <- fit_model_fn(formula(null_model_formula), data = dat)
      M1 <- fit_model_fn(formula(alt_model_formula), data = dat)
      anova_res <- anova(M0, M1, test = "F", refit = FALSE)
      p <- anova_res[2, ncol(anova_res)]
      S <- anova_res[2, ncol(anova_res) - 2]
    }
    pval_full[i] <- p
    stat_full[i] <- S
  }

  if (is.null(p_adjust_method)) {
    p_sig_adj <- pval_full
  } else {
    p_sig_adj <- p.adjust(pval_full, method = p_adjust_method)
  }

  pval_matrix[, 1] <- p_sig_adj
  stat_matrix[, 1] <- stat_full

  overall_class <- rep("sig", ncol(Z))
  if (!is.null(pval_quantile_cutoff)) {
    cutoff = sort(p_sig_adj)[length(p_sig_adj) * pval_quantile_cutoff]
    overall_class[p_sig_adj > cutoff] <- "non-sig"
    class_df <- data.frame(Sig0 = overall_class)
    sig_feature_indices <- which(overall_class == "sig")
    if (length(sig_feature_indices) > 1) {
      Z_sig <- Z[, sig_feature_indices]
    } else {
      Z_sig <- data.frame(Z[, sig_feature_indices])
      colnames(Z_sig) <- colnames(Z)[sig_feature_indices]
      rownames(Z_sig) <- rownames(Z)
    }
  } else if(pval_cutoff_full) {
    overall_class[p_sig_adj > pval_cutoff_full] <- "non-sig"
    class_df <- data.frame(Sig0 = overall_class)
    sig_feature_indices <- which(overall_class == "sig")
    if (length(sig_feature_indices) > 1) {
      Z_sig <- Z[, sig_feature_indices]
    } else {
      Z_sig <- data.frame(Z[, sig_feature_indices])
      colnames(Z_sig) <- colnames(Z)[sig_feature_indices]
      rownames(Z_sig) <- rownames(Z)
    }
  }

  # Initialize temporary variables for the interaction and main effect tests
  pval_interaction = stat_interaction <- rep(NA, ncol(Z_sig))
  pval_factor1 = stat_factor1 <- rep(NA, ncol(Z_sig))
  pval_factor2 = stat_factor2 <- rep(NA, ncol(Z_sig))

  for (i in 1:ncol(Z_sig)) {
    # Interaction test: compare formula_interaction_null vs. formula_full
    if(model_fit_function == "lm" & is.null(random_effect)) {
      calc0 <- lm_formula["int_null"]
      calc1 <- lm_formula["full"]
      dat <- data.frame(Y = Z_sig[, i], formatted_factors$factor1, formatted_factors$factor2)
      M0 <- fit_model_fn(formula(calc0), data = dat)
      M1 <- fit_model_fn(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F")
      p <- anova_res[2, ncol(anova_res)]
      S <- anova_res[2, ncol(anova_res) - 1]
    } else {
      calc0 <- lm_formula["int_null"]
      calc1 <- lm_formula["full"]
      dat <- data.frame(Y = Z_sig[, i], formatted_factors$factor1, formatted_factors$factor2,
                        random_effect = formatted_factors$random_effect)
      M0 <- fit_model_fn(formula(calc0), data = dat)
      M1 <- fit_model_fn(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F", refit = FALSE)
      p <- anova_res[2, ncol(anova_res)]
      S <- anova_res[2, ncol(anova_res) - 2]
    }
    pval_interaction[i] <- p
    stat_interaction[i] <- S

    # Main effect test for factor1: use formula_factor2 as null vs. formula_interaction_null
    if(model_fit_function == "lm" & is.null(random_effect)) {
      calc0 <- lm_formula["f2"]
      calc1 <- lm_formula["int_null"]
      dat <- data.frame(Y = Z_sig[, i], formatted_factors$factor1, formatted_factors$factor2)
      M0 <- fit_model_fn(formula(calc0), data = dat)
      M1 <- fit_model_fn(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F")
      p <- anova_res[2, ncol(anova_res)]
      S <- anova_res[2, ncol(anova_res) - 1]
    } else {
      calc0 <- lm_formula["f2"]
      calc1 <- lm_formula["int_null"]
      dat <- data.frame(Y = Z_sig[, i], formatted_factors$factor1, formatted_factors$factor2,
                        random_effect = formatted_factors$random_effect)
      M0 <- fit_model_fn(formula(calc0), data = dat)
      M1 <- fit_model_fn(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F", refit = FALSE)
      p <- anova_res[2, ncol(anova_res)]
      S <- anova_res[2, ncol(anova_res) - 2]
    }
    pval_factor1[i] <- p
    stat_factor1[i] <- S

    # Main effect test for factor2: use formula_factor1 as null vs. formula_interaction_null
    if(model_fit_function == "lm" & is.null(random_effect)) {
      calc0 <- lm_formula["f1"]
      calc1 <- lm_formula["int_null"]
      dat <- data.frame(Y = Z_sig[, i], formatted_factors$factor1, formatted_factors$factor2)
      M0 <- fit_model_fn(formula(calc0), data = dat)
      M1 <- fit_model_fn(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F")
      p <- anova_res[2, ncol(anova_res)]
      S <- anova_res[2, ncol(anova_res) - 1]
    } else {
      calc0 <- lm_formula["f1"]
      calc1 <- lm_formula["int_null"]
      dat <- data.frame(Y = Z_sig[, i], formatted_factors$factor1, formatted_factors$factor2,
                        random_effect = formatted_factors$random_effect)
      M0 <- fit_model_fn(formula(calc0), data = dat)
      M1 <- fit_model_fn(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F", refit = FALSE)
      p <- anova_res[2, ncol(anova_res)]
      S <- anova_res[2, ncol(anova_res) - 2]
    }
    pval_factor2[i] <- p
    stat_factor2[i] <- S
  }

  # Adjust p-values for the interaction and main effects
  if (is.null(p_adjust_method)) {
    pval_matrix[sig_feature_indices, 2] = pval_interaction
    pval_matrix[sig_feature_indices, 3] = pval_factor1
    pval_matrix[sig_feature_indices, 4] = pval_factor2
    stat_matrix[sig_feature_indices, 2] = stat_interaction
    stat_matrix[sig_feature_indices, 3] = stat_factor1
    stat_matrix[sig_feature_indices, 4] = stat_factor2
  } else {
    pval_interaction = p.adjust(pval_interaction, method = p_adjust_method)
    pval_factor1 = p.adjust(pval_factor1, method = p_adjust_method)
    pval_factor2 = p.adjust(pval_factor2, method = p_adjust_method)
    pval_matrix[sig_feature_indices, 2] = pval_interaction
    pval_matrix[sig_feature_indices, 3] = pval_factor1
    pval_matrix[sig_feature_indices, 4] = pval_factor2
    stat_matrix[sig_feature_indices, 2] = stat_interaction
    stat_matrix[sig_feature_indices, 3] = stat_factor1
    stat_matrix[sig_feature_indices, 4] = stat_factor2
  }

  interaction_class <- rep("Additive", ncol(Z_sig))
  interaction_class[pval_interaction < pval_cutoff_interaction] = "Interaction"
  names(interaction_class) <- colnames(Z_sig)
  overall_class <- rep(NA, ncol(Z))
  idx_add <- which(interaction_class == "Additive")
  overall_class_subset <- interaction_class[idx_add]
  pval_factor2_subset <- pval_factor2[idx_add]
  pval_factor1_subset <- pval_factor1[idx_add]
  overall_class_subset[(pval_factor2_subset <= pval_cutoff_factor2) & (pval_factor1_subset > pval_cutoff_factor1)] = "Factor2"
  overall_class_subset[(pval_factor2_subset > pval_cutoff_factor2) & (pval_factor1_subset <= pval_cutoff_factor1)] = "Factor1"
  interaction_class[idx_add] <- overall_class_subset
  overall_class[sig_feature_indices] <- interaction_class
  class_df$Sig1 <- overall_class
  class_df$varname <- colnames(Z)
  colnames(pval_matrix) = colnames(stat_matrix) <- c("Full", "Interaction", "Factor1", "Factor2")
  rownames(pval_matrix) = rownames(stat_matrix) <- colnames(Z)
  rownames(class_df) <- colnames(Z)

  return(list(pval_matrix = pval_matrix, stat_matrix = stat_matrix, class_df = class_df))
}






#' Check the iDAS_2F input factors' name
#'
#' This function will be no longer be used.
#' @param factor1_name A string for the first factor variable's name.
#' @param factor2_name A string for the second factor variable's name.
#' @param random_effect_name A string for the random effect term's name.
#' @param factor1 A vector of the first factor variable.
#' @param factor2 A vector of the second factor variable.
#' @param random_effect A vector of the random effect term variables.
#'
#' @return A list of each factor's name and values.
#' @keywords internal
#' @examples
#' # formatted_factors = check_factor_name(factor1_name, factor2_name,
#' # random_effect_name, factor1, factor2, random_effect)
check_factor_name2 = function(factor1_name, factor2_name, random_effect_name, factor1, factor2, random_effect) {
  if (is.null(factor1_name)) {
    factor1_name = "factor1"
  }
  if (is.null(factor2_name)) {
    factor2_name = "factor2"
  }
  if (is.null(random_effect_name)) {
    random_effect_name = "random_effect"
  }

  if (is.factor(factor1)) {
    factor1 <- data.frame(factor1)
    colnames(factor1) <- factor1_name
  } else if (is.vector(factor1)) {
    factor1 <- data.frame(factor1)
    colnames(factor1) <- factor1_name
  }
  factor1_name <- colnames(factor1)

  if (is.factor(factor2)) {
    factor2 <- data.frame(factor2)
    colnames(factor2) <- factor2_name
  } else if (is.vector(factor2)) {
    factor2 <- data.frame(factor2)
    colnames(factor2) <- factor2_name
  }
  factor2_name <- colnames(factor2)

  if(is.null(random_effect)) {
    random_effect = NULL
  } else {
    if (is.factor(random_effect)) {
      random_effect <- data.frame(random_effect)
      colnames(random_effect) <- random_effect_name
    } else if (is.vector(random_effect)) {
      random_effect <- data.frame(random_effect)
      colnames(random_effect) <- random_effect_name
    }
  }
  random_effect_name = colnames(random_effect)

  if (length(factor1_name) > 1) {
    factor1_tmp <- paste(factor1_name, collapse = "+")
  } else {
    factor1_tmp <- factor1_name
  }
  if (length(factor2_name) > 1) {
    factor2_tmp <- paste(factor2_name, collapse = "+")
  } else {
    factor2_tmp <- factor2_name
  }
  if (!is.null(random_effect)) {
    if (length(random_effect_name) > 1) {
      random_effect_tmp <- paste(random_effect_name, collapse = "+")
    } else {
      random_effect_tmp <- random_effect_name
    }
  } else {
    random_effect_tmp <- NULL
  }

  return(list(factor1_tmp = factor1_tmp, factor2_tmp = factor2_tmp, random_effect_tmp = random_effect_tmp,
              factor1 = factor1, factor2 = factor2, random_effect = random_effect))
}
