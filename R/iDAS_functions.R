#' Run Model Test and Return p-value and Statistic
#'
#' Fits two models using either \code{lm} or \code{lmer} (if a random effect is provided)
#' and performs an ANOVA to compare them. Returns the p-value and a test statistic.
#'
#' @param data A data frame containing the data for the models.
#' @param formula0 A character string representing the formula for the null model.
#' @param formula1 A character string representing the formula for the alternative model.
#' @param model_fit_function A character string specifying the model fitting function to use ("lm" or "lmer").
#'
#' @return A numeric vector of length 2, where the first element is the p-value and the second element is the test statistic.
#'
#' @examples
#' \dontrun{
#' # Using linear models:
#' run_test(data = my_data, "Y ~ x", "Y ~ x + z", model_fit_function = "lm")
#' }
#' @keywords internal
anova_test <- function(data, formula0, formula1, model_fit_function) {
  func=get(model_fit_function)
  if (model_fit_function == "lm") {
    M0 <- func(formula(formula0), data = data)
    M1 <- func(formula(formula1), data = data)
    anova_res <- anova(M0, M1, test = "F")
    p_val <- anova_res[2, ncol(anova_res)]
    stat_val <- anova_res[2, ncol(anova_res) - 1]
  } else if (model_fit_function == "lmer"){
    M0 <- func(formula(formula0), data = data)
    M1 <- func(formula(formula1), data = data)
    anova_res <- anova(M0, M1, test = "F", refit = FALSE)
    p_val <- anova_res[2, ncol(anova_res)]
    stat_val <- anova_res[2, ncol(anova_res) - 2]
  }
  return(c(p_val, stat_val))
}


#' Permutation Test for ANOVA Model Comparison
#'
#' This function performs a permutation test for comparing two nested models using an ANOVA F-test.
#' It fits two models (a null model and an alternative model) to the original data, extracts the observed
#' F statistic, and then generates a null distribution of F statistics by permuting the response variable.
#' Parallel execution is supported via the \pkg{BiocParallel} framework.
#'
#' @param data A data frame containing the variables used in the models. The response variable must be
#'   named \code{Y}.
#' @param formula0 A formula or character string specifying the null (reduced) model.
#' @param formula1 A formula or character string specifying the alternative (full) model.
#' @param model_fit_function A character string naming the function used to fit the models
#'   (e.g., \code{"lm"}, \code{"glm"}, \code{"lmer"}). The function is retrieved via
#'   \code{\link[base]{get}} and must be available in the current environment.
#' @param n_perm An integer specifying the number of permutations to perform. Default is \code{50}.
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object specifying the parallel
#'   backend used for permutations. Defaults to \code{\link[BiocParallel]{MulticoreParam}()}.
#'
#' @return A numeric vector of length two:
#'   \itemize{
#'     \item The first element is the permutation p-value (\code{perm_p_value}).
#'     \item The second is the observed F statistic (\code{observed_F}).
#'   }
#'
#' @details In each permutation iteration, the response variable \code{Y} in \code{data} is randomly
#'   shuffled, while the predictors remain unchanged. This generates a null distribution for the
#'   F statistic under the hypothesis of no association between the response and the predictors.
#'
#' @import BiocParallel
#'
#' @examples
#' \dontrun{
#' # Example using lm for two nested models
#' data <- data.frame(
#'   Y = rnorm(100),
#'   factor1 = gl(2, 50),
#'   factor2 = gl(2, 25, length = 100)
#' )
#'
#' formula0 <- "Y ~ factor1 + factor2"
#' formula1 <- "Y ~ factor1 * factor2"
#'
#' result <- perm_anova_test(
#'   data = data,
#'   formula0 = formula0,
#'   formula1 = formula1,
#'   model_fit_function = "lm",
#'   n_perm = 100
#' )
#'
#' # The first element is the permutation p-value
#' print(result[1])
#' # The second element is the observed F statistic
#' print(result[2])
#' }
#'
#' @keywords internal
perm_anova_test <- function(data, formula0, formula1, model_fit_function, n_perm = 50, BPPARAM = MulticoreParam()) {
  # Fit the original models
  func=get(model_fit_function)
  M0 <- func(formula(formula0), data = data)
  M1 <- func(formula(formula1), data = data)
  anova_res <- anova(M0, M1, test = "F")
  observed_F <- anova_res[2, ncol(anova_res) - 1]  # Extract observed F statistic

  # Run permutation test in parallel
  perm_F <- bplapply(1:n_perm, function(i) {
    data_perm <- data
    data_perm$Y <- sample(data$Y)

    # Refit models using permuted data
    M0_perm <- func(formula(formula0), data = data_perm)
    M1_perm <- func(formula(formula1), data = data_perm)

    # Extract F statistic
    anova_perm <- anova(M0_perm, M1_perm, test = "F")
    return(anova_perm[2, ncol(anova_perm) - 1])
  }, BPPARAM = BPPARAM)

  # Convert list to numeric vector
  perm_F <- unlist(perm_F)

  # Compute permutation p-value
  perm_p_value <- mean(perm_F >= observed_F)

  return(c(perm_p_value, observed_F))
}





#' Run Either an ANOVA or Permutation-Based ANOVA Test
#'
#' This function calls either \code{\link{anova_test}} or \code{\link{perm_anova_test}}
#' based on the specified \code{test_function}. If \code{test_function} is
#' \code{"parametric"}, a standard ANOVA is performed. If \code{test_function} is
#' \code{"permutation"}, a permutation-based ANOVA is performed.
#'
#' @param data A data frame containing the variables referenced in \code{formula0} and \code{formula1}.
#' @param formula0 A formula specifying the reduced model.
#' @param formula1 A formula specifying the full model.
#' @param model_fit_function A function used to fit the model (e.g. \code{lm}, \code{glm}, etc.).
#' @param test_function A character string specifying which test to run. Must be either
#'   \code{"parametric"} or \code{"permutation"}.
#' @param n_perm An integer specifying the number of permutations to use when
#'   \code{test_function = "permutation"}. Default is \code{1000}.
#' @param BPPARAM An optional \code{\link[BiocParallel]{BiocParallelParam}} object for parallel computation.
#'   Defaults to \code{\link[BiocParallel]{MulticoreParam}()}.
#'
#' @return The result object returned by either \code{\link{anova_test}} or
#'   \code{\link{perm_anova_test}}, depending on \code{test_function}.
#'
#' @seealso \code{\link{anova_test}}, \code{\link{perm_anova_test}}
#'
#' @examples
#' \dontrun{
#'   # Example using a linear model and standard ANOVA:
#'   run_test(
#'     data = mtcars,
#'     formula0 = mpg ~ 1,
#'     formula1 = mpg ~ cyl,
#'     model_fit_function = "lm",
#'     test_function = "parametric"
#'   )
#'
#'   # Example using a linear model and permutation-based ANOVA:
#'   run_test(
#'     data = mtcars,
#'     formula0 = mpg ~ 1,
#'     formula1 = mpg ~ cyl,
#'     model_fit_function = "lm",
#'     test_function = "permutation",
#'     n_perm = 500
#'   )
#' }
#'
#' @keywords internal
run_test<- function(data, formula0, formula1, model_fit_function, test_function, n_perm = 1000,BPPARAM = MulticoreParam()) {
  if (test_function == "parametric") {
    return(anova_test(data, formula0, formula1, model_fit_function))
  } else if (test_function == "permutation") {
    return(perm_anova_test(data, formula0, formula1, model_fit_function, n_perm,BPPARAM = MulticoreParam()))
  } else {
    stop("Invalid test function name. Choose 'parametric' or 'permutation'.")
  }
}


#' iDAS: Interpretable Differential Analysis of Genes with Two Factors
#'
#' This function implements the iDAS (Interpretable Differential Analysis Signature) framework to identify
#' features associated with two experimental factors (\code{factor1} and \code{factor2}), as well as
#' their interaction. The analysis involves an overall model test, followed by specific tests for
#' interactions and main effects. Results include adjusted p-values and test statistics for each feature.
#'
#' @param Z A numeric matrix or data frame where each column represents a feature (e.g., gene expression) and each row represents an observation.
#' @param factor1 A factor or vector representing the first categorical variable.
#' @param factor2 A factor or vector representing the second categorical variable.
#' @param random_effect An optional factor or vector representing a random effect (e.g., subject ID).
#'   Use \code{NULL} if not applicable.
#' @param model_fit_function A character string specifying the model fitting function to use.
#'   Acceptable values are \code{"lm"} for linear models or \code{"lmer"} for mixed-effects models.
#' @param p_adjust_method_for_factors_and_interation Logical indicating whether to apply p-value adjustment
#'   for follow-up tests on factors and interaction. Defaults to \code{FALSE}.
#' @param pval_quantile_cutoff Numeric; a quantile cutoff used for an overall significance filter in the analysis.
#'   Defaults to \code{0.02}.
#' @param pval_cutoff_full Numeric; the p-value threshold for the overall model test. Defaults to \code{0.05}.
#' @param test_function A character string specifying which test function to use for statistical comparison.
#'   Defaults to \code{"parametric"}. Other valid functions (if implemented) might include \code{"permutation"}, etc.
#' @param pval_cutoff_interaction Numeric; the p-value threshold for testing the interaction effect. Defaults to \code{0.01}.
#' @param pval_cutoff_factor1 Numeric; the p-value threshold for testing the main effect of \code{factor1}.
#'   Defaults to \code{0.01}.
#' @param pval_cutoff_factor2 Numeric; the p-value threshold for testing the main effect of \code{factor2}.
#'   Defaults to \code{0.01}.
#' @param p_adjust_method A character string specifying the method for p-value adjustment (e.g., \code{"BH"}).
#'   Defaults to \code{"BH"}.
#' @param factor1_name An optional character string to label \code{factor1} in model formulas; if \code{NULL},
#'   a default is used.
#' @param factor2_name An optional character string to label \code{factor2} in model formulas; if \code{NULL},
#'   a default is used.
#' @param random_effect_name An optional character string to label the random effect in model formulas;
#'   if \code{NULL}, a default is used.
#' @param \dots Additional arguments passed to internal functions (e.g., model-fitting functions or test functions).
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{pval_matrix}}{A matrix of adjusted p-values for the overall test, interaction test,
#'   and main effects.}
#'   \item{\code{stat_matrix}}{A matrix of test statistics corresponding to the computed p-values.}
#'   \item{\code{class_df}}{A data frame summarizing the significance classification for each feature
#'   based on the analysis.}
#' }
#'
#' @examples
#' \dontrun{
#' results <- twofactors(
#'   Z = my_data_matrix,
#'   factor1 = my_factor1,
#'   factor2 = my_factor2,
#'   model_fit_function = "lm",
#'   p_adjust_method_for_factors_and_interation = FALSE,
#'   pval_quantile_cutoff = 0.02,
#'   pval_cutoff_full = 0.05,
#'   test_function = "parametric",
#'   pval_cutoff_interaction = 0.01,
#'   pval_cutoff_factor1 = 0.01,
#'   pval_cutoff_factor2 = 0.01,
#'   p_adjust_method = "BH",
#'   factor1_name = "Group",
#'   factor2_name = "Treatment"
#' )
#' # Inspect the output
#' print(results$pval_matrix)
#' print(results$stat_matrix)
#' print(results$class_df)
#' }
#'
#' @export
twofactors <- function(Z, factor1, factor2, random_effect = NULL,
                    model_fit_function = "lm",p_adjust_method_for_factors_and_interation=FALSE,
                    pval_quantile_cutoff = 0.02, pval_cutoff_full = 0.05,test_function="parametric",
                    pval_cutoff_interaction = 0.01, pval_cutoff_factor1 = 0.01, pval_cutoff_factor2 = 0.01,
                    p_adjust_method = "BH", factor1_name = NULL, factor2_name = NULL, random_effect_name = NULL,...) {

  # Initialize matrices for p-values and test statistics
  pval_matrix <- stat_matrix <- matrix(NA, nrow = ncol(Z), ncol = 4)

  # Format factor variables and names (assumes check_factor_name is defined)
  formatted_factors <- check_factor_name(factor1_name=factor1_name,factor2_name= factor2_name,
                                          random_effect_name=random_effect_name,factor1= factor1,factor2= factor2,
                                          random_effect= random_effect)

  lm_formula <- build_formulas(formatted_factors, model_fit_function, random_effect)
  # Set up the base model for overall tests
  if (model_fit_function == "lm" && is.null(random_effect)) {
    null_model_formula <- "Y ~ 1"
    fit_model_fn <- get("lm")
  } else {
    null_model_formula <- "Y ~ 1 + (1|random_effect)"
    fit_model_fn <- get("lmer")
  }


  num_features <- ncol(Z)
  pval_full <- stat_full <- numeric(num_features)

  print("full model")
  # Loop over features for the overall (full) model test
  for (i in 1:num_features) {
    dat <- if (is.null(random_effect)) {
      data.frame(Y = Z[, i], formatted_factors$factor1, formatted_factors$factor2)
    } else {
      data.frame(Y = Z[, i], formatted_factors$factor1, formatted_factors$factor2,
                 random_effect = formatted_factors$random_effect)
    }
    res <- run_test(dat,null_model_formula, lm_formula$full,  model_fit_function,test_function,...)
    pval_full[i] <- res[1]
    stat_full[i] <- res[2]
  }

  # Adjust overall p-values if needed and assign to matrix
  pval_full_adj <- if (is.null(p_adjust_method)) pval_full else p.adjust(pval_full, method = p_adjust_method)
  pval_matrix[, 1] <- pval_full_adj
  stat_matrix[, 1] <- stat_full

  # Classify overall significance
  overall_class <- rep("sig", num_features)
  if (!is.null(pval_quantile_cutoff)) {
    cutoff <- sort(pval_full_adj)[ceiling(length(pval_full_adj) * pval_quantile_cutoff)] ## ceiling
    overall_class[pval_full_adj > cutoff] <- "non-sig"
  } else {
    overall_class[pval_full_adj > pval_cutoff_full] <- "non-sig"
  }
  class_df <- data.frame(varname = colnames(Z),Sig0 = overall_class)
  sig_feature_indices <- which(overall_class == "sig")

  if (length(sig_feature_indices) == 0) {
    warning("No significant features found.")
    return(list(pval_matrix = pval_matrix, stat_matrix = stat_matrix, class_df = class_df))
  }

  # Subset significant features
  Z_sig <- if (length(sig_feature_indices) > 1) {
    Z[, sig_feature_indices]
  } else {
    data.frame(Z[, sig_feature_indices], check.names = FALSE)
  }

  # Initialize follow-up test vectors in one line
  pval_interaction = stat_interaction = pval_factor1 = stat_factor1 = pval_factor2 = stat_factor2 <- numeric(ncol(Z_sig))


  print("interaction,f1,f2")
  # Loop over significant features for follow-up tests
  for (i in 1:ncol(Z_sig)) {
    dat <- if (is.null(random_effect)) {
      data.frame(Y = Z_sig[, i], formatted_factors$factor1, formatted_factors$factor2)
    } else {
      data.frame(Y = Z_sig[, i], formatted_factors$factor1, formatted_factors$factor2,
                 random_effect = formatted_factors$random_effect)
    }
    # Interaction test: compare interaction-null vs. full model
    res_int <- run_test( dat,lm_formula$int_null, lm_formula$full, model_fit_function,test_function,...)
    pval_interaction[i] <- res_int[1]
    stat_interaction[i] <- res_int[2]
    # Main effect test for factor1: compare model with only factor2 vs. interaction-null
    res_f1 <- run_test( dat,lm_formula$f2, lm_formula$int_null,  model_fit_function,test_function,...)
    pval_factor1[i] <- res_f1[1]
    stat_factor1[i] <- res_f1[2]
    # Main effect test for factor2: compare model with only factor1 vs. interaction-null
    res_f2 <- run_test( dat,lm_formula$f1, lm_formula$int_null, model_fit_function,test_function,...)
    pval_factor2[i] <- res_f2[1]
    stat_factor2[i] <- res_f2[2]
  }

  # Adjust follow-up p-values if needed
  if (!is.null(p_adjust_method)) {
    if(p_adjust_method_for_factors_and_interation){
      pval_interaction <- p.adjust(pval_interaction, method = p_adjust_method)
      pval_factor1 <- p.adjust(pval_factor1, method = p_adjust_method)
      pval_factor2 <- p.adjust(pval_factor2, method = p_adjust_method)
    }
  }

  pval_matrix[sig_feature_indices, 2] <- pval_interaction
  pval_matrix[sig_feature_indices, 3] <- pval_factor1
  pval_matrix[sig_feature_indices, 4] <- pval_factor2
  stat_matrix[sig_feature_indices, 2] <- stat_interaction
  stat_matrix[sig_feature_indices, 3] <- stat_factor1
  stat_matrix[sig_feature_indices, 4] <- stat_factor2

  # Classify follow-up effects
  interaction_class <- rep("Additive", ncol(Z_sig))
  interaction_class[pval_interaction < pval_cutoff_interaction] <- "Interaction"
  names(interaction_class) <- colnames(Z_sig)

  overall_follow_class <- rep(NA, num_features)
  idx_add <- which(interaction_class == "Additive")
  if (length(idx_add) > 0) {
    class_subset <- interaction_class[idx_add]
    class_subset[(pval_factor2[idx_add] <= pval_cutoff_factor2) & (pval_factor1[idx_add] > pval_cutoff_factor1)] <- "Factor2"
    class_subset[(pval_factor2[idx_add] >  pval_cutoff_factor2) & (pval_factor1[idx_add] <= pval_cutoff_factor1)] <- "Factor1"
    interaction_class[idx_add] <- class_subset
  }
  overall_follow_class[sig_feature_indices] <- interaction_class
  class_df$Sig1 <- overall_follow_class

  colnames(pval_matrix) <- colnames(stat_matrix) <- c("Full", "Interaction", "Factor1", "Factor2")
  rownames(pval_matrix) <- rownames(stat_matrix) <- colnames(Z)
  rownames(class_df) <- colnames(Z)

  return(list(pval_matrix = pval_matrix, stat_matrix = stat_matrix, class_df = class_df))
}




#' Build Model Formulas for Three-Factor Differential Analysis
#'
#' This function constructs various model formulas for differential analysis involving three factors.
#' It generates formulas for the full model with all interactions, alternative two-way interaction models,
#' models excluding specific interactions, a null interaction model, and individual main effects models.
#' The formulas vary depending on whether a random effect is included and the specified model fitting function.
#'
#' @param formatted_factors A list containing formatted factor names for the analysis. Expected elements include
#'   \code{factor1_tmp}, \code{factor2_tmp}, \code{factor3_tmp}, and optionally \code{random_effect_tmp}.
#' @param test_func A character string specifying the model fitting function to use (e.g., \code{"lm"}).
#' @param random_effect An optional argument representing the random effect. If \code{NULL}, formulas for fixed effects are generated.
#'
#' @return A list of character strings containing the following formulas:
#'   \item{lm_full}{The full model formula including all three-way interactions.}
#'   \item{lm_int_alt2way}{An alternative model formula including only two-way interactions.}
#'   \item{lm_int_altf1f2}{A model formula excluding the interaction between factor1 and factor2.}
#'   \item{lm_int_altf1f3}{A model formula excluding the interaction between factor1 and factor3.}
#'   \item{lm_int_altf2f3}{A model formula excluding the interaction between factor2 and factor3.}
#'   \item{lm_int_null}{A null model formula with only main effects (no interactions).}
#'   \item{lm_f1}{A model formula for the main effect of factor1 only.}
#'   \item{lm_f2}{A model formula for the main effect of factor2 only.}
#'   \item{lm_f3}{A model formula for the main effect of factor3 only.}
#'
#' @details Depending on whether a random effect is provided, the function creates formulas for either a fixed effects model or a mixed-effects model (including random effects using \code{(1|random_effect)}).
#' @keywords deprecated
#' @examples
#' \dontrun{
#'   # Example for a fixed effects model:
#'   formatted_factors <- list(factor1_tmp = "Group", factor2_tmp = "Treatment",
#'   factor3_tmp = "Time")
#'   formulas <- build_formulas_3F(formatted_factors, test_func = "lm", random_effect = NULL)
#'   print(formulas$lm_full)
#'
#'   # Example for a random effects model:
#'   formatted_factors$random_effect_tmp <- "Subject"
#'   formulas_re <- build_formulas_3F(formatted_factors, test_func = "lm",
#'   random_effect = formatted_factors$random_effect_tmp)
#'   print(formulas_re$lm_full)
#' }
#'
#' @keywords internal
build_formulas_3F <- function(formatted_factors, test_func, random_effect) {
  if (test_func == "lm" && is.null(random_effect)) {
    lm_full <- paste("Y ~ (", formatted_factors$factor1_tmp, ")*(", formatted_factors$factor2_tmp, ")*(", formatted_factors$factor3_tmp, ")", sep = "")
    lm_int_alt2way <- paste("Y ~ (", formatted_factors$factor2_tmp, "):(", formatted_factors$factor3_tmp, ") + (",
                            formatted_factors$factor1_tmp, "):(", formatted_factors$factor3_tmp, ") + (",
                            formatted_factors$factor1_tmp, "):(", formatted_factors$factor2_tmp, ") + (",
                            formatted_factors$factor1_tmp, ") + (", formatted_factors$factor2_tmp, ") + (",
                            formatted_factors$factor3_tmp, ")", sep = "")
    lm_int_altf1f2 <- paste("Y~(", formatted_factors$factor2_tmp, "):(", formatted_factors$factor3_tmp, ")+(",
                            formatted_factors$factor1_tmp, "):(", formatted_factors$factor3_tmp, ")+(",
                            formatted_factors$factor1_tmp, ")+(", formatted_factors$factor2_tmp, ")+(",
                            formatted_factors$factor3_tmp,")",
                            sep = "")

    lm_int_altf1f3 <- paste("Y~(", formatted_factors$factor2_tmp, "):(", formatted_factors$factor3_tmp, ")+(",
                            formatted_factors$factor1_tmp, "):(", formatted_factors$factor2_tmp, ")+(",
                            formatted_factors$factor1_tmp, ")+(", formatted_factors$factor2_tmp, ")+(",
                            formatted_factors$factor3_tmp,")",
                            sep = "")


    lm_int_altf2f3 <- paste("Y~(", formatted_factors$factor1_tmp, "):(", formatted_factors$factor3_tmp, ")+(",
                            formatted_factors$factor1_tmp, "):(", formatted_factors$factor2_tmp, ")+(",
                            formatted_factors$factor1_tmp, ")+(", formatted_factors$factor2_tmp, ")+(",
                            formatted_factors$factor3_tmp,")",
                            sep = "")
    lm_int_null <- paste("Y ~ (", formatted_factors$factor1_tmp, ") + (", formatted_factors$factor2_tmp, ") + (",
                         formatted_factors$factor3_tmp, ")", sep = "")
    lm_f1 <- paste("Y ~ (", formatted_factors$factor1_tmp, ")", sep = "")
    lm_f2 <- paste("Y ~ (", formatted_factors$factor2_tmp, ")", sep = "")
    lm_f3 <- paste("Y ~ (", formatted_factors$factor3_tmp, ")", sep = "")
  } else {
    lm_full <- paste("Y ~ (", formatted_factors$factor1_tmp, ")*(", formatted_factors$factor2_tmp, ")*(", formatted_factors$factor3_tmp,
                     ") + (1|", formatted_factors$random_effect_tmp, ")", sep = "")
    lm_int_alt2way <- paste("Y ~ (", formatted_factors$factor2_tmp, "):(", formatted_factors$factor3_tmp, ") + (",
                            formatted_factors$factor1_tmp, "):(", formatted_factors$factor3_tmp, ") + (",
                            formatted_factors$factor1_tmp, "):(", formatted_factors$factor2_tmp, ") + (",
                            formatted_factors$factor1_tmp, ") + (", formatted_factors$factor2_tmp, ") + (",
                            formatted_factors$factor3_tmp, ") + (1|", formatted_factors$random_effect_tmp, ")", sep = "")
    lm_int_altf1f2 <- paste("Y~(", formatted_factors$factor2_tmp, "):(", formatted_factors$factor3_tmp, ")+(",
                            formatted_factors$factor1_tmp, "):(", formatted_factors$factor3_tmp, ")+(",
                            formatted_factors$factor1_tmp, ")+(", formatted_factors$factor2_tmp, ")+(",
                            formatted_factors$factor3_tmp,")+(1|",formatted_factors$random_effect_tmp,")",sep = "")

    lm_int_altf1f3 <- paste("Y~(", formatted_factors$factor2_tmp, "):(", formatted_factors$factor3_tmp, ")+(",
                            formatted_factors$factor1_tmp, "):(", formatted_factors$factor2_tmp, ")+(",
                            formatted_factors$factor1_tmp, ")+(", formatted_factors$factor2_tmp, ")+(",
                            formatted_factors$factor3_tmp,")+(1|",formatted_factors$random_effect_tmp,")",sep = "")

    lm_int_altf2f3 <- paste("Y~(", formatted_factors$factor1_tmp, "):(", formatted_factors$factor3_tmp, ")+(",
                            formatted_factors$factor1_tmp, "):(", formatted_factors$factor2_tmp, ")+(",
                            formatted_factors$factor1_tmp, ")+(", formatted_factors$factor2_tmp, ")+(",
                            formatted_factors$factor3_tmp,")+(1|",formatted_factors$random_effect_tmp,")",sep = "")
    lm_int_null <- paste("Y ~ (", formatted_factors$factor1_tmp, ") + (", formatted_factors$factor2_tmp, ") + (",
                         formatted_factors$factor3_tmp, ") + (1|", formatted_factors$random_effect_tmp, ")", sep = "")
    lm_f1 <- paste("Y ~ (", formatted_factors$factor1_tmp, ") + (1|", formatted_factors$random_effect_tmp, ")", sep = "")
    lm_f2 <- paste("Y ~ (", formatted_factors$factor2_tmp, ") + (1|", formatted_factors$random_effect_tmp, ")", sep = "")
    lm_f3 <- paste("Y ~ (", formatted_factors$factor3_tmp, ") + (1|", formatted_factors$random_effect_tmp, ")", sep = "")
  }
  return(list(lm_full = lm_full,
              lm_int_alt2way = lm_int_alt2way,
              lm_int_altf1f2 = lm_int_altf1f2,
              lm_int_altf1f3 = lm_int_altf1f3,
              lm_int_altf2f3 = lm_int_altf2f3,
              lm_int_null = lm_int_null,
              lm_f1 = lm_f1,
              lm_f2 = lm_f2,
              lm_f3 = lm_f3))
}



#' Construct Model Formulas for Two- or Three-Factor Analysis
#'
#' This function builds model formula strings for either two or three factors,
#' with optional random effects. When \code{test_func} is \code{"lm"}, the
#' function assumes no random effect (\code{random_effect = NULL}). When
#' \code{test_func} is \code{"lmer"}, a random effect must be provided.
#' Supported factors:
#' \itemize{
#'   \item Two-factor scenario: \code{factor1_tmp} and \code{factor2_tmp}.
#'   \item Three-factor scenario: \code{factor1_tmp}, \code{factor2_tmp}, and \code{factor3_tmp}.
#' }
#'
#' @param formatted_factors A named list containing the formatted factor names
#'   (e.g., \code{factor1_tmp}, \code{factor2_tmp}, \code{factor3_tmp}) and optionally
#'   \code{random_effect_tmp} if using a mixed-effects model. Typically produced by
#'   an internal or auxiliary formatting function.
#' @param test_func A character string indicating which model-fitting approach to use.
#'   Currently supports \code{"lm"} (no random effect) and \code{"lmer"} (mixed-effects).
#' @param random_effect A factor or \code{NULL}. If \code{test_func = "lmer"}, this
#'   should be a valid factor; if \code{test_func = "lm"}, it must be \code{NULL}.
#'
#' @return A named list of formula strings. The exact list structure depends on the
#'   number of non-\code{NULL} factors in \code{formatted_factors}:
#'   \itemize{
#'     \item **Two-factor case**:
#'       \code{list(
#'         full = <full_model_formula>,
#'         int_null = <interaction_null_formula>,
#'         int_alt = <interaction_alt_formula>,
#'         f1 = <factor1_formula>,
#'         f2 = <factor2_formula>
#'       )}
#'
#'     \item **Three-factor case**:
#'       \code{list(
#'         lm_full = <full_model_formula>,
#'         lm_int_alt2way = <two_way_interaction_formula>,
#'         lm_int_altf1f2 = <omitting_f1f2_interaction_formula>,
#'         lm_int_altf1f3 = <omitting_f1f3_interaction_formula>,
#'         lm_int_altf2f3 = <omitting_f2f3_interaction_formula>,
#'         lm_int_null = <no_interaction_formula>,
#'         lm_f1 = <factor1_formula>,
#'         lm_f2 = <factor2_formula>,
#'         lm_f3 = <factor3_formula>
#'       )}
#'   }
#'
#' @details The function calculates \code{num_factors} internally by counting
#'   how many non-\code{NULL} items exist in \code{formatted_factors} (apart from
#'   the random effect) and dividing by 2. This must yield either 2 or 3
#'   (representing two- or three-factor designs). Otherwise, it raises an error.
#'   \enumerate{
#'     \item For \code{test_func = "lm"}, \code{random_effect} must be \code{NULL}.
#'     \item For \code{test_func = "lmer"}, \code{random_effect} must be a valid factor.
#'   }
#'
#' @examples
#' \dontrun{
#' # Example for a two-factor linear model:
#' my_factors <- list(
#'   factor1_tmp = "Group",
#'   factor2_tmp = "Treatment",
#'   random_effect_tmp = NULL  # Not used for 'lm'
#' )
#' formulas_2f <- build_formulas(my_factors, test_func = "lm", random_effect = NULL)
#' print(formulas_2f)
#'
#' # Example for a three-factor mixed-effects model:
#' my_factors_3 <- list(
#'   factor1_tmp = "Group",
#'   factor2_tmp = "Treatment",
#'   factor3_tmp = "Time",
#'   random_effect_tmp = "Subject"  # random effect for 'lmer'
#' )
#' formulas_3f <- build_formulas(my_factors_3, test_func = "lmer", random_effect = factor("Subject"))
#' print(formulas_3f)
#' }
#'
#' @keywords internal
build_formulas <- function(formatted_factors, test_func, random_effect) {
  num_factors <-  (length(formatted_factors)-sum(sapply(formatted_factors,is.null)))/2

  if (num_factors == 2) {
    # Two-factor model
    if (test_func == "lm" && is.null(random_effect)) {
      formula_full          <- paste("Y ~ (", formatted_factors$factor1_tmp, ")*(", formatted_factors$factor2_tmp, ")", sep = "")
      formula_interaction_alt <- paste("Y ~ (", formatted_factors$factor1_tmp, "):(", formatted_factors$factor2_tmp, ")", sep = "")
      formula_interaction_null<- paste("Y ~ (", formatted_factors$factor1_tmp, ")+(", formatted_factors$factor2_tmp, ")", sep = "")
      formula_factor1       <- paste("Y ~ (", formatted_factors$factor1_tmp, ")", sep = "")
      formula_factor2       <- paste("Y ~ (", formatted_factors$factor2_tmp, ")", sep = "")
    } else if (test_func == "lmer" && !is.null(random_effect)) {
      formula_full          <- paste("Y ~ (", formatted_factors$factor1_tmp, ")*(", formatted_factors$factor2_tmp, ") + (1|", formatted_factors$random_effect_tmp, ")", sep = "")
      formula_interaction_alt <- paste("Y ~ (", formatted_factors$factor1_tmp, "):(", formatted_factors$factor2_tmp, ") + (1|", formatted_factors$random_effect_tmp, ")", sep = "")
      formula_interaction_null<- paste("Y ~ (", formatted_factors$factor1_tmp, ")+(", formatted_factors$factor2_tmp, ") + (1|", formatted_factors$random_effect_tmp, ")", sep = "")
      formula_factor1       <- paste("Y ~ (", formatted_factors$factor1_tmp, ") + (1|", formatted_factors$random_effect_tmp, ")", sep = "")
      formula_factor2       <- paste("Y ~ (", formatted_factors$factor2_tmp, ") + (1|", formatted_factors$random_effect_tmp, ")", sep = "")
    } else {
      stop("Invalid model specification: For 'lm', random_effect must be NULL; for 'lmer', random_effect must be provided.")
    }
    return(list(
      full = formula_full,
      int_null = formula_interaction_null,
      int_alt = formula_interaction_alt,

      f1 = formula_factor1,
      f2 = formula_factor2
    ))

  } else if (num_factors == 3) {
    # Three-factor model (same as your existing function)
    if (test_func == "lm" && is.null(random_effect)) {
      lm_full <- paste("Y ~ (", formatted_factors$factor1_tmp, ")*(", formatted_factors$factor2_tmp, ")*(", formatted_factors$factor3_tmp, ")", sep = "")
      lm_int_alt2way <- paste("Y ~ (", formatted_factors$factor2_tmp, "):(", formatted_factors$factor3_tmp, ") + (",
                              formatted_factors$factor1_tmp, "):(", formatted_factors$factor3_tmp, ") + (",
                              formatted_factors$factor1_tmp, "):(", formatted_factors$factor2_tmp, ") + (",
                              formatted_factors$factor1_tmp, ") + (", formatted_factors$factor2_tmp, ") + (",
                              formatted_factors$factor3_tmp, ")", sep = "")
      lm_int_altf1f2 <- paste("Y~(", formatted_factors$factor2_tmp, "):(", formatted_factors$factor3_tmp, ")+(",
                              formatted_factors$factor1_tmp, "):(", formatted_factors$factor3_tmp, ")+(",
                              formatted_factors$factor1_tmp, ")+(", formatted_factors$factor2_tmp, ")+(",
                              formatted_factors$factor3_tmp,")",
                              sep = "")

      lm_int_altf1f3 <- paste("Y~(", formatted_factors$factor2_tmp, "):(", formatted_factors$factor3_tmp, ")+(",
                              formatted_factors$factor1_tmp, "):(", formatted_factors$factor2_tmp, ")+(",
                              formatted_factors$factor1_tmp, ")+(", formatted_factors$factor2_tmp, ")+(",
                              formatted_factors$factor3_tmp,")",
                              sep = "")


      lm_int_altf2f3 <- paste("Y~(", formatted_factors$factor1_tmp, "):(", formatted_factors$factor3_tmp, ")+(",
                              formatted_factors$factor1_tmp, "):(", formatted_factors$factor2_tmp, ")+(",
                              formatted_factors$factor1_tmp, ")+(", formatted_factors$factor2_tmp, ")+(",
                              formatted_factors$factor3_tmp,")",
                              sep = "")
      lm_int_null <- paste("Y ~ (", formatted_factors$factor1_tmp, ") + (", formatted_factors$factor2_tmp, ") + (",
                           formatted_factors$factor3_tmp, ")", sep = "")
      lm_f1 <- paste("Y ~ (", formatted_factors$factor1_tmp, ")", sep = "")
      lm_f2 <- paste("Y ~ (", formatted_factors$factor2_tmp, ")", sep = "")
      lm_f3 <- paste("Y ~ (", formatted_factors$factor3_tmp, ")", sep = "")
    } else {
      lm_full <- paste("Y ~ (", formatted_factors$factor1_tmp, ")*(", formatted_factors$factor2_tmp, ")*(", formatted_factors$factor3_tmp,
                       ") + (1|", formatted_factors$random_effect_tmp, ")", sep = "")
      lm_int_alt2way <- paste("Y ~ (", formatted_factors$factor2_tmp, "):(", formatted_factors$factor3_tmp, ") + (",
                              formatted_factors$factor1_tmp, "):(", formatted_factors$factor3_tmp, ") + (",
                              formatted_factors$factor1_tmp, "):(", formatted_factors$factor2_tmp, ") + (",
                              formatted_factors$factor1_tmp, ") + (", formatted_factors$factor2_tmp, ") + (",
                              formatted_factors$factor3_tmp, ") + (1|", formatted_factors$random_effect_tmp, ")", sep = "")
      lm_int_altf1f2 <- paste("Y~(", formatted_factors$factor2_tmp, "):(", formatted_factors$factor3_tmp, ")+(",
                              formatted_factors$factor1_tmp, "):(", formatted_factors$factor3_tmp, ")+(",
                              formatted_factors$factor1_tmp, ")+(", formatted_factors$factor2_tmp, ")+(",
                              formatted_factors$factor3_tmp,")+(1|",formatted_factors$random_effect_tmp,")",sep = "")

      lm_int_altf1f3 <- paste("Y~(", formatted_factors$factor2_tmp, "):(", formatted_factors$factor3_tmp, ")+(",
                              formatted_factors$factor1_tmp, "):(", formatted_factors$factor2_tmp, ")+(",
                              formatted_factors$factor1_tmp, ")+(", formatted_factors$factor2_tmp, ")+(",
                              formatted_factors$factor3_tmp,")+(1|",formatted_factors$random_effect_tmp,")",sep = "")

      lm_int_altf2f3 <- paste("Y~(", formatted_factors$factor1_tmp, "):(", formatted_factors$factor3_tmp, ")+(",
                              formatted_factors$factor1_tmp, "):(", formatted_factors$factor2_tmp, ")+(",
                              formatted_factors$factor1_tmp, ")+(", formatted_factors$factor2_tmp, ")+(",
                              formatted_factors$factor3_tmp,")+(1|",formatted_factors$random_effect_tmp,")",sep = "")
      lm_int_null <- paste("Y ~ (", formatted_factors$factor1_tmp, ") + (", formatted_factors$factor2_tmp, ") + (",
                           formatted_factors$factor3_tmp, ") + (1|", formatted_factors$random_effect_tmp, ")", sep = "")
      lm_f1 <- paste("Y ~ (", formatted_factors$factor1_tmp, ") + (1|", formatted_factors$random_effect_tmp, ")", sep = "")
      lm_f2 <- paste("Y ~ (", formatted_factors$factor2_tmp, ") + (1|", formatted_factors$random_effect_tmp, ")", sep = "")
      lm_f3 <- paste("Y ~ (", formatted_factors$factor3_tmp, ") + (1|", formatted_factors$random_effect_tmp, ")", sep = "")
    }

    return(list(lm_full = lm_full,
                lm_int_alt2way = lm_int_alt2way,
                lm_int_altf1f2 = lm_int_altf1f2,
                lm_int_altf1f3 = lm_int_altf1f3,
                lm_int_altf2f3 = lm_int_altf2f3,
                lm_int_null = lm_int_null,
                lm_f1 = lm_f1,
                lm_f2 = lm_f2,
                lm_f3 = lm_f3))
  } else {
    stop("Unsupported number of factors. Only 2 or 3 factors are supported.")
  }
}




#' iDAS: Interpretable Differential Analysis of Genes with Three Factors
#'
#' This function implements the iDAS (Interpretable Differential Analysis Signature) framework to identify
#' features associated with three experimental factors (\code{factor1}, \code{factor2}, and \code{factor3}),
#' as well as their interactions. The analysis involves an overall model test, interaction tests
#' (two-way and three-way), and main effects tests. Results include adjusted p-values and test
#' statistics for each feature.
#'
#' @param Z A numeric matrix or data frame where each column represents a feature and each row represents an observation.
#' @param factor1 A factor or vector representing the primary experimental factor.
#' @param factor2 A factor or vector representing the secondary experimental factor.
#' @param factor3 A factor or vector representing the tertiary experimental factor.
#' @param random_effect An optional factor or vector for random effects (e.g., subject ID). Use \code{NULL} if not applicable.
#' @param model_fit_function A character string specifying the model-fitting function (e.g., \code{"lm"} or \code{"lmer"}). Defaults to \code{"lm"}.
#' @param test_function A character string specifying the testing function to use (e.g., \code{"parametric"} or \code{"permutation"}). Defaults to \code{"parametric"}.
#' @param pval_quantile_cutoff A numeric threshold for the quantile-based filtering of overall p-values. Defaults to \code{0.02}.
#' @param pval_cutoff_full A numeric p-value cutoff for the overall model test. Defaults to \code{0.05}.
#' @param pval_cutoff_interaction A numeric p-value cutoff for the omnibus interaction test. Defaults to \code{0.01}.
#' @param pval_cutoff_factor1 A numeric p-value cutoff for testing the main effect of \code{factor1}. Defaults to \code{0.01}.
#' @param pval_cutoff_factor2 A numeric p-value cutoff for testing the main effect of \code{factor2}. Defaults to \code{0.01}.
#' @param pval_cutoff_factor3 A numeric p-value cutoff for testing the main effect of \code{factor3}. Defaults to \code{0.01}.
#' @param pval_cutoff_int12 Numeric p-value cutoff for the two-way interaction between factor1 and factor2. Defaults to \code{0.01}.
#' @param pval_cutoff_int13 Numeric p-value cutoff for the two-way interaction between factor1 and factor3. Defaults to \code{0.01}.
#' @param pval_cutoff_int23 Numeric p-value cutoff for the two-way interaction between factor2 and factor3. Defaults to \code{0.01}.
#' @param pval_cutoff_int123 Numeric p-value cutoff for the three-way interaction among all factors. Defaults to \code{0.01}.
#' @param p_adjust_method A character string specifying the method used to adjust p-values (e.g., \code{"BH"}). Defaults to \code{"BH"}.
#' @param factor1_name Optional label for \code{factor1}.
#' @param factor2_name Optional label for \code{factor2}.
#' @param factor3_name Optional label for \code{factor3}.
#' @param random_effect_name Optional label for the random effect.
#' @param \dots Additional arguments passed to internal functions, model-fitting routines, or test functions.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{pval_matrix}}{A matrix of adjusted p-values for each gene, including main effects and
#'   interactions.}
#'   \item{\code{stat_matrix}}{A matrix of corresponding test statistics.}
#'   \item{\code{class_df}}{A data frame classifying each gene based on the significance of main effects
#'   and interactions.}
#' }
#'
#' @details
#' Internally, the function:
#' \enumerate{
#'   \item Builds the appropriate model formulas for each gene, depending on \code{model_fit_function}
#'     (e.g., \code{lm} vs. \code{lmer}) and whether \code{random_effect} is provided.
#'   \item Performs an overall significance test for each gene (the \code{pval_cutoff_full} threshold).
#'   \item For those genes passing the overall test, conducts an omnibus interaction test and further
#'     specific tests (main effects or two-way/three-way interactions) controlled by the respective
#'     p-value cutoffs.
#' }
#' Multiple testing corrections are applied based on \code{p_adjust_method}.
#'
#' @examples
#' \dontrun{
#' # Generate sample data
#' set.seed(123)
#' Z <- matrix(rnorm(1000), ncol = 10)
#' colnames(Z)=paste0("gene",1:10)
#' factor1 <- as.factor(rep(1:2, each = 5))
#' factor2 <- as.factor(rep(1:2, times = 5))
#' factor3 <- as.factor(rep(1:2, length.out = 10))
#'
#' # Run the differential analysis using iDAS
#' result <- threefactors(
#'   Z, factor1, factor2, factor3,
#'   model_fit_function = "lm",
#'   test_function = "parametric",
#'   pval_quantile_cutoff = 0.02,
#'   pval_cutoff_full = 0.05,
#'   pval_cutoff_interaction = 0.01,
#'   pval_cutoff_factor1 = 0.01,
#'   pval_cutoff_factor2 = 0.01,
#'   pval_cutoff_factor3 = 0.01,
#'   pval_cutoff_int12 = 0.01,
#'   pval_cutoff_int13 = 0.01,
#'   pval_cutoff_int23 = 0.01,
#'   pval_cutoff_int123 = 0.01,
#'   p_adjust_method = "BH"
#' )
#'
#' # Inspect results
#' head(result$pval_matrix)
#' head(result$stat_matrix)
#' head(result$class_df)
#' }
#'
#' @export
threefactors <- function(Z, factor1, factor2, factor3, random_effect = NULL, model_fit_function = "lm",test_function="parametric",
                    pval_quantile_cutoff = 0.02, pval_cutoff_full = 0.05, pval_cutoff_interaction = 0.01,
                    pval_cutoff_factor1 = 0.01, pval_cutoff_factor2 = 0.01, pval_cutoff_factor3 = 0.01,
                    pval_cutoff_int12 = 0.01, pval_cutoff_int13 = 0.01, pval_cutoff_int23 = 0.01,
                    pval_cutoff_int123 = 0.01,
                    p_adjust_method = "BH", factor1_name = NULL, factor2_name = NULL,
                    factor3_name = NULL, random_effect_name = NULL,...) {

  # Initialize matrices for p-values and test statistics (9 columns)
  pval_matrix <- stat_matrix <- matrix(NA, nrow = ncol(Z), ncol = 9)

  # Format factors with the shared current helper.
  formatted_factors <- check_factor_name(factor1_name=factor1_name, factor2_name=factor2_name,
                                  factor3_name=factor3_name, random_effect_name=random_effect_name,
                                  factor1=factor1, factor2=factor2, factor3=factor3, random_effect=random_effect)

  # Build all necessary formulas using helper function
  formulas <- build_formulas(formatted_factors, model_fit_function, random_effect)

  # Select base model formula and model fitting function
  if (model_fit_function == "lm" && is.null(random_effect)) {
    calc0_overall <- "Y ~ 1"
  } else if (model_fit_function == "lmer" && !is.null(random_effect)) {
    calc0_overall <- "Y ~ 1 + (1|random_effect)"
  } else {
    stop("Mismatch: 'lm' requires random_effect=NULL; 'lmer' requires non-NULL random_effect.")
  }
  calc1_overall <- formulas$lm_full

  nfeat <- ncol(Z)
  p_sig <- numeric(nfeat)
  s_sig <- numeric(nfeat)

  # Loop over features: Overall (full) model test
  for (i in 1:nfeat) {
    dat <- if (is.null(random_effect)) {
      data.frame(Y = Z[, i], factor1 = formatted_factors$factor1, factor2 = formatted_factors$factor2, factor3 = formatted_factors$factor3)
    } else {
      data.frame(Y = Z[, i], factor1= formatted_factors$factor1, factor2 = formatted_factors$factor2, factor3 = formatted_factors$factor3, random_effect = formatted_factors$random_effect)
    }
    res <- run_test(dat, calc0_overall,calc1_overall,model_fit_function,test_function,...)
    p_sig[i] <- res[1]
    s_sig[i] <- res[2]
  }
  p_sig_adj <- if (is.null(p_adjust_method)) p_sig else p.adjust(p_sig, method = p_adjust_method)
  pval_matrix[, 1] <- p_sig_adj
  stat_matrix[, 1] <- s_sig

  # Classify overall significance
  if (!is.null(pval_quantile_cutoff)) {
    cutoff <- sort(p_sig_adj)[ceiling(length(p_sig_adj) * pval_quantile_cutoff)]
    overall_class <- ifelse(p_sig_adj <= cutoff, "sig", "non-sig")
  } else {
    overall_class <- ifelse(p_sig_adj <= pval_cutoff_full, "sig", "non-sig")
  }
  class_df <- data.frame(varname = colnames(Z),Sig0 = overall_class)
  idx_sig <- which(class_df$Sig0 == "sig")

  if (length(idx_sig) == 0) {
    warning("No significant features found.")
    colnames(pval_matrix) <- colnames(stat_matrix) <- c("Sig0", "Intornotint", "F1", "F2", "F3",
                                                        "twowaysorthreeways", "F1F2", "F2F3", "F1F3")
    rownames(pval_matrix) <- rownames(stat_matrix) <- class_df$varname
    return(list(pval_matrix = pval_matrix, stat_matrix = stat_matrix, class_df = class_df))
  }

  # Subset Z to significant features
  if (length(idx_sig) > 1) {
    Z_sig <- Z[, idx_sig]
  } else {
    Z_sig <- data.frame(Z[, idx_sig], check.names = FALSE)
    colnames(Z_sig) <- colnames(Z)[idx_sig]
    rownames(Z_sig) <- rownames(Z)
  }

  # Interaction test on significant features: compare lm_int_null vs. lm_full
  n_sig <- ncol(Z_sig)
  p_int <- numeric(n_sig)
  s_int <- numeric(n_sig)
  for (i in 1:n_sig) {
    dat <- if (is.null(random_effect)) {
      data.frame(Y = Z_sig[, i], factor1 = formatted_factors$factor1, factor2 = formatted_factors$factor2, factor3 = formatted_factors$factor3)
    } else {
      data.frame(Y = Z_sig[, i],factor1= formatted_factors$factor1, factor2 = formatted_factors$factor2, factor3 = formatted_factors$factor3, random_effect = formatted_factors$random_effect)
    }
    res <- run_test(dat, formulas$lm_int_null, formulas$lm_full, model_fit_function,test_function,...)
    p_int[i] <- res[1]
    s_int[i] <- res[2]
  }
  p_int_adj <- if (is.null(p_adjust_method)) p_int else p.adjust(p_int, method = p_adjust_method)
  pval_matrix[idx_sig, 2] <- p_int_adj
  stat_matrix[idx_sig, 2] <- s_int
  class_df[idx_sig, "Sig1"] <- ifelse(pval_matrix[idx_sig, 2] > pval_cutoff_interaction, "notInt", "Int")

  # Split significant features into interaction and non-interaction groups
  idx_int <- which(class_df$Sig1 == "Int")
  idx_notint <- which(class_df$Sig1 == "notInt")

  ## Main Effects Tests on Non-interaction Features
  if (length(idx_notint) > 0) {
    if (length(idx_notint) > 1) {
      Z_notint <- Z[, idx_notint]
    } else {
      Z_notint <- data.frame(Z[, idx_notint], check.names = FALSE)
      colnames(Z_notint) <- colnames(Z)[idx_notint]
      rownames(Z_notint) <- rownames(Z)
    }
    n_notint <- ncol(Z_notint)
    p_f1 <- numeric(n_notint); s_f1 <- numeric(n_notint)
    p_f2 <- numeric(n_notint); s_f2 <- numeric(n_notint)
    p_f3 <- numeric(n_notint); s_f3 <- numeric(n_notint)

    for (i in 1:n_notint) {
      dat <- if (is.null(random_effect)) {
        data.frame(Y = Z_notint[, i], factor1=formatted_factors$factor1, factor2 = formatted_factors$factor2, factor3 = formatted_factors$factor3)
      } else {
        data.frame(Y = Z_notint[, i], factor1=formatted_factors$factor1, factor2 = formatted_factors$factor2, factor3 = formatted_factors$factor3, random_effect = formatted_factors$random_effect)
      }
      # Main effect test for factor1: compare lm_f1 vs. lm_int_null
      res1 <- run_test(dat, formulas$lm_f1, formulas$lm_int_null, model_fit_function,test_function,...)
      p_f1[i] <- res1[1]; s_f1[i] <- res1[2]
      # Main effect test for factor2: compare lm_f2 vs. lm_int_null
      res2 <- run_test(dat, formulas$lm_f2, formulas$lm_int_null,model_fit_function,test_function,...)
      p_f2[i] <- res2[1]; s_f2[i] <- res2[2]
      # Main effect test for factor3: compare lm_f3 vs. lm_int_null
      res3 <- run_test(dat, formulas$lm_f3, formulas$lm_int_null,model_fit_function,test_function,...)
      p_f3[i] <- res3[1]; s_f3[i] <- res3[2]
    }
    p_f1_adj <- if (is.null(p_adjust_method)) p_f1 else p.adjust(p_f1, method = p_adjust_method)
    p_f2_adj <- if (is.null(p_adjust_method)) p_f2 else p.adjust(p_f2, method = p_adjust_method)
    p_f3_adj <- if (is.null(p_adjust_method)) p_f3 else p.adjust(p_f3, method = p_adjust_method)

    pval_matrix[idx_notint, 3] <- p_f1_adj
    pval_matrix[idx_notint, 4] <- p_f2_adj
    pval_matrix[idx_notint, 5] <- p_f3_adj
    stat_matrix[idx_notint, 3] <- s_f1
    stat_matrix[idx_notint, 4] <- s_f2
    stat_matrix[idx_notint, 5] <- s_f3

    # Classify non-interaction features based on main effects:
    cls_notint <- rep("Add", length(idx_notint))
    cls_notint[pval_matrix[idx_notint, 3] > pval_cutoff_factor1 & pval_matrix[idx_notint, 4] <= pval_cutoff_factor2 & pval_matrix[idx_notint, 5] <= pval_cutoff_factor3] <- "F1"
    cls_notint[pval_matrix[idx_notint, 3] <= pval_cutoff_factor1 & pval_matrix[idx_notint, 4] > pval_cutoff_factor2 & pval_matrix[idx_notint, 5] <= pval_cutoff_factor3] <- "F2"
    cls_notint[pval_matrix[idx_notint, 3] <= pval_cutoff_factor1 & pval_matrix[idx_notint, 4] <= pval_cutoff_factor2 & pval_matrix[idx_notint, 5] > pval_cutoff_factor3] <- "F3"
    class_df[idx_notint, "notInt"] <- cls_notint
  }

  ## Two-way Tests on Interaction Features
  if (length(idx_int) > 0) {
    if (length(idx_int) > 1) {
      Z_int <- Z[, idx_int]
    } else {
      Z_int <- data.frame(Z[, idx_int], check.names = FALSE)
      colnames(Z_int) <- colnames(Z)[idx_int]
      rownames(Z_int) <- rownames(Z)
    }
    n_int <- ncol(Z_int)
    p_two <- numeric(n_int)
    s_two <- numeric(n_int)
    for (i in 1:n_int) {
      dat <- if (is.null(random_effect)) {
        data.frame(Y = Z_int[, i], factor1 = formatted_factors$factor1, factor2 = formatted_factors$factor2, factor3 = formatted_factors$factor3)
      } else {
        data.frame(Y = Z_int[, i], factor1= formatted_factors$factor1, factor2 = formatted_factors$factor2, factor3 = formatted_factors$factor3, random_effect = formatted_factors$random_effect)
      }
      # Two-way test: compare lm_int_alt2way vs. lm_full
      res_two <- run_test(dat, formulas$lm_int_alt2way, formulas$lm_full, model_fit_function,test_function,...)
      p_two[i] <- res_two[1]
      s_two[i] <- res_two[2]
    }
    p_two_adj <- if (is.null(p_adjust_method)) p_two else p.adjust(p_two, method = p_adjust_method)
    pval_matrix[idx_int, 6] <- p_two_adj
    stat_matrix[idx_int, 6] <- s_two
    class_df[idx_int, "Int"] <- ifelse(pval_matrix[idx_int, 6] < pval_cutoff_int123, "threeway", "twoways")

    # Further two-way tests on features classified as "twoways"
    idx_twoways <- which(class_df$Int == "twoways")
    if (length(idx_twoways) > 0) {
      if (length(idx_twoways) > 1) {
        Z_twoways <- Z[, idx_twoways]
      } else {
        Z_twoways <- data.frame(Z[, idx_twoways], check.names = FALSE)
        colnames(Z_twoways) <- colnames(Z)[idx_twoways]
        rownames(Z_twoways) <- rownames(Z)
      }
      n_twoways <- ncol(Z_twoways)
      p_f1f2 <- numeric(n_twoways); s_f1f2 <- numeric(n_twoways)
      p_f1f3 <- numeric(n_twoways); s_f1f3 <- numeric(n_twoways)
      p_f2f3 <- numeric(n_twoways); s_f2f3 <- numeric(n_twoways)

      for (i in 1:n_twoways) {
        dat <- if (is.null(random_effect)) {
          data.frame(Y = Z_twoways[, i], factor1 = formatted_factors$factor1, factor2 = formatted_factors$factor2, factor3 = formatted_factors$factor3)
        } else {
          data.frame(Y = Z_twoways[, i],  factor1= formatted_factors$factor1, factor2 = formatted_factors$factor2, factor3 = formatted_factors$factor3, random_effect = formatted_factors$random_effect)
        }
        # Test pval_cutoff_int12: compare lm_int_altf1f2 vs. lm_int_alt2way
        res_f1f2 <- run_test(dat, formulas$lm_int_altf1f2, formulas$lm_int_alt2way,model_fit_function,test_function,...)
        p_f1f2[i] <- res_f1f2[1]; s_f1f2[i] <- res_f1f2[2]
        # Test pval_cutoff_int13: compare lm_int_altf1f3 vs. lm_int_alt2way
        res_f1f3 <- run_test(dat, formulas$lm_int_altf1f3, formulas$lm_int_alt2way,model_fit_function,test_function,...)
        p_f1f3[i] <- res_f1f3[1]; s_f1f3[i] <- res_f1f3[2]
        # Test pval_cutoff_int23: compare lm_int_altf2f3 vs. lm_int_alt2way
        res_f2f3 <- run_test(dat, formulas$lm_int_altf2f3, formulas$lm_int_alt2way,model_fit_function,test_function,...)
        p_f2f3[i] <- res_f2f3[1]; s_f2f3[i] <- res_f2f3[2]
      }
      p_f1f2_adj <- if (is.null(p_adjust_method)) p_f1f2 else p.adjust(p_f1f2, method = p_adjust_method)
      p_f1f3_adj <- if (is.null(p_adjust_method)) p_f1f3 else p.adjust(p_f1f3, method = p_adjust_method)
      p_f2f3_adj <- if (is.null(p_adjust_method)) p_f2f3 else p.adjust(p_f2f3, method = p_adjust_method)
      pval_matrix[idx_twoways, 7] <- p_f1f2_adj
      pval_matrix[idx_twoways, 8] <- p_f2f3_adj
      pval_matrix[idx_twoways, 9] <- p_f1f3_adj
      stat_matrix[idx_twoways, 7] <- s_f1f2
      stat_matrix[idx_twoways, 8] <- s_f2f3
      stat_matrix[idx_twoways, 9] <- s_f1f3

      cls_twoways <- rep("twowaycombinations", length(idx_twoways))
      cls_twoways[pval_matrix[idx_twoways, 7] <= pval_cutoff_int12 & pval_matrix[idx_twoways, 8] > pval_cutoff_int23 & pval_matrix[idx_twoways, 9] > pval_cutoff_int13] <- "F1F2"
      cls_twoways[pval_matrix[idx_twoways, 7] > pval_cutoff_int12 & pval_matrix[idx_twoways, 8] <= pval_cutoff_int23 & pval_matrix[idx_twoways, 9] > pval_cutoff_int13] <- "F2F3"
      cls_twoways[pval_matrix[idx_twoways, 7] > pval_cutoff_int12 & pval_matrix[idx_twoways, 8] > pval_cutoff_int23 & pval_matrix[idx_twoways, 9] <= pval_cutoff_int13] <- "F1F3"
      class_df[idx_twoways, "twoways"] <- cls_twoways
    }
  }

  # Set column names and return results
  colnames(pval_matrix) <- colnames(stat_matrix) <- c("Sig0", "Intornotint", "F1", "F2", "F3",
                                                      "twowaysorthreeways", "F1F2", "F2F3", "F1F3")
  rownames(pval_matrix) <- rownames(stat_matrix) <- class_df$varname
  return(list(pval_matrix = pval_matrix, stat_matrix = stat_matrix, class_df = class_df))
}



#' iDAS: Interpretable Differential Abundance Analysis
#'
#' This function implements the iDAS (Interpretable Differential Abundance Analysis Signature)
#' framework for analyzing differential abundance gene signatures. It serves as a
#' comprehensive wrapper supporting both two-factor and three-factor experimental designs.
#' When a third factor is provided, a three-factor analysis is performed via the
#' \code{threefactors} function; otherwise, a two-factor analysis is executed via
#' the \code{twofactors} function.
#'
#' @param Z A numeric matrix or data frame where each column represents a feature (e.g., microbial taxa, metabolites) to be analyzed.
#' @param factor1 A vector or factor representing the first experimental factor.
#' @param factor2 A vector or factor representing the second experimental factor.
#' @param factor3 An optional vector or factor representing the third experimental factor. If provided, a three-factor analysis is performed. Default is \code{NULL}.
#' @param random_effect An optional vector or factor representing a random effect (e.g., subject ID). Default is \code{NULL}.
#' @param model_fit_function A character string indicating the model fitting function to use (e.g., \code{"lm"} for linear models or \code{"lmer"} for mixed-effects models). Default is \code{"lm"}.
#' @param p_adjust_method_for_factors_and_interation Logical or character, specifying whether p-values for factors and interactions should be adjusted. Default is \code{FALSE}.
#' @param pval_quantile_cutoff Numeric value representing the quantile cutoff for overall significance testing. Default is \code{0.02}.
#' @param pval_cutoff_full Numeric p-value cutoff for the overall (full) model test. Default is \code{0.05}.
#' @param pval_cutoff_interaction Numeric p-value cutoff for the interaction test. Default is \code{0.01}.
#' @param pval_cutoff_factor1 Numeric p-value cutoff for testing the main effect of the first factor. Default is \code{0.01}.
#' @param pval_cutoff_factor2 Numeric p-value cutoff for testing the main effect of the second factor. Default is \code{0.01}.
#' @param pval_cutoff_factor3 Numeric p-value cutoff for testing the main effect of the third factor (if \code{factor3} is provided). Default is \code{NULL}.
#' @param p_adjust_method Character string specifying the method for p-value adjustment (e.g., \code{"BH"}). Default is \code{"BH"}.
#' @param factor1_name Optional character string for naming the first factor. Default is \code{NULL}.
#' @param factor2_name Optional character string for naming the second factor. Default is \code{NULL}.
#' @param random_effect_name Optional character string for naming the random effect. Default is \code{NULL}.
#' @param ... Additional arguments passed to the \code{twofactors} or \code{threefactors} functions.
#'
#' @return A list containing the results from the differential abundance analysis. The output typically includes matrices of p-values and test statistics, and a data frame classifying features based on significance.
#'
#' @details The function distinguishes between a two-factor design and a three-factor design based on whether \code{factor3} is provided.
#' It prints \code{"Running three-factor model"} when executing a three-factor analysis and \code{"Running two-factor model"} for a two-factor analysis (for debugging purposes).
#' @import BiocParallel
#'
#' @examples
#' \dontrun{
#'   # Example using two factors
#'   result_two <- iDAS(Z = my_feature_matrix, factor1 = group1, factor2 = group2)
#'
#'   # Example using three factors
#'   result_three <- iDAS(Z = my_feature_matrix, factor1 = group1, factor2 = group2,
#'   factor3 = timepoint)
#' }
#'
#' @export

iDAS <- function(Z, factor1, factor2,factor3=NULL, random_effect = NULL,
                 model_fit_function = "lm",p_adjust_method_for_factors_and_interation=FALSE,
                 pval_quantile_cutoff = 0.02, pval_cutoff_full = 0.05,
                 pval_cutoff_interaction = 0.01, pval_cutoff_factor1 = 0.01, pval_cutoff_factor2 = 0.01,pval_cutoff_factor3=NULL,
                 p_adjust_method = "BH", factor1_name = NULL, factor2_name = NULL, random_effect_name = NULL,...){
  if(!is.null(factor3)){
    print("Running three-factor model")
    result=threefactors(Z,  factor1,
                        factor2,
                        factor3,random_effect = random_effect, model_fit_function = model_fit_function,
                        pval_quantile_cutoff = pval_quantile_cutoff,pval_cutoff_interaction = pval_cutoff_interaction,pval_cutoff_factor1 = pval_cutoff_factor1,
                        pval_cutoff_factor2=pval_cutoff_factor2,pval_cutoff_factor3=pval_cutoff_factor3,p_adjust_method = p_adjust_method,...)
  }else{
    print("Running two-factor model")
    result=twofactors(Z, factor1, factor2, random_effect = random_effect,
               model_fit_function = model_fit_function,p_adjust_method_for_factors_and_interation=p_adjust_method_for_factors_and_interation,
               pval_quantile_cutoff = pval_quantile_cutoff, pval_cutoff_full =pval_cutoff_full,
               pval_cutoff_interaction = pval_cutoff_interaction, pval_cutoff_factor1 = pval_cutoff_factor1, pval_cutoff_factor2 =pval_cutoff_factor2,
               p_adjust_method = p_adjust_method, factor1_name = factor1_name, factor2_name = factor2_name, random_effect_name = random_effect_name,...)
  }
  return(result)
}
