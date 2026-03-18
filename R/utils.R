#' Set Default Name
#'
#' Returns the default name if the provided name is NULL.
#'
#' @param name A character value representing the current name.
#' @param default A character value representing the default name to use if \code{name} is NULL.
#' @return The original \code{name} if not NULL, otherwise \code{default}.
#' @keywords internal
set_default_name <- function(name, default) {
  ifelse(is.null(name), default, name)
}

#' Process Input
#'
#' Converts a factor or vector to a data frame and assigns a specified column name.
#'
#' @param x A factor or vector that should be converted.
#' @param name A character value to assign as the column name.
#' @return A data frame with one column named \code{name} containing \code{x},
#'   or the original \code{x} if it is not a factor or vector.
#' @keywords internal
process_input <- function(x, name) {
  if (is.factor(x) || is.vector(x)) {
    df <- data.frame(x)
    colnames(df) <- name
    return(df)
  }
  return(x)
}

#' Collapse Names
#'
#' Collapses a vector of names into a single string, concatenated with a plus sign.
#'
#' @param name_vec A character vector of names.
#' @return A single string with names separated by "+" if there is more than one name,
#' otherwise returns the single name.
#' @keywords internal
collapse_names <- function(name_vec) {
  if (length(name_vec) > 1) {
    paste(name_vec, collapse = "+")
  } else {
    name_vec
  }
}

#' Check Factor Name
#'
#' Processes input factors and a random effect by converting them into data frames with specified column names.
#' It also collapses the column names into single strings if multiple columns are present.
#' Supports an optional third factor.
#'
#' @param factor1_name A character string representing the name for factor1. Defaults to "factor1" if NULL.
#' @param factor2_name A character string representing the name for factor2. Defaults to "factor2" if NULL.
#' @param factor3_name A character string representing the name for factor3. Defaults to "factor3" if \code{factor3} is provided and factor3_name is NULL.
#' @param random_effect_name A character string representing the name for the random effect. Defaults to "random_effect" if NULL.
#' @param factor1 A factor or vector representing the first factor.
#' @param factor2 A factor or vector representing the second factor.
#' @param factor3 A factor or vector representing the third factor. Optional; defaults to \code{NULL}.
#' @param random_effect A factor or vector representing the random effect. May be \code{NULL}.
#' @return A list containing:
#' \item{factor1_tmp}{Collapsed name for factor1.}
#' \item{factor2_tmp}{Collapsed name for factor2.}
#' \item{factor3_tmp}{Collapsed name for factor3 (or \code{NULL} if \code{factor3} is not provided).}
#' \item{random_effect_tmp}{Collapsed name for the random effect (or \code{NULL} if \code{random_effect} is \code{NULL}).}
#' \item{factor1}{Data frame for factor1.}
#' \item{factor2}{Data frame for factor2.}
#' \item{factor3}{Data frame for factor3 (or \code{NULL}).}
#' \item{random_effect}{Data frame for random_effect (or \code{NULL}).}
check_factor_name <- function(factor1_name, factor2_name, factor3_name = NULL, random_effect_name,
                              factor1, factor2, factor3 = NULL, random_effect=NULL) {
  # Set default names if needed
  factor1_name <- set_default_name(factor1_name, "factor1")
  factor2_name <- set_default_name(factor2_name, "factor2")
  if (!is.null(factor3)) {
    factor3_name <- set_default_name(factor3_name, "factor3")
  }
  random_effect_name <- set_default_name(random_effect_name, "random_effect")

  # Process factor1 and factor2 inputs
  factor1 <- process_input(factor1, factor1_name)
  factor1_name <- colnames(factor1)

  factor2 <- process_input(factor2, factor2_name)
  factor2_name <- colnames(factor2)

  # Process factor3 only if provided
  if (!is.null(factor3)) {
    factor3 <- process_input(factor3, factor3_name)
    factor3_name <- colnames(factor3)
  }

  # Process random_effect only if it's not NULL
  if (!is.null(random_effect)) {
    random_effect <- process_input(random_effect, random_effect_name)
    random_effect_name <- colnames(random_effect)
  }

  # Collapse names if necessary
  factor1_tmp <- collapse_names(factor1_name)
  factor2_tmp <- collapse_names(factor2_name)
  if (!is.null(factor3)) {
    factor3_tmp <- collapse_names(factor3_name)
  } else {
    factor3_tmp <- NULL
  }

  if (!is.null(random_effect)) {
    random_effect_tmp <- collapse_names(random_effect_name)
  } else {
    random_effect_tmp <- NULL
  }

  # Return a list with the processed names and data frames
  return(list(factor1_tmp = factor1_tmp,
       factor2_tmp = factor2_tmp,
       factor3_tmp = factor3_tmp,
       random_effect_tmp = random_effect_tmp,
       factor1 = factor1,
       factor2 = factor2,
       factor3 = factor3,
       random_effect = random_effect))
}
