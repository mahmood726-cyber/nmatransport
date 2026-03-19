#' nmatransport: Population Transportability and Transitivity Assessment in Network Meta-Analysis
#'
#' @keywords internal
"_PACKAGE"

#' @import stats
#' @importFrom utils head tail
#' @importFrom metafor rma
#' @importFrom netmeta netmeta
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_errorbarh theme_minimal labs
#' @importFrom igraph graph_from_data_frame
NULL

# Suppress R CMD check notes
utils::globalVariables(c("Treatment", "Effect", "Lower", "Upper", "from", "to", "weight"))

#' Assess Transitivity in Network Meta-Analysis
#'
#' Evaluates the overlap of covariate distributions (effect modifiers) across different
#' direct comparisons (edges) in a network meta-analysis. Violations of transitivity
#' occur when effect modifiers differ substantially across comparisons.
#'
#' @param data A data frame containing NMA data. Must include 'treat1', 'treat2', and covariates.
#' @param modifiers A character vector of column names representing effect modifiers.
#' @return A list with pairwise distance matrices for covariates and a transitivity score.
#' @export
assess_transitivity <- function(data, modifiers) {
  if (!all(c("treat1", "treat2") %in% names(data))) {
    stop("Data must contain 'treat1' and 'treat2' columns identifying the comparisons.")
  }
  if (!all(modifiers %in% names(data))) {
    stop("Not all specified modifiers are present in the data.")
  }
  
  # Standardize comparison names so A vs B is same as B vs A
  comparisons <- apply(data[, c("treat1", "treat2")], 1, function(x) paste(sort(x), collapse = "_vs_"))
  data$comparison <- comparisons
  
  # Compute mean of modifiers per comparison
  comp_means <- list()
  for (mod in modifiers) {
    comp_means[[mod]] <- tapply(data[[mod]], data$comparison, mean, na.rm = TRUE)
  }
  
  # Calculate distance between comparisons for each modifier
  distances <- list()
  transitivity_score <- 0
  
  unique_comps <- unique(data$comparison)
  if (length(unique_comps) < 2) {
    warning("Only one unique comparison in data. Transitivity assessment requires a network.")
    return(NULL)
  }
  
  for (mod in modifiers) {
    means <- comp_means[[mod]]
    dist_mat <- as.matrix(dist(means))
    distances[[mod]] <- dist_mat
    
    # Simple transitivity violation score: variance of means scaled by overall mean
    transitivity_score <- transitivity_score + var(means, na.rm = TRUE) / (abs(mean(means, na.rm = TRUE)) + 1e-5)
  }
  
  list(
    comparison_means = comp_means,
    distances = distances,
    global_transitivity_violation = transitivity_score,
    n_comparisons = length(unique_comps)
  )
}

#' Compute NMA Transport Weights (Entropy Balancing)
#'
#' Calculates study-level weights to match the overall network to a target
#' population profile, resolving generalizability gaps.
#'
#' @param data A data frame containing NMA data and covariates.
#' @param target_pop A named list or vector of target population covariate means.
#' @return A numeric vector of transport weights for each study.
#' @export
compute_nma_weights <- function(data, target_pop) {
  covariates <- names(target_pop)
  if (!all(covariates %in% names(data))) {
    stop("Target population variables must be present in the data.")
  }
  
  X <- data[, covariates, drop = FALSE]
  
  # Target moments
  target_moments <- as.numeric(target_pop)
  
  # Simple optimization for entropy balancing (simplified version for package)
  init_weights <- rep(1 / nrow(X), nrow(X))
  
  obj_function <- function(lambda) {
    eta <- as.numeric(as.matrix(X) %*% lambda)
    eta <- pmin(eta, 700) # prevent overflow
    w <- init_weights * exp(eta)
    w <- w / sum(w)
    sum((colSums(as.matrix(X) * w) - target_moments)^2)
  }
  
  opt_result <- tryCatch({
    optim(par = rep(0, ncol(X)), fn = obj_function, control = list(maxit = 1000))
  }, error = function(e) NULL)
  
  if (!is.null(opt_result) && opt_result$convergence == 0) {
    weights <- init_weights * exp(as.matrix(X) %*% opt_result$par)
    weights <- as.numeric(weights / sum(weights))
    return(weights)
  } else {
    warning("Weight optimization failed. Returning uniform weights.")
    return(init_weights)
  }
}

#' Run Transported Network Meta-Analysis
#'
#' Runs a network meta-analysis reweighted to represent a target population.
#'
#' @param data Data frame with columns: TE (treatment effect), seTE (standard error), treat1, treat2, studlab.
#' @param weights Numeric vector of study weights (from compute_nma_weights). If NULL, unweighted standard NMA is run.
#' @param sm Summary measure (e.g., "OR", "RR", "MD"). Default is "MD".
#' @return An object of class 'transported_nma' containing the netmeta fit and results.
#' @export
transported_nma <- function(data, weights = NULL, sm = "MD") {
  req_cols <- c("TE", "seTE", "treat1", "treat2", "studlab")
  if (!all(req_cols %in% names(data))) {
    stop(paste("Data must contain all required columns:", paste(req_cols, collapse = ", ")))
  }
  
  # Adjust standard errors based on weights if provided (pseudo-weighting for frequentist NMA)
  # A study with lower weight has a higher effective standard error
  if (!is.null(weights)) {
    if (length(weights) != nrow(data)) stop("Weights vector must match number of rows in data.")
    # Normalize weights to mean 1
    w_norm <- weights / mean(weights)
    # Inflate standard errors for low-weight studies to downweight them in meta-analysis
    data$seTE_adj <- data$seTE / sqrt(w_norm + 1e-8)
  } else {
    data$seTE_adj <- data$seTE
  }
  
  if (!requireNamespace("netmeta", quietly = TRUE)) {
    stop("The 'netmeta' package is required to run transported_nma.")
  }
  
  # Run frequentist NMA
  net_fit <- netmeta::netmeta(TE = data$TE, seTE = data$seTE_adj,
                              treat1 = data$treat1, treat2 = data$treat2,
                              studlab = data$studlab,
                              data = data, sm = sm, 
                              common = FALSE, random = TRUE)
  
  result <- list(
    netmeta_model = net_fit,
    data = data,
    weights_used = !is.null(weights),
    league_table = net_fit$TE.random,
    treatments = net_fit$trts
  )
  class(result) <- "transported_nma"
  return(result)
}

#' Print Transported NMA Results
#' @param x transported_nma object
#' @param ... Additional arguments
#' @export
print.transported_nma <- function(x, ...) {
  cat("
==========================================================
")
  cat("Transported Network Meta-Analysis
")
  cat("==========================================================
")
  if (x$weights_used) {
    cat("Population: Target Population (Weighted)
")
  } else {
    cat("Population: Source Network (Unweighted)
")
  }
  cat("Treatments compared:", length(x$treatments), "
")
  cat("Number of studies:", length(unique(x$data$studlab)), "

")
  
  cat("Heterogeneity (tau):", round(x$netmeta_model$tau, 3), "

")
  
  cat("Treatment Rankings (P-scores):
")
  ranks <- netmeta::netrank(x$netmeta_model)
  print(round(ranks$ranking.random, 3))
  cat("==========================================================
")
  invisible(x)
}
