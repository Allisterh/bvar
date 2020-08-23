#' Function for tidying up Bayesian VAR impulse responses
#'
#' Function for tidy summarizing of impulse responses obtained from
#' \code{\link{irf.bvar}} in tibble format.
#' Impulse responses of all or a subset of the available variables can be
#' summarized.
#'
#' @param x A \code{bvar_irf} object, obtained from \code{\link{irf.bvar}}
#' or a \code{bvar} object, obtained from \code{\link{bvar}}. In the latter
#' case impulse reponses will be calculated ex-post using arguments from the
#' ellipsis.
#' @param vars_impulse,vars_response Optional numeric or character vector. Used
#' to subset the impulses / responses to certain variables by position or name
#' (must be available). Defaults to \code{NULL}, i.e. all variables.
#' @param variables Optional character vector. Names of all variables in the
#' object. Used to subset. Taken from \code{x$variables} if available.
#' @param ... Other parameters for calculating impulse responses if \emph{x}
#' is not a \code{bvar_irf} object.
#'
#' @return Returns a \code{\link[tibble]{tibble}} with relevant information
#' for further processing.
#'
#' @seealso \code{\link{bvar}}; \code{\link{irf.bvar}}
#'
#' @keywords BVAR irf impulse responses tidy
#'
#' @export
tidy.bvar_irf <- function(
  x,
  vars_response = NULL,
  vars_impulse = NULL,
  variables = NULL,
  ...) {

  if(!inherits(x, "bvar") && !inherits(x, "bvar_irf")) {
    stop("Please provide a `bvar` or `bvar_irf` object.")
  }

  if(inherits(x, "bvar")) {
    if(is.null(x[["irf"]])) {message("No IRFs found. Calculating...")}
    x <- irf(x, ...)
  }

  df <- gg_df_irf(x)
  df <- tidyr::pivot_wider(df, names_from = quant)
  M <- length(levels(df[["response"]]))

  variables <- name_deps(variables = if(is.null(variables)) {
    x[["variables"]]} else {variables}, M = M)

  pos_imp <- pos_vars(vars_impulse, variables, M)
  pos_res <- pos_vars(vars_response, variables, M)

  df <- df[intersect(which(df[["impulse"]] %in% variables[pos_imp]),
                     which(df[["response"]] %in% variables[pos_res])), ]

  return(tibble::as_tibble(df))
}


#' Function for tidying up Bayesian VAR forecast error variance decompositions
#'
#' Function for tidy summarizing of forecast error variance decompositions
#' obtained from \code{\link{fevd.bvar}} in tibble format.
#' Forecast error variance decompositions of all or a subset of the available
#' variables can be summarized.
#'
#' @param x A \code{bvar_fevd} object, obtained from \code{\link{fevd.bvar}}
#' or either a \code{bvar} object, obtained from \code{\link{bvar}} or a
#' \code{bvar_irf} object, obtained from \code{\link{irf.bvar}}. For the latter
#' cases forecast error variance decompositions will be calculated ex-post using
#' arguments from the ellipsis.
#' @param vars_impulse,vars_response Optional numeric or character vector. Used
#' to subset the impulses / responses to certain variables by position or name
#' (must be available). Defaults to \code{NULL}, i.e. all variables.
#' @param variables Optional character vector. Names of all variables in the
#' object. Used to subset. Taken from \code{x$variables} if available.
#' @param ... Other parameters for calculating impulse responses if \emph{x}
#' is not a \code{bvar_irf} object.
#'
#' @return Returns a \code{\link[tibble]{tibble}} with relevant information
#' for further processing.
#'
#' @seealso \code{\link{bvar}}; \code{\link{irf.bvar}}
#'
#' @keywords BVAR irf impulse responses tidy
#'
#' @export
tidy.bvar_fevd <- function(
  x,
  vars_response = NULL,
  vars_impulse = NULL,
  variables = NULL,
  ...) {

  if(!inherits(x, "bvar") && !inherits(x, "bvar_irf") && !inherits(x, "bvar_fevd")) {
    stop("Please provide a `bvar`, `bvar_irf` or `bvar_fevd` object.")
  }

  if(inherits(x, "bvar")) {
    if(is.null(x[["irf"]][["fevd"]])) {
      message("No FEVDs found. Calculating...")
    }
    x <- irf(x, ..., fevd = TRUE)[["fevd"]]
  }

  if(inherits(x, "bvar_irf")) {
    if(is.null(x[["fevd"]])) {
      message("No FEVDs found. Calculating...")
    }
    x <- fevd(x, ...)
  }

  df <- gg_df_irf(x)
  df <- tidyr::pivot_wider(df, names_from = quant)
  M <- length(levels(df[["response"]]))

  variables <- name_deps(variables = if(is.null(variables)) {
    x[["variables"]]} else {variables}, M = M)

  pos_imp <- pos_vars(vars_impulse, variables, M)
  pos_res <- pos_vars(vars_response, variables, M)

  df <- df[intersect(which(df[["impulse"]] %in% variables[pos_imp]),
                     which(df[["response"]] %in% variables[pos_res])), ]
  names(df)[1:2] <- c("Contribution", "Variation")

  return(tibble::as_tibble(df))
}


#' Function for tidying up Bayesian VAR predictions
#'
#' Function for tidy summarizing of predictions obtained from
#' \code{\link{predict.bvar}} in tibble format.
#' Predictions of all or a subset of the available variables can be
#' summarized.
#'
#' @param x A \code{bvar_fcast} object, obtained from \code{\link{predict.bvar}}
#' or a \code{bvar} object, obtained from \code{\link{bvar}}. In the latter
#' case predictions will be calculated ex-post using arguments from the
#' ellipsis.
#' @param vars Optional numeric or character vector. Used to subset the set
#' of predictions to certain variables by position or name (must be available).
#' Defaults to \code{NULL}, i.e. all variables.
#' @param variables Optional character vector. Names of all variables in the
#' object. Used to subset. Taken from \code{x$variables} if available.
#' @param ... Other parameters for calculating impulse responses if \emph{x}
#' is not a \code{bvar_fcast} object.
#'
#' @return Returns a \code{\link[tibble]{tibble}} with relevant information
#' for further processing.
#'
#' @seealso \code{\link{bvar}}; \code{\link{predict.bvar}}
#'
#' @keywords BVAR predictions forecasts tidy data
#'
#' @export
tidy.bvar_fcast <- function(
  x,
  vars = NULL,
  variables = NULL,
  ...) {

  if(!inherits(x, "bvar") && !inherits(x, "bvar_fcast")) {
    stop("Please provide a `bvar` or `bvar_fcast` object.")
  }

  if(inherits(x, "bvar")) {
    if(is.null(x[["fcast"]])) {message("No forecasts found. Calculating...")}
    x <- predict(x, ...)
  }

    df <- gg_df_fcast(x, t_back = 0)
  df <- tidyr::pivot_wider(df, names_from = quant)

  M <- length(levels(df[["vars"]]))

  variables <- name_deps(variables = if(is.null(variables)) {
    x[["variables"]]} else {variables}, M = M)

  pos <- pos_vars(vars, variables, M)

  df <- df[which(df[["vars"]] %in% variables[pos]), ]

  return(tibble::as_tibble(df))
}


#' Function for tidying up Bayesian VAR outputs
#'
#' Function for tidy summarizing output obtained from \code{\link{bvar}}
#' in tibble format.
#' All or a subset of the available variables can be summarized.
#'
#' @param x A \code{bvar} object, obtained from \code{\link{bvar}}.
#' @param vars Character vector used to select variables. Elements are matched
#' to hyperparameters or coefficients. Coefficients may be matched based on
#' the dependent variable (by providing the name or position) or the
#' explanatory variables (by providing the name and the desired lag). See the
#' example section for a demonstration. Defaults to \code{NULL}, i.e. all
#' hyperparameters.
#' @param vars_impulse,vars_response Optional character or integer vectors used
#' to select coefficents. Dependent variables are specified with
#' \emph{vars_response}, explanatory ones with \emph{vars_impulse}. Default to
#' \code{NULL}, indicating that no coefficients will be processed.
#' @param variables Optional character vector. Names of all variables in the
#' object. Used to subset and title. Taken from \code{x$variables} if available.
#' @param quants Optional numeric vector for computing quantiles of parameter
#' draws.
#' @param chains List of \code{bvar} objects. Contents are then added to trace
#' and density plots to help assessing covergence.
#'
#' @return Returns a \code{\link[tibble]{tibble}} with relevant information
#' for further processing.
#'
#' @seealso \code{\link{bvar}}
#'
#' @keywords BVAR MCMC analysis tidy data
#'
#' @export
tidy.bvar <- function(
  x,
  vars = NULL,
  vars_response = NULL,
  vars_impulse = NULL,
  variables = NULL,
  quants = NULL,
  chains = list()) {

  if(!inherits(x, "bvar")) {
    if(inherits(x[[1]], "bvar")) { # Allow chains to x
      chains <- x
      x <- x[[1]]
      chains[[1]] <- NULL
    } else {stop("Please provide a `bvar` object.")}
  }

  if(inherits(chains, "bvar")) {chains <- list(chains)}
  lapply(chains, function(x) {if(!inherits(x, "bvar")) {
    stop("Please provide `bvar` objects to the chains parameter.")
  }})

  if(!is.null(quants)) {
    quants <- quantile_check(quants)
  }

  prep <- prep_data(x, vars = vars, vars_response = vars_response,
                    vars_impulse = vars_impulse,
                    chains = chains, check_chains = FALSE)

  df <- gg_df_bvar(prep)

  out <- expand.grid(unique(df[["vars"]]), unique(df[["chain"]]))
  names(out) <- c("vars", "chain")
  if(!is.null(quants) && quants != 0.5) {
    quants <- quantile_check(quants)
    out <- cbind(out, do.call("rbind",
                              tapply(X = df[["value"]],
                                     INDEX = list(df[["vars"]],
                                                  df[["chain"]]),
                                     FUN = quantile, probs = quants)))
  } else {
    out <- cbind(out, as.vector(tapply(X = df[["value"]],
                                       INDEX = list(df[["vars"]],
                                                    df[["chain"]]),
                                       FUN = quantile, probs = 0.5)))
    names(out)[3] <- "50%"
  }

  return(tibble::as_tibble(out))
}
