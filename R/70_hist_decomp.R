#' Historical decomposition
#'
#' Function to compute a historical variance decomposition of a VAR.
#'
#' @param x,object A \code{bvar} object, obtained from \code{\link{bvar}}.
#' @param type Character scalar. Whether to use mean or quantile values.
#' @param vars_response,vars_shock Optional numeric or character vector. Used
#' to subset the summary method's outputs to certain response / shock variables
#' by position or name (must be available). Defaults to \code{NULL}, i.e. all
#' variables.
#' @param ... Not used.
#'
#' @return Returns a list of class \code{bvar_hd} with elements \code{hd}
#' (a numerical array of dimensions time, variable, shock), \code{variables}
#' and \code{type}. The summary method returns a list of class
#' \code{bvar_hd_summary}, which can be used with \code{print}.
#'
#' @keywords BVAR historical decomposition hd
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Access a subset of the fred_qd dataset
#' data <- fred_qd[, c("CPIAUCSL", "UNRATE", "FEDFUNDS")]
#' # Transform it to be stationary
#' data <- fred_transform(data, codes = c(5, 5, 1), lag = 4)
#'
#' # Estimate a BVAR using one lag, default settings and very few draws
#' x <- bvar(data, lags = 1, n_draw = 600L, n_burn = 100L, verbose = FALSE)
#'
#' # Compute historical decomposition
#' hist_decomp(x, type = "mean")
#'
#' # Summarise the contributions to the first variable
#' summary(hist_decomp(x, type = "mean"), vars_response = 1)
#' }
hist_decomp.bvar <- function(x, type = c("mean", "quantile"), ...) {
  type <- match.arg(type)

  # Necessary quantities
  shock <- matrix(0, x[["meta"]][["K"]] - 1, x[["meta"]][["M"]])
  shock[seq(x[["meta"]][["M"]]), ] <- t(chol(vcov(x, type = type)))
  comp <- companion(x, type = type)
  eps <- solve(shock[seq(x[["meta"]][["M"]]), ], t(residuals(x, type = type)))

  hd <- compute_hd(x, shock, comp, eps)

  M <- x[["meta"]][["M"]]
  vars <- name_deps(x[["variables"]], M = M)
  time_names <- rownames(x[["meta"]][["Y"]])
  if (is.null(time_names) || length(time_names) != dim(hd)[1]) {
    time_names <- as.character(seq_len(dim(hd)[1]))
  }

  dimnames(hd) <- list(
    "time" = time_names,
    "variable" = vars,
    "shock" = vars
  )

  out <- structure(
    list(
      "hd" = hd,
      "variables" = vars,
      "type" = type
    ),
    class = "bvar_hd"
  )

  return(out)
}


#' @rdname hist_decomp.bvar
#' @export
hist_decomp <- function(x, ...) {
  UseMethod("hist_decomp", x)
}


#' @rdname hist_decomp.bvar
#' @export
hist_decomp.default <- function(x, ...) {
  stop("No methods for class ", paste0(class(x), collapse = " / "), " found.")
}


#' @export
print.bvar_hd <- function(x, ...) {
  cat(
    "Numeric array (dimensions ",
    paste0(dim(x[["hd"]]), collapse = ", "),
    ") of historical decomposition values from a BVAR.\n",
    sep = ""
  )
  cat("Type: ", x[["type"]], "\n", sep = "")
  cat("Variables: ", paste(x[["variables"]], collapse = ", "), "\n", sep = "")

  return(invisible(x))
}


#' @rdname hist_decomp.bvar
#' @export
summary.bvar_hd <- function(
  object,
  vars_response = NULL,
  vars_shock = NULL,
  ...
) {
  if (!inherits(object, "bvar_hd")) {
    stop("Please provide a `bvar_hd` object.")
  }

  hd <- object[["hd"]]
  M <- dim(hd)[2]
  variables <- name_deps(object[["variables"]], M = M)
  pos_res <- pos_vars(vars_response, variables, M)
  pos_shock <- pos_vars(vars_shock, variables, M)

  out <- structure(
    list(
      "hd" = object,
      "variables" = variables,
      "pos_res" = pos_res,
      "pos_shock" = pos_shock
    ),
    class = "bvar_hd_summary"
  )

  return(out)
}


#' @noRd
#' @export
print.bvar_hd_summary <- function(x, digits = 2L, ...) {
  if (!inherits(x, "bvar_hd_summary")) {
    stop("Please provide a `bvar_hd_summary` object.")
  }

  print.bvar_hd(x[["hd"]])

  hd <- x[["hd"]][["hd"]]
  time_names <- dimnames(hd)[["time"]]
  if (is.null(time_names)) {
    time_names <- as.character(seq_len(dim(hd)[1]))
  }
  shocks <- x[["variables"]][x[["pos_shock"]]]

  for (i in x[["pos_res"]]) {
    cat("\nHistorical decomposition of ", x[["variables"]][i], ":\n", sep = "")
    contrib <- hd[, i, x[["pos_shock"]], drop = FALSE]
    contrib <- matrix(
      contrib,
      nrow = dim(hd)[1],
      ncol = length(x[["pos_shock"]]),
      dimnames = list(time_names, shocks)
    )
    print(round(contrib, digits = digits))
  }

  return(invisible(x))
}


#' Historical variance decompostions draws
#'
#' Computes a historical variance decomposition (HD) using draws of the
#' orthogonalized vcov-matrix (currently retrieved with Cholesky decomposition).
#'
#' @param x Numeric matrix. Contains a draw of impulse responses
#' obtained from \code{\link{compute_irf}}.
#' @param shock Numeric matrix with the shocks.
#' @param comp Numeric matrix with the coefficients in companion form.
#' @param eps Numeric matrix of the inverse shock times residuals.
#' @param lags,N,M,K Integer scalar. Dimensions of the VAR.
#'
#' @return Returns a numeric array of HDs.
#'
#' @noRd
compute_hd <- function(
  x,
  shock,
  comp,
  eps,
  lags = x[["meta"]][["lags"]],
  N = x[["meta"]][["N"]],
  M = x[["meta"]][["M"]],
  K = x[["meta"]][["K"]]
) {
  kind_of_I <- cbind(diag(M), matrix(0, M, (lags - 1) * M))
  hd_decomp <- array(0, dim = c(M, N + 1, M))
  tmp <- array(0, dim = c(K - 1, N + 1, M))

  for (j in seq_len(dim(hd_decomp)[3])) {
    # Over variables
    eps_tmp <- matrix(0, M, N + 1)
    eps_tmp[j, seq(2, N + 1)] <- eps[j, ]

    for (i in seq.int(2, dim(hd_decomp)[2])) {
      # Over observations
      tmp[, i, j] <- shock %*% eps_tmp[, i] + comp %*% tmp[, i - 1, j]
      hd_decomp[, i, j] <- kind_of_I %*% tmp[, i, j]
    }
  }

  # (variable, time, shock) -> (time, variable, shock); shock stays last to
  # match the ordering used for impulse responses and forecasts
  return(aperm(hd_decomp, c(2, 1, 3)))
}
