#' Plotting method for priors
#'
#' Plot prior densities configured via \code{\link{bv_priors}}. The method
#' draws panels for the Minnesota prior first, followed by all configured
#' dummy priors.
#'
#' @param x A \code{bv_priors} object, obtained from \code{\link{bv_priors}}.
#' @param mar Numeric vector. Margins for \code{\link[graphics]{par}}.
#' @param ... Other graphical parameters for \code{\link[graphics]{par}}.
#'
#' @return Returns \emph{x} invisibly.
#'
#' @export
#'
#' @examples
#' priors <- bv_priors(hyper = "auto", sur = bv_sur())
#' plot(priors)
plot.bv_priors <- function(x, mar = c(2, 2, 2, 0.5), ...) {
  if (!inherits(x, "bv_priors")) {
    stop("Please provide a `bv_priors` object.")
  }

  dummy_pos <- names(x) %in% x[["dummy"]]
  dummies <- names(x)[dummy_pos]
  psi_plot_n <- if (!is.null(x[["psi"]])) {
    sum(as.character(x[["psi"]][["mode"]]) != "auto")
  } else {
    0L
  }
  n_plots <- if (!is.null(x[["lambda"]])) {
    2L + psi_plot_n + length(dummies)
  } else {
    0L + length(dummies)
  }

  if (n_plots == 0) {
    return(invisible(x))
  }

  op <- par(
    mfrow = if (n_plots > 4) {
      c(ceiling(n_plots / 2), 2)
    } else {
      c(n_plots, 1)
    },
    mar = mar,
    ...
  )
  if (!is.null(x[["lambda"]])) {
    .plot_bv_minnesota(x)
  }
  if (length(dummies) > 0) {
    for (dummy in dummies) {
      .plot_bv_dummy(x[[dummy]], name = dummy)
    }
  }
  par(op)

  return(invisible(x))
}


#' @noRd
.plot_bv_minnesota <- function(x) {
  .plot_gamma_prior(x[["lambda"]], title = "Minnesota prior: lambda")
  .plot_gamma_prior(x[["alpha"]], title = "Minnesota prior: alpha")
  .plot_bv_psi(x[["psi"]], title = "Minnesota prior: psi")

  return(invisible(x))
}


#' @noRd
.plot_bv_dummy <- function(x, name = NULL) {
  .plot_gamma_prior(
    x,
    title = paste0(
      "Dummy prior",
      if (!is.null(name)) {
        paste0(": ", name)
      }
    )
  )

  return(invisible(x))
}


#' @noRd
.plot_bv_psi <- function(x, title = "psi") {
  mode_chr <- as.character(x[["mode"]])
  if (all(mode_chr == "auto")) {
    message("Skipping psi prior plot(s): modes are set to `auto`.")
    return(invisible(x))
  }

  pos <- which(mode_chr != "auto")
  for (i in pos) {
    .plot_ig_prior(
      mode = x[["mode"]][i],
      min = x[["min"]][i],
      max = x[["max"]][i],
      shape = x[["coef"]][["k"]],
      scale = x[["coef"]][["theta"]],
      title = paste0(title, " #", i)
    )
  }

  return(invisible(x))
}


#' @noRd
.plot_gamma_prior <- function(x, title = NULL) {
  x_seq <- seq(x[["min"]], x[["max"]], length.out = 500L)
  y_seq <- dgamma(
    x_seq,
    shape = x[["coef"]][["k"]],
    scale = x[["coef"]][["theta"]]
  )

  plot(
    x_seq,
    y_seq,
    type = "l",
    xlab = "Value",
    ylab = "Density",
    main = title
  )
  abline(v = c(x[["min"]], x[["max"]]), lty = "dashed", col = "darkgray")
  abline(v = x[["mode"]], lty = "dotted", col = "gray40")

  return(invisible(x))
}


#' @noRd
.plot_ig_prior <- function(mode, min, max, shape, scale, title = NULL) {
  x_seq <- seq(min, max, length.out = 500L)
  y_seq <- exp(p_log_ig(x_seq, shape = shape, scale = scale))

  plot(
    x_seq,
    y_seq,
    type = "l",
    xlab = "Value",
    ylab = "Density",
    main = title
  )
  abline(v = c(min, max), lty = "dashed", col = "darkgray")
  abline(v = mode, lty = "dotted", col = "gray40")

  return(invisible(x_seq))
}
