#' Plotting method for Bayesian VAR forecast error variance decompositions
#'
#' Plotting method for forecast error variance decompositions obtained from
#' \code{\link{fevd.bvar}}. Forecast error variance decompositions of all or a
#' subset of the available variables can be plotted.
#'
#' @param x A \code{bvar_fevd} object, obtained from \code{\link{fevd.bvar}}.
#' @param vars Optional numeric or character vector. Used to subset the plot's
#' forecast error variance decompositions to certain variables by position
#' or name (must be available). Defaults to \code{NULL}, i.e. all variables.
#' @param variables Optional character vector. Names of all variables in the
#' object. Used to subset and title. Taken from \code{x$variables} if available.
#' @param col Character vector. Colours used for shocks (recycled as needed).
#' @param mar Numeric vector. Margins for \code{\link[graphics]{par}}.
#' @param ... Other graphical parameters for \code{\link[graphics]{par}}.
#'
#' @return Returns \emph{x} invisibly.
#'
#' @seealso \code{\link{bvar}}; \code{\link{fevd.bvar}}
#'
#' @keywords BVAR fevd irf analysis plot
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
#' x <- bvar(data, lags = 1, n_draw = 1000L, n_burn = 200L, verbose = FALSE)
#'
#' # Store FEVDs ex-post
#' irf(x) <- irf(x, fevd = TRUE)
#'
#' # Plot FEVDs for all available variables
#' plot(fevd(x))
#' }
plot.bvar_fevd <- function(
  x,
  vars = NULL,
  variables = NULL,
  col = c(
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7"
  ),
  mar = c(3, 3, 1, 0.5),
  ...
) {
  if (!inherits(x, "bvar_fevd")) {
    stop("Please provide a `bvar_fevd` object.")
  }
  plot_fevd(
    x = x,
    vars = vars,
    variables = variables,
    col = col,
    mar = mar,
    ...
  )
}


#' @noRd
plot_fevd <- function(
  x,
  vars = NULL,
  variables = NULL,
  col,
  mar = c(3, 3, 1, 0.5),
  ...
) {
  if (!inherits(x, "bvar") && !inherits(x, "bvar_fevd")) {
    stop("Please provide a `bvar` or `bvar_fevd` object.")
  }
  if (inherits(x, "bvar")) {
    x <- fevd(x)
  }

  M <- dim(x[["fevd"]])[2]
  fevd_mean <- apply(x[["fevd"]], c(2, 3, 4), mean)

  variables <- name_deps(
    variables = if (is.null(variables)) {
      x[["variables"]]
    } else {
      variables
    },
    M = M
  )
  pos <- pos_vars(vars, variables, M)

  .plot_fevd(
    x = fevd_mean,
    variables = variables,
    pos = pos,
    col = col,
    mar = mar,
    mfrow = c(length(pos), 1),
    ...
  )

  return(invisible(x))
}


#' Forecast error variance decomposition plot
#'
#' @param x Numeric array (3-dimensional) with data to plot. The first
#' dimension contains responses, the second paths and the third impulses.
#' @param variables Character vector with the names of variables.
#' @param pos Integer vector. Positions of the forecast error variance
#' decompositions to plot.
#' @param mar Numeric vector. Margins for \code{\link[graphics]{par}}.
#' @param mfrow Numeric vector. Layout for \code{\link[graphics]{par}}.
#' @param ... Other graphical parameters for \code{\link[graphics]{par}}.
#'
#' @importFrom graphics par barplot
#'
#' @noRd
.plot_fevd <- function(
  x,
  variables,
  pos,
  col,
  mar,
  mfrow,
  ...
) {
  M <- dim(x)[1]
  if (M > 5L) {
    M <- 5L
  }

  op <- par(mfrow = mfrow, mar = mar, ...)
  for (i in pos) {
    barplot(
      t(as.matrix(x[i, , ])),
      col = rep(col, length.out = dim(x)[3]),
      border = NA,
      main = paste("Forecast error variance decomposition of", variables[i]),
      legend.text = variables,
      args.legend = list(x = "bottom", inset = c(0, -0.25), bty = "n", ncol = M)
    )
  }
  par(op)
}


#' Plotting method for Bayesian VAR historical decompositions
#'
#' Plot historical decompositions obtained from
#' \code{\link{hist_decomp.bvar}}. Contributions are shown as stacked areas,
#' with an optional line for the reconstructed total.
#'
#' @param x A \code{bvar_hd} object, obtained from
#' \code{\link{hist_decomp.bvar}}.
#' @param vars_response,vars_shock Optional numeric or character vector. Used
#' to subset responses / shocks by position or name.
#' @param include_total Logical scalar. Whether to overlay the sum of all shock
#' contributions as a line.
#' @param total_col Character scalar. Colour of the total line.
#' @param col Character vector. Colours used for shocks (recycled as needed).
#' @param border Character scalar. Border colour for stacked polygons.
#' @param legend Logical scalar. Whether to draw a legend.
#' @param drop_zeros Logical scalar. Whether to drop the first initialized
#' period from the plot.
#' @param mar Numeric vector. Margins for \code{\link[graphics]{par}}.
#' @param ... Other graphical parameters for \code{\link[graphics]{par}}.
#'
#' @return Returns \emph{x} invisibly.
#'
#' @keywords BVAR historical decomposition hd plot
#'
#' @export
#'
#' @importFrom graphics plot.new plot.window axis box polygon lines legend title
plot.bvar_hd <- function(
  x,
  vars_response = NULL,
  vars_shock = NULL,
  include_total = TRUE,
  total_col = "#000000",
  col = c(
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7"
  ),
  border = NA,
  legend = TRUE,
  drop_zeros = TRUE,
  mar = c(3, 3, 1, 0.5),
  ...
) {
  if (!inherits(x, "bvar_hd")) {
    stop("Please provide a `bvar_hd` object.")
  }

  hd <- x[["hd"]]
  M <- dim(hd)[2]
  variables <- name_deps(x[["variables"]], M = M)

  pos_res <- pos_vars(vars_response, variables, M)
  pos_shock <- pos_vars(vars_shock, variables, M)

  if (isTRUE(drop_zeros) && dim(hd)[1] > 1L) {
    hd <- hd[-1, , , drop = FALSE]
  }

  plot_hd(
    x = hd,
    variables = variables,
    pos_res = pos_res,
    pos_shock = pos_shock,
    include_total = include_total,
    total_col = total_col,
    col = col,
    border = border,
    legend = legend,
    mar = mar,
    ...
  )

  return(invisible(x))
}


#' @noRd
plot_hd <- function(
  x,
  variables,
  pos_res,
  pos_shock,
  include_total = TRUE,
  total_col,
  col,
  border = NA,
  legend = TRUE,
  mar = c(3, 3, 1, 0.5),
  ...
) {
  .plot_hd(
    x = x,
    variables = variables,
    pos_res = pos_res,
    pos_shock = pos_shock,
    include_total = include_total,
    total_col = total_col,
    col = col,
    border = border,
    legend = legend,
    mar = mar,
    ...
  )

  return(invisible(x))
}


#' @noRd
.plot_hd <- function(
  x,
  variables,
  pos_res,
  pos_shock,
  include_total = TRUE,
  total_col = "#1A1A1A",
  col,
  border = NA,
  legend = TRUE,
  mar = c(2, 2, 2, 0.5),
  ...
) {
  Tn <- dim(x)[1]
  t_index <- seq_len(Tn)
  shock_cols <- rep(col, length.out = length(pos_shock))

  op <- par(mfrow = c(length(pos_res), 1), mar = mar, ...)
  for (i in pos_res) {
    contrib <- x[, i, pos_shock, drop = FALSE]
    contrib <- matrix(contrib, nrow = Tn, ncol = length(pos_shock))
    total <- rowSums(contrib)

    pos_part <- pmax(contrib, 0)
    neg_part <- pmin(contrib, 0)

    ymax <- max(
      0,
      max(rowSums(pos_part), na.rm = TRUE),
      max(total, na.rm = TRUE)
    )
    ymin <- min(
      0,
      min(rowSums(neg_part), na.rm = TRUE),
      min(total, na.rm = TRUE)
    )

    plot.new()
    plot.window(xlim = c(1, Tn), ylim = c(ymin, ymax))
    axis(1)
    axis(2)
    grid()
    box()

    base <- rep(0, Tn)
    for (j in seq_along(pos_shock)) {
      y <- pos_part[, j]
      top <- base + y
      polygon(
        x = c(t_index, rev(t_index)),
        y = c(base, rev(top)),
        col = shock_cols[j],
        border = border
      )
      base <- top
    }

    base <- rep(0, Tn)
    for (j in seq_along(pos_shock)) {
      y <- neg_part[, j]
      bot <- base + y
      polygon(
        x = c(t_index, rev(t_index)),
        y = c(base, rev(bot)),
        col = shock_cols[j],
        border = border
      )
      base <- bot
    }

    lines(t_index, rep(0, Tn), lty = "dashed", col = "gray")
    if (include_total) {
      lines(t_index, total, col = total_col, lwd = 2)
    }

    title(main = paste("Historical decomposition of", variables[i]))

    if (legend) {
      legend(
        "topright",
        legend = c(variables[pos_shock]),
        fill = c(shock_cols),
        border = NA,
        bty = "n"
      )
    }
  }
  par(op)

  return(invisible(x))
}
