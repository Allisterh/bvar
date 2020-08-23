#' \code{\link[ggplot2]{ggplot}} plotting function for Bayesian VAR impulse responses
#'
#' Plotting function for impulse responses obtained from \code{\link{irf.bvar}}
#' using \code{\link[ggplot2]{ggplot}}.
#' Impulse responses of all or a subset of the available variables can be
#' plotted.
#'
#' @param x A \code{bvar_irf} object, obtained from \code{\link{irf.bvar}}
#' or a \code{bvar} object, obtained from \code{\link{bvar}}. In the latter
#' case impulse reponses will be calculated ex-post using arguments from the
#' ellipsis.
#' @param vars_impulse,vars_response Optional character or integer vectors used
#' to select coefficents. Dependent variables are specified with
#' \emph{vars_response}, explanatory ones with \emph{vars_impulse}. Default to
#' \code{NULL}, indicating that no coefficients will be displayed.
#' @param col Character vector. Colour(s) of the lines delineating credible
#' intervals. Single values will be recycled if necessary. Recycled HEX color
#' codes are varied in transparency if not provided (e.g. "#737373FF").
#' @param area Logical scalar. Whether to fill the credible intervals.
#' @param fill Character vector. Colour(s) to fill the credible intervals with.
#' See \emph{col} for more information.
#' @param variables Optional character vector. Names of all variables in the
#' object. Used to subset and title. Taken from \code{x$variables} if available.
#' @param scales String indicating whether the scales of the plotted impulse
#' responses are \code{"free"}, i.e. each plot has its own y-axis, or
#' \code{"fixed"}, i.e. the y-axis is common across plots. Defaults to
#' \code{"fixed"}.
#' @param ... Other parameters for calculating impulse responses if \emph{x}
#' is not a \code{bvar_irf} object.
#'
#' @return Returns an object of class \code{\link[ggplot2]{ggplot}}.
#'
#' @seealso \code{\link{bvar}}; \code{\link{irf.bvar}}
#'
#' @keywords BVAR irf fevd analysis ggplot
#'
#' @export
irf_plot_gg <- function(
  x,
  vars_response = NULL,
  vars_impulse = NULL,
  col = "#737373",
  area = FALSE,
  fill = "#808080",
  variables = NULL,
  scales = c("fixed", "free"),
  ...) {

  if(!inherits(x, "bvar") && !inherits(x, "bvar_irf")) {
    stop("Please provide a `bvar` or `bvar_irf` object.")
  }

  if(inherits(x, "bvar")) {
    if(is.null(x[["irf"]])) {message("No IRFs found. Calculating...")}
    x <- irf(x, ...)
  }

  scales <- match.arg(scales)

  df <- gg_df_irf(x)
  M <- length(levels(df[["response"]]))
  P_n <- levels(df[["quant"]])
  P <- length(P_n)
  df <- tidyr::pivot_wider(df, names_from = quant)

  variables <- name_deps(variables = if(is.null(variables)) {
    x[["variables"]]} else {variables}, M = M)

  pos_imp <- pos_vars(vars_impulse, variables, M)
  pos_res <- pos_vars(vars_response, variables, M)

  df <- df[intersect(which(df[["impulse"]] %in% variables[pos_imp]),
                     which(df[["response"]] %in% variables[pos_res])), ]

  df[["impulse"]] <- factor(paste("Shock:", df[["impulse"]]))
  df[["response"]] <- factor(paste("Response of:", df[["response"]]))

  p <- ggplot(df, aes(x = time - 1)) +
    scale_x_continuous(breaks = function(x) {
      unique(floor(pretty(seq(0, (max(x)) * 1.1))))
    }) +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(strip.background = element_blank())

  if(scales == "free") {
    p <- p + facet_wrap(response ~ impulse, scales = "free_y")
  } else {
    p <- p + facet_grid(response ~ impulse)
  }

  if(P > 1) {
    if(area) {
      fill <- fill_ci_col(x = integer(), y = fill, P = P)
      for(pp in seq_len(P - 1)) {
        p <- p + geom_ribbon(aes(ymin = .data[[P_n[pp]]],
                                 ymax = .data[[P_n[pp + 1]]]),
                             fill = fill[pp])
      }
    } else {
      col <- fill_ci_col(x = integer(), y = col, P = P)
      names(col) <- P_n[which(P_n != "50%")]
      for(pp in P_n[which(P_n != "50%")]) {
        p <- p + geom_line(aes(y = .data[[pp]]),
                           color = col[pp])
      }
    }
  }

  p <- p +
    geom_hline(yintercept = 0, colour = "red", lty = 2) +
    geom_line(aes(y = `50%`))

  return(p)
}


#' \code{\link[ggplot2]{ggplot}} plotting function for Bayesian VAR forecast
#' error variance decompositions
#'
#' Plotting function for forecast error variance decompositions obtained from
#' \code{\link{fevd.bvar}} using \code{\link[ggplot2]{ggplot}}.
#' Forecast error variance decompositions of all or a subset of the available
#' variables can be plotted.
#'
#' @param x A \code{bvar_fevd} object, obtained from \code{\link{fevd.bvar}}
#' or either a \code{bvar} object, obtained from \code{\link{bvar}} or a
#' \code{bvar_irf} object, obtained from \code{\link{irf.bvar}}. For the latter
#' cases forecast error variance decompositions will be calculated ex-post using
#' arguments from the ellipsis.
#' @param vars_impulse,vars_response Optional numeric or character vector. Used
#' to subset the plot's impulses / responses to certain variables by position
#' or name (must be available). Defaults to \code{NULL}, i.e. all variables.
#' @param fill Character vector. Colour(s) to fill the bars.
#' @param variables Optional character vector. Names of all variables in the
#' object. Used to subset and title. Taken from \code{x$variables} if available.
#' @param ... Other parameters for calculating forecast error variance
#' decompositions if \emph{x} is not a \code{bvar_fevd} object.
#'
#' @return Returns an object of class \code{\link[ggplot2]{ggplot}}.
#'
#' @seealso \code{\link{bvar}}; \code{\link{fevd.bvar}}; \code{\link{irf.bvar}}
#'
#' @keywords BVAR irf fevd analysis ggplot
#'
#' @export
fevd_plot_gg <- function(
  x,
  vars_response = NULL,
  vars_impulse = NULL,
  fill = "#808080",
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
  M <- length(levels(df[["response"]]))
  P_n <- levels(df[["quant"]])
  P <- length(P_n)
  df <- tidyr::pivot_wider(df, names_from = quant)

  variables <- name_deps(variables = if(is.null(variables)) {
    x[["variables"]]} else {variables}, M = M)

  pos_imp <- pos_vars(vars_impulse, variables, M)
  pos_res <- pos_vars(vars_response, variables, M)

  df <- df[intersect(which(df[["impulse"]] %in% variables[pos_imp]),
                     which(df[["response"]] %in% variables[pos_res])), ]

  df[["response"]] <- factor(paste("Variation in:", df[["response"]]))

  p <- ggplot(df, aes(x = time - 1)) +
    scale_x_continuous(breaks = function(x) {
      unique(floor(pretty(seq(0, (max(x)) * 1.1))))
    }) +
    scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +
    labs(x = NULL, y = NULL, fill = "Contribution of:") +
    theme_bw() +
    theme(strip.background = element_blank()) +
    facet_grid(response~.) +
    geom_col(aes(y = `50%`, fill = impulse, group = impulse)) +
    theme(legend.position = "bottom")

  # if(P > 1) {
  #   if(area) {
  #     fill <- fill_ci_col(x = integer(), y = fill, P = P)
  #     for(pp in seq_len(P - 1)) {
  #       p <- p + geom_ribbon(aes(ymin = .data[[P_n[pp]]],
  #                                ymax = .data[[P_n[pp + 1]]]),
  #                            fill = fill[pp])
  #     }
  #   } else {
  #     col <- fill_ci_col(x = integer(), y = col, P = P)
  #     names(col) <- P_n[which(P_n != "50%")]
  #     for(pp in P_n[which(P_n != "50%")]) {
  #       p <- p + geom_line(aes(y = .data[[pp]]),
  #                          color = col[pp])
  #     }
  #   }
  # }
  #
  # p <- p +
  #   geom_hline(yintercept = 0, colour = "red", lty = 2) +
  #   geom_line(aes(y = `50%`))

  return(p)
}


#' \code{\link[ggplot2]{ggplot}} plotting function for Bayesian VAR predictions
#'
#' Plotting function for predictions obtained from \code{\link{predict.bvar}}
#' using \code{\link[ggplot2]{ggplot}}.
#' Predictions of all or a subset of the available variables can be
#' plotted.
#'
#' @param x A \code{bvar_fcast} object, obtained from \code{\link{predict.bvar}}
#' or a \code{bvar} object, obtained from \code{\link{bvar}}. In the latter
#' case predictions will be calculated ex-post using arguments from the ellipsis.
#' @param vars Optional numeric or character vector. Used to subset the plot to
#' certain variables by position or name (must be available). Defaults to
#' \code{NULL}, i.e. all variables.
#' @param col Character vector. Colour(s) of the lines delineating credible
#' intervals. Single values will be recycled if necessary. Recycled HEX color
#' codes are varied in transparency if not provided (e.g. "#737373FF").
#' @param t_back Integer scalar. Number of observed datapoints to plot ahead of
#' the forecast.
#' @param area Logical scalar. Whether to fill the credible intervals.
#' @param fill Character vector. Colour(s) to fill the credible intervals with.
#' See \emph{col} for more information.
#' @param variables Optional character vector. Names of all variables in the
#' object. Used to subset and title. Taken from \code{x$variables} if available.
#' @param orientation String indicating the orientation of the plots. Defaults
#' to \code{"vertical"}; may be set to \code{"horizontal"}.
#' @param ... Other parameters for calculating predictions if \code{x}
#' is not a \code{bvar_fcast} object.
#'
#' @return Returns an object of class \code{\link[ggplot2]{ggplot}}.
#'
#' @seealso \code{\link{bvar}}; \code{\link{irf.bvar}}
#'
#' @keywords BVAR fcast analysis ggplot
#'
#' @export
fcast_plot_gg <- function(
  x,
  vars = NULL,
  col = "#737373",
  t_back = 16,
  area = FALSE,
  fill = "#808080",
  variables = NULL,
  orientation = c("vertical", "horizontal"),
  ...) {

  if(!inherits(x, "bvar") && !inherits(x, "bvar_fcast")) {
    stop("Please provide a `bvar` or `bvar_fcast` object.")
  }

  if(inherits(x, "bvar")) {
    if(is.null(x[["fcast"]])) {message("No forecasts found. Calculating...")}
    x <- predict(x, ...)
  }

  orientation <- match.arg(orientation)

  t_back <- int_check(t_back, min = 0, max = nrow(x[["data"]]),
                      msg = "t_back too large.")

  df <- gg_df_fcast(x, t_back)
  M <- length(levels(df[["vars"]]))
  P_n <- levels(df[["quant"]])
  P <- length(P_n)
  df <- tidyr::pivot_wider(df, names_from = quant)

  variables <- name_deps(variables = if(is.null(variables)) {
    x[["variables"]]} else {variables}, M = M)

  pos <- pos_vars(vars, variables, M)

  df <- df[which(df[["vars"]] %in% variables[pos]), ]

  df[["vars"]] <- factor(paste("Forecast of:", df[["vars"]]))

  p <- ggplot(df, aes(x = time)) +
    scale_x_continuous(breaks = function(x) {
      unique(floor(pretty(seq((-t_back + 1), (max(x)) * 1.1))))
    }) +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(strip.background = element_blank()) +
    facet_wrap(.~vars, scales = "free_y",
               dir = ifelse(orientation == "vertical", "v", "h"))

  if(P > 1) {
    if(area) {
      fill <- fill_ci_col(x = integer(), y = fill, P = P)
      for(pp in seq_len(P - 1)) {
        p <- p + geom_ribbon(aes(ymin = .data[[P_n[pp]]],
                                 ymax = .data[[P_n[pp + 1]]]),
                             fill = fill[pp])
      }
    } else {
      col <- fill_ci_col(x = integer(), y = col, P = P)
      names(col) <- P_n[which(P_n != "50%")]
      for(pp in P_n[which(P_n != "50%")]) {
        p <- p + geom_line(aes(y = .data[[pp]]),
                           color = col[pp], na.rm = TRUE)
      }
    }
  }
  p <- p +
    geom_vline(xintercept = 1, colour = "red", lty = 2) +
    geom_line(aes(y = `50%`))

  return(p)
}


#' \code{\link[ggplot2]{ggplot}} plotting function for Bayesian VARs
#'
#' Function to plot trace and densities of coefficient, hyperparameter and
#' marginalised draws obtained from \code{\link{bvar}}. Several types of plot
#' are available via the argument \emph{type}, including traces, densities,
#' plots of forecasts and impulse responses.
#'
#' @param x A \code{bvar} object, obtained from \code{\link{bvar}}.
#' @param type A string with the type of plot desired. The default option
#' \code{"full"} plots both densities and traces.
#' @param vars Character vector used to select variables. Elements are matched
#' to hyperparameters or coefficients. Coefficients may be matched based on
#' the dependent variable (by providing the name or position) or the
#' explanatory variables (by providing the name and the desired lag). See the
#' example section for a demonstration. Defaults to \code{NULL}, i.e. all
#' hyperparameters.
#' @param vars_impulse,vars_response Optional character or integer vectors used
#' to select coefficents. Dependent variables are specified with
#' \emph{vars_response}, explanatory ones with \emph{vars_impulse}. Default to
#' \code{NULL}, indicating that no coefficients will be displayed.
#' @param quants Optional numeric vector for displaying quantiles of parameter
#' draws.
#' @param orientation String indicating the orientation of the plots. Defaults
#' to \code{"vertical"}; may be set to \code{"horizontal"}.
#' @param chains List of \code{bvar} objects. Contents are then added to trace
#' and density plots to help assessing covergence.
#' @param ... Other parameters for calculating predictions or impulse responses
#' if \code{type} is \code{"irf"} or \code{"fcast"}, respectively.
#'
#' @return Returns an object of class \code{\link[ggplot2]{ggplot}}.
#'
#' @seealso \code{\link{bvar}}.
#'
#' @keywords BVAR MCMC ggplot analysis hyperparameters
#'
#' @export
bvar_plot_gg <- function(
  x,
  type = c("full", "trace", "density", "irf", "fcast"),
  vars = NULL,
  vars_response = NULL,
  vars_impulse = NULL,
  quants = NULL,
  orientation = c("vertical", "horizontal"),
  chains = list(),
  ...) {

  if(!inherits(x, "bvar")) {
    if(inherits(x[[1]], "bvar")) { # Allow chains to x
      chains <- x
      x <- x[[1]]
      chains[[1]] <- NULL
    } else {stop("Please provide a `bvar` object.")}
  }

  type <- match.arg(type)
  orientation <- match.arg(orientation)

  # Forward and return if "irf" or "fcast"
  if(type == "irf") {
    if(is.null(x[["irf"]])) {message("No `bvar_irf` found. Calculating...")}
    return(gg_plot.bvar_irf(
      irf(x), vars_response = vars_response, vars_impulse = vars_impulse,
      variables = x[["variables"]], ...))
  }
  if(type == "fcast") {
    if(is.null(x[["fcast"]])) {message("No `bvar_fcast` found. Calculating...")}
    return(gg_plot.bvar_fcast(
      predict(x), vars = vars, variables = x[["variables"]], ...))
  }

  if(inherits(chains, "bvar")) {chains <- list(chains)}
  lapply(chains, function(x) {if(!inherits(x, "bvar")) {
    stop("Please provide `bvar` objects to the chains parameter.")
  }})

  prep <- prep_data(x, vars = vars, vars_response = vars_response,
                    vars_impulse = vars_impulse,
                    chains = chains, check_chains = FALSE)

  df <- gg_df_bvar(prep)

  if(!is.null(quants)) {
    quants <- quantile_check(quants)
  }

  plot_list <- list()

  if(type != "density") {
    p_trace <- ggplot(df, aes(y = value, x = x_axis)) +
      labs(x = NULL, y = NULL, color = "MCMC chain") +
      theme_bw() +
      theme(strip.background = element_blank()) +
      facet_wrap(.~vars, scales = "free",
                 ncol = ifelse(orientation == "vertical",
                               1,
                               length(unique(df[["vars"]])))) +
      theme(legend.position = "bottom")

    if(length(prep[["chains"]]) == 0) {
      p_trace <- p_trace + geom_line()
    } else {
      p_trace <- p_trace + geom_line(aes(color = chain), alpha = 0.5)
    }

    if(length(quants) != 0) {
      for(i in quants) {
        df_quants <- df %>% group_by(chain, vars) %>%
          summarise(quant = quantile(value, i))
        if(length(prep[["chains"]]) == 0) {
          p_trace <- p_trace + geom_hline(data = df_quants,
                                          aes(yintercept = quant),
                                          color = "red", lty = 2)
        } else{
          p_trace <- p_trace + geom_hline(data = df_quants,
                                          aes(yintercept = quant,
                                              color = chain), lty = 2)
        }
      }
    }

    if(type == "full") {
      p_trace <- p_trace + theme(legend.position = "none")
    }

    plot_list[["trace"]] <- p_trace
  }

  if(type != "trace") {

    p_dens <- ggplot(df, aes(x = value)) +
      labs(x = NULL, y = NULL, fill = "MCMC chain", color = "MCMC chain") +
      theme_bw() +
      theme(strip.background = element_blank()) +
      facet_wrap(.~vars, scales = "free",
                 ncol = ifelse(orientation == "vertical",
                               1,
                               length(unique(df[["vars"]])))) +
      theme(legend.position = "bottom")

    if(length(prep[["chains"]]) == 0) {
      p_dens <- p_dens +
        geom_density(fill = "grey", alpha = 0.5)
    } else {
      p_dens <- p_dens +
        geom_density(aes(color = chain, fill = chain), alpha = 0.5)
    }

    if(length(quants) != 0) {
      for(i in quants) {
        df_quants <- df %>% group_by(chain, vars) %>%
          summarise(quant = quantile(value, i))
        if(length(prep[["chains"]]) == 0) {
          p_dens <- p_dens +
            geom_vline(data = df_quants, aes(xintercept = quant),
                       color = "red", lty = 2)
        } else{
          p_dens <- p_dens +
            geom_vline(data = df_quants,
                       aes(xintercept = quant, color = chain), lty = 2)
        }
      }
    }

    if(type == "full") {
      p_legend <- cowplot::get_legend(p_dens)
      p_dens <- p_dens + theme(legend.position = "none")
    }

    plot_list[["dens"]] <- p_dens
  }

  if(type == "full") {
    return(cowplot::plot_grid(
      cowplot::plot_grid(plotlist = plot_list,
                         ncol = ifelse(orientation == "vertical", 2, 1),
                         align = "v", axis = "b"),
      p_legend, ncol = 1, rel_heights = c(10, 1)))
  } else {
    return(plot_list[[1]])
  }
}
