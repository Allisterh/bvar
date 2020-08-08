files_source <- list.files(path = "./R", pattern = ".R")
sapply(paste0("./R/", files_source), source, echo = FALSE)

gg_df_irf <- function(x) {

  stopifnot(inherits(x, "bvar_irf"))

  has_quants <- length(dim(x[["quants"]])) == 4L
  if(!has_quants) {
    temp <- array(NA, c(1, dim(x[["quants"]])))
    temp[1, , , ] <- x[["quants"]]
    dimnames(temp)[[1]] <- list("50%")
    x[["quants"]] <- temp
  }

  # So we know what's what
  dimnames(x[["quants"]])[[2]] <- x[["variables"]]
  dimnames(x[["quants"]])[[4]] <- x[["variables"]]
  dimnames(x[["quants"]])[[3]] <- seq(dim(x[["quants"]])[3])

  out <- as.data.frame.table(x[["quants"]]) # Magic base R
  colnames(out) <- c("quant", "response", "time", "impulse", "value")
  out[["time"]] <- as.integer(out[["time"]]) # Can't be a factor for the x axis

  return(out)
}


gg_df_fcast <- function(x, t_back) {

  stopifnot(inherits(x, "bvar_fcast"))

  # So we know what's what
  has_quants <- length(dim(x[["quants"]])) == 3L
  if(!has_quants) {
    temp <- array(NA, c(1, dim(x[["quants"]])))
    temp[1, , ] <- x[["quants"]]
    dimnames(temp)[[1]] <- list("50%")
    x[["quants"]] <- temp
  }

  time_x <- seq(-t_back + 1, dim(x[["quants"]])[2])
  out <- array(NA, c(dim(x[["quants"]])[1],
                         length(time_x),
                         dim(x[["quants"]])[3]))
  rownames(out) <- rownames(x[["quants"]])
  dimnames(out)[[2]] <- time_x
  dimnames(out)[[3]] <- x[["variables"]]
  out["50%", 1:t_back, ] <- x[["data"]][(nrow(x[["data"]]) - t_back + 1):nrow(x[["data"]]),]
  out[ , (t_back + 1):length(time_x), ] <- x[["quants"]]

  out <- as.data.frame.table(out) # Magic base R
  colnames(out) <- c("quant", "time", "vars", "value")
  out[["time"]] <- as.numeric(as.character(out[["time"]])) # Can't be a factor for the x axis

  return(out)
}


irf_plot_gg <- function(
  x,
  vars_response = NULL,
  vars_impulse = NULL,
  variables = NULL,
  area = FALSE,
  col = "#737373",
  fill = "#808080",
  scales = "fixed") {

  if(!inherits(x, "bvar") && !inherits(x, "bvar_irf")) {
    stop("Please provide a `bvar` or `bvar_irf` object.")
  }

  if(inherits(x, "bvar")) {
    if(is.null(x[["irf"]])) {message("No `bvar_irf` found. Calculating...")}
    x <- irf(x)
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
      for(pp in P_n[which(P_n != "50%")]) {
        p <- p + geom_line(aes(y = .data[[pp]]),
                           colour = col[which(P_n == pp)])
      }
    }
  }

  p <- p +
    geom_hline(yintercept = 0, colour = "red", lty = 2) +
    geom_line(aes(y = `50%`))

  return(p)
}



fcast_plot_gg <- function(
  x,
  vars = NULL,
  variables = NULL,
  orientation = "horizontal",
  t_back = 16,
  area = FALSE,
  col = "#737373",
  fill = "#808080"){
  if(!inherits(x, "bvar") && !inherits(x, "bvar_fcast")) {
    stop("Please provide a `bvar` or `bvar_fcast` object.")
  }
  if(inherits(x, "bvar")) {
    if(is.null(x[["fcast"]])) {message("No `bvar_fcast` found. Calculating...")}
    x <- predict(x)
  }

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
               dir = ifelse(orientation == "horizontal", "h", "v"))

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
      for(pp in P_n[which(P_n != "50%")]) {
        p <- p + geom_line(aes(y = .data[[pp]]),
                           colour = col[which(P_n == pp)])
      }
    }
  }
  p <- p +
    geom_vline(xintercept = 1, colour = "red", lty = 2) +
    geom_line(aes(y = `50%`))

  return(p)
}

bvar_plot_gg <- function(x, type, vars, vars_response, vars_impulse, chains, ...) {


}


#####
# Testing

data <- matrix(rnorm(3000, 1, 1), 1000, 3)
colnames(data) <- c("gdp", "inf", "irst")

#library(BVAR)
library(tidyverse)
library(ggplot2)

y <- bvar(data, lags = 5)
x <- irf(y, conf_bands = c(0.05, 0.16))
z <- irf(y, conf_bands = 0.5)


irf_plot_gg(x, area = TRUE, vars_response = c(1,2), vars_impulse = c("irst"))

irf_plot_gg(y, area = TRUE, vars_response = c(2,3), vars_impulse = c("gdp", "irst"), scales = "free")

irf_plot_gg(z, area = TRUE, vars_response = c(1,2), vars_impulse = c(2,3))


xx <- predict(y, conf_bands = c(0.05, 0.16))
zz <- predict(y, conf_bands = 0.5)

fcast_plot_gg(y, area = TRUE, vars = c("irst"))

fcast_plot_gg(xx, area = TRUE, vars = c("gdp", "irst"))

fcast_plot_gg(zz, area = TRUE, vars = c(1,2))
