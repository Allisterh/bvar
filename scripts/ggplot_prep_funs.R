#' Create dataframe for plotting impulse responses with \code{\link[ggplot2]{ggplot}}
#'
#' Helper function to create dataframe from a \code{\link{bvar_irf}} object
#' suitable for further processing with \code{\link[ggplot2]{ggplot}}.
#'
#'
#' @param x A \code{\link{bvar_irf}} object.
#'
#' @return Returns a dataframe with data in suitable format to be passed on
#' to \code{\link[ggplot2]{ggplot}}
#'
#' @noRd
gg_df_irf <- function(x) {

  has_quants <- length(dim(x[["quants"]])) == 4L
  if(!has_quants) {
    temp <- array(NA, c(1, dim(x[["quants"]])))
    temp[1, , , ] <- x[["quants"]]
    dimnames(temp)[[1]] <- list("50%")
    x[["quants"]] <- temp
  }

  dimnames(x[["quants"]])[[2]] <- x[["variables"]]
  dimnames(x[["quants"]])[[4]] <- x[["variables"]]
  dimnames(x[["quants"]])[[3]] <- seq(dim(x[["quants"]])[3])

  out <- as.data.frame.table(x[["quants"]])
  out <- out[, c(4, 2, 3, 5, 1)]
  colnames(out) <- c("impulse", "response", "time", "value", "quant")
  out[["time"]] <- as.integer(out[["time"]])
  out <- out[order(out[["impulse"]], out[["response"]], out[["time"]]), ]

  return(out)
}


#' Create dataframe for plotting predictions with \code{\link[ggplot2]{ggplot}}
#'
#' Helper function to create dataframe from a \code{\link{bvar_fcast}} object
#' suitable for further processing with \code{\link[ggplot2]{ggplot}}.
#'
#'
#' @param x A \code{\link{bvar_fcast} object.
#'
#' @return Returns a dataframe with data in suitable format to be passed on
#' to \code{\link[ggplot2]{ggplot}}
#'
#' @noRd
gg_df_fcast <- function(x, t_back) {

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
  if(t_back > 0) {
    out["50%", 1:t_back, ] <- tail(x[["data"]], t_back)
  }
  out[ , (t_back + 1):length(time_x), ] <- x[["quants"]]

  out <- as.data.frame.table(out)
  out <- out[, c(3, 2, 4, 1)]
  colnames(out) <- c("vars", "time", "value", "quant")
  out[["time"]] <- as.numeric(levels(out[["time"]]))[out[["time"]]]
  out <- out[order(out[["vars"]], out[["time"]]), ]

  return(out)
}


#' Create dataframe for plotting convergence plots with \code{\link[ggplot2]{ggplot}}
#'
#' Helper function to create dataframe from a \code{\link{bvar}} object
#' suitable for further processing with \code{\link[ggplot2]{ggplot}}.
#'
#'
#' @param x A \code{\link{bvar}} object.
#'
#' @return Returns a dataframe with data in suitable format to be passed on
#' to \code{\link[ggplot2]{ggplot}}
#'
#' @noRd
gg_df_bvar <- function(x) {

  out <- as.data.frame.table(x[["data"]])
  names(out) <- c("chain", "vars", "value")
  out[["chain"]] <- 1

  for(ll in seq_len(length(x[["chains"]]))) {
    temp <- as.data.frame.table(x[["chains"]][[ll]])
    names(temp) <- c("chain", "vars", "value")
    temp[["chain"]] <- ll + 1
    out <- rbind(out, temp)
  }
  out[["chain"]] <- factor(out[["chain"]])
  out[["x_axis"]] <- seq_len(nrow(x[["data"]]))

  return(out)
}
