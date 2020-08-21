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

tidy.bvar <- function(
  x,
  vars = NULL,
  vars_response = NULL,
  vars_impulse = NULL,
  variables = NULL,
  quants = NULL,
  chains = list(),
  ...) {

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

  return(as_tibble(out))
}
