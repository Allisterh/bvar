#' FRED-MD and FRED-QD: Databases for Macroeconomic Research
#'
#' FRED-MD and FRED-QD are large macroeconomic databases, containing monthly
#' and quarterly time series that are frequently used in the literature. They
#' are intended to facilitate the reproduction of empirical work and simplify
#' data related tasks.
#' Included datasets are provided as is - transformation codes are available
#' in \code{system.file("fred_trans.rds", package = "BVAR")}. These can be
#' applied automatically with \code{\link{fred_transform}}.
#'
#' The versions of FRED-MD and FRED-QD that are provided here are licensed
#' under a modified ODC-BY 1.0 license that can be found in the provided
#' \emph{LICENSE} file. The provided versions are subset to variables that are
#' either in public domain or for which we were given permission to use.
#' For further details see McCracken and Ng (2016) or
#' \url{https://research.stlouisfed.org/econ/mccracken/fred-databases/}.
#' We would like to thank Michael McCracken and Serena Ng, Adrienne Brennecke
#' and the Federal Reserve Bank of St. Louis for creating, updating and making
#' available the datasets and many of the contained time series. We also thank
#' all other owners of included time series that permitted their use.
#'
#' @docType data
#'
#' @format A \code{data.frame} object with dates as rownames.
#'
#' @seealso \code{\link{fred_transform}}, \code{\link{fred_download}}
#'
#' @keywords FRED macroeconomics
#'
#' @references
#'   McCracken, M. W. and Ng, S. (2016) FRED-MD: A Monthly Database for
#'   Macroeconomic Research. \emph{Journal of Business & Economic Statistics},
#'   \bold{34:4}, 574-589, \doi{10.1080/07350015.2015.1086655}.
#'   McCracken, M. W., & Ng, S. (2020). FRED-QD: A Quarterly Database for
#'   Macroeconomic Research \bold{w26872}. National Bureau of Economic Research.
#'
#' @source \url{https://research.stlouisfed.org/econ/mccracken/fred-databases/}
"fred_qd"

#' @rdname fred_qd
"fred_md"

#' Download a FRED-MD or FRED-QD vintage
#'
#' Downloads a vintage of the FRED-MD or FRED-QD database as a numeric
#' \code{data.frame}, matching the structure of \code{\link{fred_md}}
#' and \code{\link{fred_qd}}. See
#' \url{https://research.stlouisfed.org/econ/mccracken/fred-databases/}
#'
#' Metadata rows are removed. Unlike the packaged datasets, the result
#' is not restricted to series that BVAR may redistribute.
#'
#' @param type Character; either \code{"qd"} (quarterly) or \code{"md"}
#'   (monthly).
#' @param vintage Character; defaults to the latest release, or
#'   a specific vintage in \code{"YYYY-MM"} format. Only recent
#'   vintages may be available.
#' @param ... Passed on to \code{\link[utils]{read.csv}}.
#'
#' @return A numeric \code{data.frame} with \code{Date} rownames.
#'
#' @seealso \code{\link{fred_transform}}, \code{\link{fred_md}},
#'   \code{\link{fred_qd}}
#'
#' @examples
#' \dontrun{
#' fred_download("md")
#' fred_download("md", "2026-09")
#' fred_download("qd", "2026-08")
#' }
#'
#' @export
fred_download <- function(
  type = c("qd", "md"),
  vintage = format(Sys.Date() - 31L, "%Y-%m"),
  ...
) {
  type <- match.arg(type)
  if (!grepl("^[0-9]{4}-[0-9]{2}$", vintage)) {
    stop("Please provide a vintage in 'YYYY-MM' format, e.g. '2026-09'.")
  }
  link <- paste0(
    "https://www.stlouisfed.org/-/media/project/frbstl/stlouisfed/research/",
    "fred-md/",
    if (type == "qd") "quarterly/" else "monthly/",
    vintage,
    "-",
    type,
    ".csv"
  )
  # Temporary CSV to download
  file <- tempfile(fileext = ".csv")
  on.exit(unlink(file), add = TRUE)

  status <- tryCatch(
    suppressWarnings(
      utils::download.file(link, file, quiet = TRUE, mode = "wb")
    ),
    error = function(e) 1L
  )
  if (!identical(status, 0L)) {
    stop("Failed to download FRED data from:\n\t", link)
  }

  data <- utils::read.csv(file, ...)

  n_labels <- if (type == "qd") 2L else 1L
  if (
    !"sasdate" %in% names(data) ||
      ncol(data) < 2L ||
      nrow(data) <= n_labels
  ) {
    stop("Unexpected FRED CSV format.")
  }
  # Remove metadata and empty rows
  data <- data[-seq_len(n_labels), , drop = FALSE]
  data <- data[!is.na(data$sasdate) & nzchar(data$sasdate), , drop = FALSE]

  dates <- as.Date(data$sasdate, "%m/%d/%Y")
  if (anyNA(dates) || anyDuplicated(dates)) {
    stop("Invalid dates in FRED CSV.")
  }
  rownames(data) <- dates
  data$sasdate <- NULL

  return(data)
}
