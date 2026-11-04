library(data.table)
library(dplyr)
library(magrittr)
library(tidyr)
library(NLP)
library(here)
library(lubridate)
library(runner)
library(suncalc)
library(purrr)
library(slider)
library(rlang)


collapse_taxa_cols <- function(
  df,
  start_name = "appendicularian",
  end_name   = "tunicate_salp",
  ignore_case = FALSE,
  detritus_handling = c("keep", "remove")  # <-- NEW
) {
  detritus_handling <- match.arg(detritus_handling)

  nm <- names(df)
  if (ignore_case) {
    nm_l <- tolower(nm)
    s <- match(tolower(start_name), nm_l)
    e <- match(tolower(end_name),   nm_l)
  } else {
    s <- match(start_name, nm)
    e <- match(end_name,   nm)
  }
  if (is.na(s) || is.na(e) || s > e) {
    stop(
      "Couldn't locate a valid taxa block between '",
      start_name, "' and '", end_name, "'."
    )
  }

  # Columns to collapse (between anchors)
  data_cols <- nm[s:e]

  # Define category patterns
  patterns <- list(
    detritus             = "detritus|art[ie]fact",  # allows both
    protist              = "protist",
    phytoplankton        = "phytoplankton|diatom",
    copepod              = "copepod",
    crustacean           = "crustacean",
    predatory_gelatinous = "(hydromedusae|chaetognath|siphonophore)",
    grazer_gelatinous    = "(appendicularian|tunicate|pyrosome)",
    fish                 = "fish"
  )

  # Logical match vectors
  idx <- lapply(patterns, function(p) grepl(p, data_cols, ignore.case = TRUE))

  # Taxa buckets
  detritus_cols <- data_cols[idx$detritus]

  matched_non_det <- Reduce(
    `|`,
    idx[c(
      "protist","phytoplankton","copepod","crustacean",
      "predatory_gelatinous","grazer_gelatinous","fish"
    )]
  )

  matched_cols <- data_cols[matched_non_det]

  # Any remaining go to "other", NEVER detritus
  other_cols <- setdiff(data_cols, c(matched_cols, detritus_cols))

  # Safe row sums
  rs <- function(cols) {
    if (!length(cols)) return(rep(0, nrow(df)))
    m <- lapply(df[, cols, drop = FALSE],
                function(z) suppressWarnings(as.numeric(z)))
    rowSums(as.data.frame(m), na.rm = TRUE)
  }

  # Build grouped columns
  grouped <- tibble::tibble(
    protist              = rs(data_cols[idx$protist]),
    phytoplankton        = rs(data_cols[idx$phytoplankton]),
    copepod              = rs(data_cols[idx$copepod]),
    crustacean           = rs(data_cols[idx$crustacean]),
    predatory_gelatinous = rs(data_cols[idx$predatory_gelatinous]),
    grazer_gelatinous    = rs(data_cols[idx$grazer_gelatinous]),
    fish                 = rs(data_cols[idx$fish]),
    other                = rs(other_cols)
  )

  # LEFT & RIGHT dataframe segments
  keep_left  <- if (s > 1)        df[, 1:(s-1), drop = FALSE]        else df[, 0, drop = FALSE]
  keep_right <- if (e < ncol(df)) df[, (e+1):ncol(df), drop = FALSE] else df[, 0, drop = FALSE]

  # detritus output option
  if (detritus_handling == "keep" && length(detritus_cols)) {
    keep_det <- df[, detritus_cols, drop = FALSE]
  } else {
    keep_det <- NULL
  }

  dplyr::bind_cols(keep_left, grouped, keep_det, keep_right)
}
