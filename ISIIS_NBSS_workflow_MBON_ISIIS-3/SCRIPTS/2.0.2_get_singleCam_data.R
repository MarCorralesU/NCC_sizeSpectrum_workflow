library(dplyr)        # bind_rows, arrange, mutate, filter
library(purrr)        # map, map_dfr, keep
library(stringr)
library(fs)
library(optparse)
library(lubridate)
library(data.table)

base_dir <- getwd()

option_list <- list(
  make_option(c("-c","--seasonYear"), type = "character",
              help = "Cruise id: w22/s22/w23/s23", metavar = "character"),
  make_option(c("-t","--Transect_ID"), type = "character",
              help = "Transect id: GH/OB/CR/CM/NH/HH/RR", metavar = "character"),
  make_option(c("-m", "--chunk_minutes"), type = "integer", default = 30,
            help = "Max minutes per output interval file (default: 30).")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$seasonYear) || is.null(opt$Transect_ID)) {
  stop("Error: Both --seasonYear and --Transect_ID are required.")
}

seasonYear  <- opt$seasonYear
Transect_ID <- opt$Transect_ID
chunk_minutes <- as.integer(opt$chunk_minutes)

message("▶ seasonYear = ", seasonYear, " | Transect_ID = ", Transect_ID)

# ---- Paths and filters ----
combined_dir <- file.path(base_dir, "data", "combined_data")
out_dir <- file.path(base_dir, "data", "combined_data_2Cams", seasonYear, Transect_ID)
if (!dir_exists(out_dir)) {
  dir_create(out_dir, recurse = TRUE)
  message("Created: ", out_dir)
}

# Helper: normalize camera labels across fields
clean_camera_labels <- function(df) {
  if ("tar_id" %in% names(df)) {
    df$tar_id <- str_replace_all(df$tar_id, "Camera_Top", "CamTop")
    df$tar_id <- str_replace_all(df$tar_id, "Camera_Bottom", "CamBottom")
  }

  if ("image" %in% names(df)) {
    df$image <- str_replace_all(df$image, "Camera_Top", "CamTop")
    df$image <- str_replace_all(df$image, "Camera_Bottom", "CamBottom")
  }

  if ("img" %in% names(df)) {
    df$img <- str_replace_all(df$img, "Camera_Top", "CamTop")
    df$img <- str_replace_all(df$img, "Camera_Bottom", "CamBottom")
  }

  df
}
#      x
#
#
#
#
#  df

# Read just the time column, fix quirks, return min/max
summarize_times <- function(file_path) {
  dt <- fread(file_path, select = "time",
              colClasses = list(character = "time"),
              showProgress = FALSE)

  # Fix "0022" year glitch if present
  if (any(grepl("^0022", dt$time))) dt$time <- gsub("^0022", "2022", dt$time)

  # Strip timezone offset at end (e.g., +00:00 / -07:00)
  dt$time <- sub("([-+]\\d{2}:\\d{2})$", "", dt$time)

  # Parse as Pacific (to match your pipeline)
  tt <- as.POSIXct(dt$time, format = "%Y-%m-%d %H:%M:%OS",
                   tz = "America/Los_Angeles")

  tibble(
    source_file = file_path,
    start = min(tt, na.rm = TRUE),
    end   = max(tt, na.rm = TRUE)
  )
}

# Discover candidate files
all_csvs <- list.files(combined_dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)

# Must include both cruise & transect substrings anywhere in the path
filtered <- keep(all_csvs, ~ all(str_detect(., c(seasonYear, Transect_ID))))

# Also ensure parent dir starts with "<cruise>_<transect>"
valid_parent <- paste0("^", seasonYear, "_", Transect_ID, ".*")
filtered <- keep(filtered, ~ str_detect(basename(dirname(.)), valid_parent))

camtop_files     <- filtered[str_detect(filtered, "CamTop")]
cambottom_files  <- filtered[str_detect(filtered, "CamBottom(?![[:alnum:]])")]  # exactly CamBottom
cambottom2_files <- filtered[str_detect(filtered, "CamBottom2")]

if (length(camtop_files) == 0 && length(c(cambottom_files, cambottom2_files)) == 0) {
  stop("No matching files found for ", seasonYear, " ", Transect_ID, " in ", combined_dir)
}

message(
  "Found: ",
  length(camtop_files), " CamTop | ",
  length(cambottom_files), " CamBottom | ",
  length(cambottom2_files), " CamBottom2"
)

# Summaries
camtop_times     <- if (length(camtop_files))     map_dfr(camtop_files, summarize_times)     else tibble()
cambottom_times  <- if (length(cambottom_files))  map_dfr(cambottom_files, summarize_times)  else tibble()
cambottom2_times <- if (length(cambottom2_files)) map_dfr(cambottom2_files, summarize_times) else tibble()

# Combine Bottom + Bottom2 into one group (BottomAll)
bottom_all_times <- bind_rows(cambottom_times, cambottom2_times)

# ---- Interval helpers (use numeric seconds for simple math) ----
to_num <- function(x) as.numeric(x)  # POSIXct -> seconds

merge_intervals <- function(df) {
  # df: columns start, end (POSIXct)
  if (nrow(df) == 0) return(df)
  x <- df %>% arrange(start) %>% mutate(s = to_num(start), e = to_num(end))
  out <- list()
  cur_s <- x$s[1]; cur_e <- x$e[1]
  for (i in seq_len(nrow(x))[-1]) {
    if (x$s[i] <= cur_e) {
      cur_e <- max(cur_e, x$e[i])
    } else {
      out[[length(out) + 1]] <- c(cur_s, cur_e)
      cur_s <- x$s[i]; cur_e <- x$e[i]
    }
  }
  out[[length(out) + 1]] <- c(cur_s, cur_e)
  as_tibble(do.call(rbind, out)) %>%
    rename(s = V1, e = V2) %>%
    mutate(
      start = as.POSIXct(s, origin = "1970-01-01", tz = "America/Los_Angeles"),
      end   = as.POSIXct(e, origin = "1970-01-01", tz = "America/Los_Angeles")
    ) %>%
    select(start, end)
}

subtract_union <- function(A_merged, B_merged) {
  # Return portions of A not covered by the union B
  if (nrow(A_merged) == 0) return(A_merged[0, ])
  if (nrow(B_merged) == 0) return(A_merged)

  res <- list()
  for (k in seq_len(nrow(A_merged))) {
    a_s <- to_num(A_merged$start[k]); a_e <- to_num(A_merged$end[k])
    cur <- a_s
    for (j in seq_len(nrow(B_merged))) {
      b_s <- to_num(B_merged$start[j]); b_e <- to_num(B_merged$end[j])
      if (b_e <= cur || b_s >= a_e) next  # no overlap
      if (b_s > cur) {
        res[[length(res) + 1]] <- c(cur, min(b_s, a_e))
      }
      cur <- max(cur, b_e)
      if (cur >= a_e) break
    }
    if (cur < a_e) res[[length(res) + 1]] <- c(cur, a_e)
  }
  if (!length(res)) return(A_merged[0, ])
  as_tibble(do.call(rbind, res)) %>%
    rename(s = V1, e = V2) %>%
    mutate(
      start = as.POSIXct(s, origin = "1970-01-01", tz = "America/Los_Angeles"),
      end   = as.POSIXct(e, origin = "1970-01-01", tz = "America/Los_Angeles")
    ) %>%
    select(start, end) %>%
    arrange(start)
}

# ---- Build unions ----
top_union    <- merge_intervals(select(camtop_times, start, end))
bottom_union <- merge_intervals(select(bottom_all_times, start, end))

# Exclusive coverage
top_only    <- subtract_union(top_union, bottom_union)
bottom_only <- subtract_union(bottom_union, top_union)

message("Intervals → Top-only: ", nrow(top_only), " | Bottom-only: ", nrow(bottom_only))

# ---- Read rows within an interval, from all relevant files in a camera group ----
read_rows_for_interval <- function(files, interval) {
  # Use file summaries to prefilter by [min,max] overlap before reading fully
  if (!length(files)) return(tibble())

  # Build a lookup of min/max for quick prefilter
  sums <- map_dfr(files, summarize_times)
  sums <- sums %>% filter(!(end < interval$start | start > interval$end))
  if (!nrow(sums)) return(tibble())

  parts <- lapply(sums$source_file, function(fp) {
    dt <- fread(fp, colClasses = list(character = "time"), showProgress = FALSE)

    if (any(grepl("^0022", dt$time))) dt$time <- gsub("^0022", "2022", dt$time)
    dt$time <- sub("([-+]\\d{2}:\\d{2})$", "", dt$time)
    dt$time <- as.POSIXct(dt$time, format = "%Y-%m-%d %H:%M:%OS",
                          tz = "America/Los_Angeles")

    dt <- as_tibble(dt) %>%
      filter(time >= interval$start & time <= interval$end)

    if (nrow(dt)) clean_camera_labels(dt) else dt
  })

  bind_rows(parts) %>% arrange(time)
}

# Split any interval longer than 'chunk_minutes' into contiguous sub-intervals.
# To avoid duplicate rows on boundaries, all but the last chunk end at (next_break - 1 ms).
chunk_large_intervals <- function(intervals, chunk_minutes) {
  if (!nrow(intervals)) return(intervals)
  out <- vector("list", 0L)

  for (i in seq_len(nrow(intervals))) {
    s <- intervals$start[i]
    e <- intervals$end[i]
    dur_min <- as.numeric(difftime(e, s, units = "mins"))

    if (dur_min <= chunk_minutes) {
      out[[length(out) + 1L]] <- tibble::tibble(start = s, end = e)
      next
    }

    # Build breaks every 'chunk_minutes'
    breaks <- seq(s, e, by = paste0(chunk_minutes, " mins"))
    if (tail(breaks, 1) < e) breaks <- c(breaks, e)

    for (j in seq_len(length(breaks) - 1L)) {
      s1 <- breaks[j]
      e1 <- breaks[j + 1L]
      # Make intermediate chunks right-open by backing off 1 ms to avoid boundary dupes
      if (j < length(breaks) - 1L) {
        e1 <- e1 - lubridate::milliseconds(1)
      }
      out[[length(out) + 1L]] <- tibble::tibble(start = s1, end = e1)
    }
  }

  dplyr::bind_rows(out) %>% dplyr::arrange(start)
}

# Bind two data frames/tables, harmonize columns, drop exact dup rows, sort by time if present
harmonize_bind_unique <- function(old_df, new_df) {
  old_dt <- as.data.table(old_df)
  new_dt <- as.data.table(new_df)

  # Ensure same columns
  cols <- union(names(old_dt), names(new_dt))
  for (nm in setdiff(cols, names(old_dt))) old_dt[, (nm) := NA]
  for (nm in setdiff(cols, names(new_dt))) new_dt[, (nm) := NA]
  setcolorder(old_dt, cols)
  setcolorder(new_dt, cols)

  # Character columns that should stay character
  char_cols <- intersect(c(
    "img", "image", "tar_id", "predicted_taxon", "group",
    "Cruise_name", "Path"
  ), cols)

  for (nm in char_cols) {
    old_dt[, (nm) := as.character(get(nm))]
    new_dt[, (nm) := as.character(get(nm))]
  }

  # Force time to POSIXct in both objects
  if ("time" %in% cols) {
    old_dt[, time := as.POSIXct(time, format = "%Y-%m-%d %H:%M:%OS", tz = "America/Los_Angeles")]
    new_dt[, time := as.POSIXct(time, format = "%Y-%m-%d %H:%M:%OS", tz = "America/Los_Angeles")]
  }

  merged <- unique(rbindlist(list(old_dt, new_dt), use.names = TRUE, fill = TRUE))

  if ("time" %in% names(merged)) {
    setorder(merged, time)
  }

  merged
}

#
#  # Ensure same columns (fill missing with NA) and same order
#
#  merged


# ---- Emit files for exclusives ----
write_exclusive <- function(intervals, group_name, src_files) {
  if (!nrow(intervals)) return(invisible(NULL))

  # Split long intervals into smaller windows (default: chunk_minutes from your CLI)
  intervals <- chunk_large_intervals(intervals, chunk_minutes)

  for (i in seq_len(nrow(intervals))) {
    iv <- intervals[i, ]
    message(sprintf("⏱  %s-only %s → %s", group_name, iv$start, iv$end))

    rows <- read_rows_for_interval(src_files, iv)

    bad_rows <- sum(is.na(rows$img))

    if (bad_rows > 0) {
      message("Removed ", bad_rows, " rows with missing img in interval ",
              iv$start, " → ", iv$end)
    }

    rows <- rows %>% dplyr::filter(!is.na(img))

    if (!nrow(rows)) {
      message("   (no rows found in files; skipping)")
      next
    }

    start_str <- gsub(" ", "_", as.character(iv$start), fixed = TRUE)
    end_str   <- gsub(" ", "_", as.character(iv$end),   fixed = TRUE)
    out_file  <- file.path(
      out_dir,
      paste0(seasonYear, "_", Transect_ID, "_joined_", start_str, "_to_", end_str, ".csv")
    )
    tmp_file <- paste0(out_file, ".tmp")

    if (fs::file_exists(out_file)) {
      # Merge new rows with what’s already there (drop exact duplicates)
      old <- data.table::fread(out_file, showProgress = FALSE)
      merged <- harmonize_bind_unique(old, rows)

      data.table::fwrite(merged, tmp_file)
      fs::file_move(tmp_file, out_file)
      message("🔄 Merged into: ", out_file, "  (", nrow(merged), " rows)")
    } else {
      data.table::fwrite(rows, out_file)
      message("✅ Wrote: ", out_file, "  (", nrow(rows), " rows)")
    }
  }
}

# CamTop-only
write_exclusive(top_only, "CamTop", camtop_files)

# Bottom-only (CamBottom + CamBottom2)
write_exclusive(bottom_only, "CamBottomAll", c(cambottom_files, cambottom2_files))

message("Done.")

