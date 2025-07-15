#' Retrieve the peptide metadata table into DuckDB, forcing atomic types
#'
#' * Downloads the RDS once, sanitizes types (logical, character, numeric),
#'   and writes into a DuckDB cache on disk.
#' * Subsequent calls return a lazy tbl_dbi without loading into R memory.
#'
#' @param force_refresh Logical. If TRUE, re-downloads and rebuilds the cache.
#' @return A dplyr tbl_dbi pointing to the `peptide_meta` table.
#' @export
get_peptide_meta <- function(force_refresh = FALSE) {
  # check if dependencies installed --> maybe hard dep in the future?
  rlang::check_installed(c("duckdb", "DBI", "dplyr", "withr"))

  # 1. Prep cache dir & DuckDB connection
  cache_dir <- tools::R_user_dir("phiper", "cache")
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
  duckdb_file <- file.path(cache_dir, "phip_cache.duckdb")
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = duckdb_file)

  # ## needed to suppress the messages on this as the defer function always
  # ## prints out a message to the console
  # suppressMessages(
  #   withr::defer(DBI::dbDisconnect(con, shutdown = TRUE), parent.frame())
  # )

  # 2. fast path: already cached? --> return
  ## the user can also force the evaluation by force_refresh arg
  if (!force_refresh && DBI::dbExistsTable(con, "peptide_meta")) {
    peptides_tbl <- dplyr::tbl(con, "peptide_meta")
    attr(peptides_tbl, "duckdb_con") <- con
    return(peptides_tbl)
  }

  # 3. download raw RDS from the github repo
  url <- paste0(
    "https://raw.githubusercontent.com/Polymerase3/phiper/",
    "master/library-metadata/all_libraries_with_important_info.rds"
  )
  tmp <- tempfile(fileext = ".rds")
  sha <- "cbfe17929c953a730375075335a0a58643650b8dcd62ebb00b562b1a9a36285c"

  message("Downloading peptide metadata…")

  ## safe download (fallbacks if file changed, or if download does not succeed)
  .safe_download(url, tmp, sha)

  ## reading the raw RDS file --> it needs a lot of polishin (is prolly
  ## python generated, see attributes)
  raw_meta <- readRDS(tmp)
  unlink(tmp)

  # 4. capture rownames and clear them --> we want a separate column for the
  # peptide_id
  peptide_ids <- rownames(raw_meta) %||% rep(NA_character_, nrow(raw_meta))
  rownames(raw_meta) <- NULL

  # 5. sanitize each column
  clean_list <- lapply(raw_meta, function(col) {
    # drop all attributes
    attributes(col) <- NULL

    # unlist it if column is a list
    if (typeof(col) == "list") {
      col <- unlist(col)
    }

    # replace NaN with NA's, also the NaNs as characters
    col[is.nan(col)] <- NA
    col[col == "NaN"] <- NA

    # 1) numeric 0/1 -> logical (in the Carlos's metadata binary variables are
    # saved as doubles; here a fallback for other types as well)
    if (all(col %in% c(0, 1, NA)) || all(col %in% c("0", "1", NA))) {
      return(as.logical(col))
    }

    # 2) character "TRUE"/"FALSE" → logical
    if (is.character(col) &&
      all(tolower(col[!is.na(col)]) %in% c("true", "false", NA))) {
      return(as.logical(col))
    }

    # 3) character column that really holds numeric values (but not just 0/1)
    if (is.character(col)) {
      # remove NAs for testing
      non_na <- col[!is.na(col)]
      # detect pure numeric strings (optionally scientific)
      is_num_str <- grepl(
        "[0-9]+",
        non_na
      )

      # only proceed if all non-NA entries are numeric strings
      # and not all of them are "0" or "1"
      if (length(non_na) > 0 &&
        all(is_num_str) &&
        !all(non_na %in% c("0", "1"))
      ) {
        return(suppressWarnings(as.numeric(col)))
      }
    }

    # 4) otherwise leave atomic as it was (logical/integer/double/character)
    col
  })

  # 6. prepend peptide_id column
  clean_list <- c(list(peptide_id = peptide_ids), clean_list)

  meta_df <- data.frame(
    clean_list,
    stringsAsFactors = FALSE,
    check.names = TRUE
  )

  # 7. write into DuckDB
  if (DBI::dbExistsTable(con, "peptide_meta")) {
    DBI::dbRemoveTable(con, "peptide_meta")
  }
  message("Importing sanitized metadata into DuckDB cache…")
  DBI::dbWriteTable(con, "peptide_meta", meta_df, overwrite = TRUE)

  # 8. return lazy handle --> the whole dataframe in the memory takes ~ 1GB
  peptides_tbl <- dplyr::tbl(con, "peptide_meta")
  attr(peptides_tbl, "duckdb_con") <- con
  return(peptides_tbl)
}

#' @keywords internal
.safe_download <- function(url,
                           dest,
                           sha_expected = NULL) {
  ## setup
  fs::dir_create(fs::path_dir(dest))
  methods <- c("", "libcurl", "curl")
  ok <- FALSE

  ## perform the actual download with given method (or at least try)
  for (m in methods) {
    status <- tryCatch(
      utils::download.file(url, dest,
        mode   = "wb",
        quiet  = TRUE,
        method = if (nzchar(m)) m else getOption("download.file.method")
      ),
      error = function(e) e,
      warning = function(w) w
    )
    if (identical(status, 0L)) {
      ok <- TRUE
      break
    }
  }

  ## Print out error if download did not succeed
  if (!ok || !fs::file_exists(dest) || fs::file_info(dest)$size == 0) {
    chk::abort_chk(sprintf("Failed to download %s", url))
  }

  ## Compare the checksums for the whole file, generate a warning if the file
  ## changed
  if (!is.null(sha_expected)) {
    sha_actual <- tryCatch(
      {
        strsplit(system2("sha256sum", dest, stdout = TRUE), "\\s+")[[1]][1]
      },
      error = function(e) NA_character_
    )

    if (is.na(sha_actual) ||
      !identical(tolower(sha_actual), tolower(sha_expected))) {
      cli::cli_warn(
        c("!" = paste0("Checksum mismatch: expected {.val {sha_expected}},",
                       "got {.val {sha_actual %||% 'NA'}}"))
      )
    }
  }

  invisible(dest)
}
