# ==============================================================================
# phiper checks + additional helpers (ASCII-only, unified with phiperio logger)
# logging utilities (.ph_opt, .ph_warn, .ph_abort, etc.) live in phiperio
# ==============================================================================

# ---- Back-compatibility (vecmatch) -------------------------------------------
# drop-in replacement for older .msg() funct
.msg <- function(verbose, ...) {
  if (isTRUE(verbose)) .ph_log_info(headline = paste(...))
  invisible(NULL)
}

# original conditional helper to not break down older code
# upgraded to the unified phiperio style
.chk_cond <- function(condition,
                      error_message,
                      error = TRUE,
                      step = NULL,
                      bullets = NULL,
                      ...) {
  # log nothing
  if (!isTRUE(condition)) {
    return(invisible(FALSE))
  }

  # print error and abort exec
  if (isTRUE(error)) {
    .ph_abort(headline = error_message, step = step, bullets = bullets, ...)
  } else {
    # print warning and go on
    .ph_warn(headline = error_message, step = step, bullets = bullets, ...)
  }
  invisible(TRUE)
}

# ==============================================================================
# phiper checks (ASCII-only, unified with phiperio logger)
# depends on: .ph_abort(), .ph_warn(), .chk_cond(), word_list(), add_quotes()
# ==============================================================================

# -- check if filename has given extension ------------------------------------
# comes in handy when loading .csv or .parquet; provide filename and vector of
# extensions to check (eg c(".csv", ".parquet"))
.chk_extension <- function(name,
                           x_name,
                           ext_vec) {
  if (is.null(ext_vec) || !length(ext_vec)) {
    return(invisible(TRUE))
  }

  base <- basename(name %||% "") # extracting filename from paths
  parts <- strsplit(base, "\\.", fixed = FALSE)[[1]] # names + complex ext

  # taking last ext after . (eg .tar.gz --> .gz)
  ext <- if (length(parts) > 1L) {
    tolower(paste(parts[-1L], collapse = "."))
  } else {
    ""
  }

  norm <- function(x) sub("^\\.+", "", tolower(x)) # normalize ext
  got <- if (nzchar(ext)) norm(ext) else "<none>"
  ok <- nzchar(ext) && got %in% norm(ext_vec) # final ext check

  if (!ok) {
    .ph_abort(
      headline = sprintf("Invalid file extension for `%s`.", x_name),
      step = sprintf("validating path: %s", name),
      bullets = c(
        sprintf("got: %s", add_quotes(got, 2L)),
        sprintf(
          "allowed: %s",
          word_list(add_quotes(norm(ext_vec), 2L), and_or = "or")
        )
      )
    )
  }
  invisible(TRUE)
}

# -- check if NULL and replace with default when TRUE (warn in unified style) --
.chk_null_default <- function(x,
                              x_name,
                              method,
                              default) {
  if (is.null(x)) {
    # format the default for print
    fmt <- function(v) {
      if (is.character(v) && length(v) == 1L) {
        return(add_quotes(v, 2L))
      }
      if (is.atomic(v) && length(v) == 1L) {
        return(as.character(v))
      }
      sprintf("<%s>", paste(class(v), collapse = "/"))
    }

    # generate warning and the replace
    .ph_warn(
      headline = sprintf("Argument `%s` missing; using default.", x_name),
      step     = sprintf("method: %s", add_quotes(method, 2L)),
      bullets  = sprintf("default: %s", fmt(default))
    )
    x <- default
  }
  x
}

# -- validate path to file -----------------------------------------------------
.chk_path <- function(path,
                      arg_name,
                      extension,
                      is_dir = FALSE) {
  ## error when path not a string
  .chk_cond(
    !chk::vld_string(path),
    sprintf("`%s` must be a character scalar.", arg_name),
    step    = "path validation",
    bullets = sprintf("got class: %s", paste(class(path), collapse = "/"))
  )

  ## error when path does not exist
  if (is_dir) {
    .chk_cond(
      !chk::vld_dir(path),
      sprintf("Folder for `%s` does not exist.", arg_name),
      step    = "path validation",
      bullets = sprintf("path: %s", path)
    )
  } else {
    .chk_cond(
      !chk::vld_file(path),
      sprintf("File for `%s` does not exist.", arg_name),
      step    = "path validation",
      bullets = sprintf("path: %s", path)
    )
  }

  # optionally extension check if provided
  if (!missing(extension) && length(extension)) {
      ## error when both is_dir and extension are given
      .chk_cond(
        is_dir,
        sprintf("Can't check if `%s` is both a valid direcotry and has a certain extension", arg_name),
        step    = "path validation",
        bullets = sprintf("path: %s", path)
      )

    .chk_extension(
      path,
      arg_name,
      extension
    )
  }

  invisible(TRUE)
}

# -- clean wordlists for message generation ------------------------------------
# for multiple arguments/values
word_list <- function(word_list = NULL,
                      and_or = "and",
                      is_are = FALSE,
                      quotes = FALSE) {
  # Make "a and b" / "a, b, and c"; optionally append "is/are".
  word_list <- setdiff(word_list, c(NA_character_, ""))

  if (is.null(word_list)) {
    out <- ""
    attr(out, "plural") <- FALSE
    return(out)
  }

  word_list <- add_quotes(word_list, quotes)

  len_wl <- length(word_list)

  if (len_wl == 1L) {
    out <- word_list
    if (is_are) out <- paste(out, "is")
    attr(out, "plural") <- FALSE
    return(out)
  }

  if (is.null(and_or) || isFALSE(and_or)) {
    out <- paste(word_list, collapse = ", ")
  } else {
    and_or <- match.arg(and_or, c("and", "or"))
    if (len_wl == 2L) {
      out <- sprintf("%s %s %s", word_list[1L], and_or, word_list[2L])
    } else {
      out <- sprintf(
        "%s, %s %s",
        paste(word_list[-len_wl], collapse = ", "),
        and_or, word_list[len_wl]
      )
    }
  }

  if (is_are) out <- sprintf("%s are", out)
  attr(out, "plural") <- TRUE
  out
}

# -- quoting helper (unified error style) --------------------------------------
# define number of quotes you want --> for printing logs/messages/warnings
# or define the quotes itself as a string
add_quotes <- function(x,
                       quotes = 2L) {
  if (isFALSE(quotes)) {
    return(x)
  }
  if (isTRUE(quotes)) quotes <- '"'

  if (chk::vld_string(quotes)) {
    return(paste0(quotes, x, quotes))
  }

  if (!chk::vld_count(quotes) || quotes > 2) {
    .ph_abort(
      headline = "Invalid `quotes` argument.",
      step = "formatting add_quotes()",
      bullets = c(
        "allowed: FALSE, TRUE, 0, 1, 2, or a single-character string",
        sprintf("got class: %s", paste(class(quotes), collapse = "/"))
      )
    )
  }

  if (quotes == 0L) {
    return(x)
  }
  if (quotes == 1L) {
    return(sprintf("'%s'", x))
  }
  sprintf('"%s"', x)
}


#' @title Path to example PhIP-Seq datasets shipped with phiper
#'
#' @param name Character scalar. Name of the example dataset.
#'   Currently supported: `"phip_mixture"`.
#'
#' @return A character scalar with an absolute path to the file.
#'
#' @examples
#' sim_path <- phip_example_path("phip_mixture")
#' # phip_obj <- phip_convert(sim_path)
#'
#' @export
phip_example_path <- function(name = c("phip_mixture")) {
  name <- match.arg(name)
  fname <- switch(
    name,
    phip_mixture = "phip_mixture.parquet"
  )

  path <- system.file("extdata", fname, package = "phiper")
  if (path == "") {
    stop("File ", fname, " not found in extdata/", call. = FALSE)
  }
  path
}

#' @title Load Example PhIP-Seq Dataset as <phip_data>
#'
#' @description
#' Convenience helper to quickly load a shipped example dataset ("phip_mixture") into a `<phip_data>` object,
#' suitable for downstream analysis and visualization. This function wraps \code{\link{phip_convert}},
#' automatically supplying the correct parameters for the included example data.
#'
#' @param name Character scalar. Name of the shipped example dataset.
#'  Currently supported: \code{"phip_mixture"}, \code{"small_mixture"}.
#'
#' @return A `<phip_data>` object created from the chosen example dataset.
#'
#' @examples
#' # Load the example data shipped with the package:
#' ex <- phip_load_example_data()
#' # ex is now a <phip_data> object ready for analysis
#'
#' # Specify the dataset name explicitly
#' ex2 <- phip_load_example_data("small_mixture")
#'
#' # Use with plotting functions
#' p = plot_enrichment_counts(ex, group_cols = "timepoint")
#'
#' @export
phip_load_example_data <- local({
  cache_env <- new.env(parent = emptyenv())
  cache_env$loaded <- list()

  function(name = c("phip_mixture", "small_mixture")) {
    name <- match.arg(name)

    # Check if already in cache
    if (name %in% names(cache_env$loaded)) return(cache_env$loaded[[name]])

    if (name == "small_mixture") {
      ps <- phip_load_example_data(name = "phip_mixture")

      # small subset for speed: 5 peptides at time t1
      keep_pep <- c("16627", "5243", "24799", "16196", "18003")
      dat_cols <- dplyr::tbl_vars(ps$data_long)
      tp_col <- "timepoint"

      ps <- ps |>
        dplyr::filter(
          peptide_id %in% keep_pep,
          !!rlang::sym(tp_col) == "T1"
        ) |>
        dplyr::collect()

    } else {
      ps <- phip_convert(
        data_long_path = phip_example_path(name),
        peptide_library = TRUE,
        subject_id = "subject_id",
        peptide_id = "peptide_id",
        sample_id  = "sample_id",
        exist      = "exist",
        timepoint  = "time",
        fold_change= "fold_change",
        materialise_table = TRUE,
        auto_expand = FALSE,
        n_cores = 5
      )
    }

    cache_env$loaded[[name]] <- ps
    ps
  }
})
