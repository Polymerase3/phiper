#--check custom condition ------------------------------------------------------
.chk_cond <- function(condition, error_message, error = TRUE, ...) {
  if (condition && error) {
    chk::abort_chk(
      strwrap(error_message, prefix = " ", initial = ""),
      ...
    )
  } else if (condition && !error) {
    chk::wrn(
      strwrap(error_message, prefix = " ", initial = ""),
      ...
    )
  }
}

#--check if filename has given extension ---------------------------------------
.chk_extension <- function(name, x_name, ext_vec) {
  base <- basename(name)

  # split on dots, drop first part (filename), paste the rest
  ext <- strsplit(base, "\\.", fixed = FALSE)[[1]][-1]
  ext <- paste0(tolower(ext), collapse = ".")
  cond <- ext %in% ext_vec

  if (!cond) {
    chk::abort_chk(strwrap(sprintf(
      "The argument `%s` is allowed to have the following extensions: %s",
      x_name, word_list(add_quotes(ext_vec))
    ), prefix = " ", initial = ""))
  }
}

## --check null and replace with default ---------------------------------------
# NULL-coalescing helper
`%||%` <- function(x, y) if (!is.null(x)) x else y

.chk_null_default <- function(x, x_name, method, default) {
  .chk_cond(is.null(x),
    error = FALSE,
    sprintf(
      "The `%s` argument for the method %s was not provided
                    and will default to `%s`.",
      x_name,
      add_quotes(method),
      default
    )
  )

  if (is.null(x)) x <- default

  x
}

## --validate path to file -----------------------------------------------------
.chk_path <- function(path,
                      arg_name,
                      extension) {
  # validate if path is string of length 1
  .chk_cond(
    !chk::vld_string(path),
    sprintf(
      "The argument `%s` has to be a valid R string of length 1",
      arg_name
    )
  )

  # validate if file exists
  .chk_cond(
    !chk::vld_file(path),
    sprintf(
      "The filepath provided in the argument `%s`
                    does not exist: %s",
      arg_name, path
    )
  )

  # validate file extension if provided
  .chk_extension(
    path,
    arg_name,
    extension
  )
}

## --wordlists for error generation---------------------------------------------
word_list <- function(word_list = NULL, and_or = "and", is_are = FALSE,
                      quotes = FALSE) {
  # When given a vector of strings, creates a string of the form "a and b"
  # or "a, b, and c"
  # If is_are, adds "is" or "are" appropriately

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
      out <- sprintf(
        "%s %s %s",
        word_list[1L],
        and_or,
        word_list[2L]
      )
    } else {
      out <- sprintf(
        "%s, %s %s",
        paste(word_list[-len_wl], collapse = ", "),
        and_or,
        word_list[len_wl]
      )
    }
  }

  if (is_are) out <- sprintf("%s are", out)

  attr(out, "plural") <- TRUE

  out
}

add_quotes <- function(x, quotes = 2L) {
  if (isFALSE(quotes)) {
    return(x)
  }

  if (isTRUE(quotes)) {
    quotes <- '"'
  }

  if (chk::vld_string(quotes)) {
    return(paste0(quotes, x, quotes))
  }

  if (!chk::vld_count(quotes) || quotes > 2) {
    stop("`quotes` must be boolean, 1, 2, or a string.")
  }

  if (quotes == 0L) {
    return(x)
  }

  x <- {
    if (quotes == 1) {
      sprintf("'%s'", x)
    } else {
      sprintf('"%s"', x)
    }
  }

  x
}

`%nin%` <- function(x, inx) {
  !(x %in% inx)
}
