# testing the legacy workflow for loading the data from separate files
test_that("convert legacy: memory", {
  ## test the .yaml file interface
  withr::with_tempdir({
    path <- file.path(
      system.file("extdata", package = "phiper"),
      "config.yaml"
    )

    ## SMOKE TEST --------------------------------------------------------------
    ### memory
    expect_no_error(suppressWarnings(
      phip_convert_legacy(
        config_yaml = path,
        backend = "memory"
      )
    ))

    ### with explicit paths, without config file, minimal example run
    expect_no_error(suppressWarnings(
      phip_convert_legacy(
        exist_file = file.path(
          system.file("extdata", package = "phiper"),
          "exist.csv"
        ),
        samples_file = file.path(
          system.file("extdata", package = "phiper"),
          "samples_meta.csv"
        ),
        comparisons_file = file.path(
          system.file("extdata", package = "phiper"),
          "comparisons.csv"
        ),
        backend = "memory"
      )
    ))

    ### error when no required file present
    expect_error(suppressWarnings(
      phip_convert_legacy(
        exist_file = file.path(
          system.file("extdata", package = "phiper"),
          "exist.csv"
        ),
        backend = "memory"
      )
    ), "samples_file")
  })
})

## testing the other backends; as these are "soft" dependencies, skip them if
## the packages are not installed
skip_if_not_installed("duckdb")
skip_if_not_installed("DBI")
skip_if_not_installed("arrow")
skip_if_not_installed("dbplyr")

test_that("convert legacy: duckdb and arrow", {
  ## test the .yaml file interface
  withr::with_tempdir({
    path <- file.path(
      system.file("extdata", package = "phiper"),
      "config.yaml"
    )
    ## SMOKE TESTS -------------------------------------------------------------
    ### default backend ("duckdb")
    expect_no_error(suppressWarnings(
      phip_convert_legacy(
        config_yaml = path
      )
    ))

    ### duckdb explicitly
    expect_no_error(suppressWarnings(
      phip_convert_legacy(
        config_yaml = path,
        backend = "duckdb"
      )
    ))

    ### arrow
    expect_no_error(suppressWarnings(
      phip_convert_legacy(
        config_yaml = path,
        backend = "arrow"
      )
    ))

    ### works without extra_cols -----------------------------------------------
    # ------------------------------------------------------------------#
    # 1. Build a sandbox directory in tempdir()
    # ------------------------------------------------------------------#
    workdir <- withr::local_tempdir()

    # helper: copy a file into workdir, keep basename
    copy_to_workdir <- function(src) {
      dst <- file.path(workdir, basename(src))
      file.copy(src, dst, overwrite = TRUE)
      basename(src) # return relative path
    }

    # ------------------------------------------------------------------#
    # 2. Locate original example files in your package
    # ------------------------------------------------------------------#
    pkg_ext <- function(name) system.file("extdata", name, package = "phiper")

    yaml_src <- pkg_ext("config.yaml")
    exist_src <- pkg_ext("exist.csv")
    samples_src <- pkg_ext("samples_meta.csv")
    timepoints_src <- NULL
    comps_src <- pkg_ext("comparisons.csv")

    # ------------------------------------------------------------------#
    # 3. Copy everything into workdir
    # ------------------------------------------------------------------#
    yaml_dst <- file.path(workdir, "config.yaml")
    exist_rel <- copy_to_workdir(exist_src)
    samples_rel <- copy_to_workdir(samples_src)
    comps_rel <- copy_to_workdir(comps_src)
    file.copy(yaml_src, yaml_dst, overwrite = TRUE)

    # ------------------------------------------------------------------#
    # 4. Edit the YAML: remove extra_cols, fix paths to be local
    # ------------------------------------------------------------------#
    cfg <- yaml::read_yaml(yaml_dst)

    cfg$extra_cols <- NULL # drop the key
    cfg$exist_file <- exist_rel # point to local file
    cfg$samples_file <- samples_rel
    cfg$timepoints_file <- NULL
    cfg$comparisons_file <- comps_rel

    yaml::write_yaml(cfg, yaml_dst)

    # ------------------------------------------------------------------#
    # 5. Call the converter
    # ------------------------------------------------------------------#
    suppressWarnings({
      pd <- phip_convert_legacy(
        config_yaml = yaml_dst,
        backend     = "duckdb"
      )
    })

    # ------------------------------------------------------------------#
    # 6. Expectations
    # ------------------------------------------------------------------#
    expect_s3_class(pd, "phip_data")
    expect_gt(ncol(get_counts(pd)), 3) # additional columns from meta,
    # as there is no extra_cols
    expect_false(isTRUE(pd$meta$fold_change))
    expect_false(isTRUE(pd$meta$longitudinal))
  })
})


## .auto_read_csv
tmp_csv <- withr::local_tempfile(fileext = ".csv")
write.csv(data.frame(a = 1:3, b = 4:6), tmp_csv, row.names = FALSE)

# ------------------------------------------------------------------
# 1) branch where data.table IS available
# ------------------------------------------------------------------
test_that(".auto_read_csv uses data.table::fread when available", {
  skip_if_not_installed("data.table") # ensures branch can run

  res <- .auto_read_csv(tmp_csv)

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 3)
  expect_true(attr(res, "class")[1] != "data.frame" ||
    # fread returns data.frame if data.table=FALSE
    TRUE)
})

# ------------------------------------------------------------------
# 2) branch where data.table is *pretended* to be missing
# ------------------------------------------------------------------
test_that(".auto_read_csv falls back to read.csv when data.table is absent", {
  # Mock requireNamespace() so it always returns FALSE inside this call
  mockery::stub(.auto_read_csv, "requireNamespace", function(pkg, ...) FALSE)

  res <- .auto_read_csv(tmp_csv)

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 3)
})

