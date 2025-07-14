# testing the legacy workflow for loading the data from separate files
test_that("convert legacy", {
  ## test the .yaml file interface
  withr::with_tempdir({
    path <- file.path(
      system.file("extdata", package = "phiper"),
      "config.yaml"
    )

    ## basic run - smoke test for all 3 backends
    ### memory
    expect_no_error(suppressWarnings(
      x <- phip_convert_legacy(
        config_yaml = path,
        backend = "arrow"
      )
    ))

    ### duckdb
    expect_no_error(suppressWarnings(
      x <- phip_convert_legacy(
        config_yaml = path,
        backend = "duckdb"
      )
    ))

    ### arrow
    expect_no_error(suppressWarnings(
      x <- phip_convert_legacy(
        config_yaml = path,
        backend = "arrow"
      )
    ))
  })
})
# x <- phip_convert_legacy(
#   config_yaml = path,
#   backend = "memory"
# )
# x
