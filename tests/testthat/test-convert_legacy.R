# testing the legacy workflow for loading the data from separate files
test_that("convert legacy", {

  ## test the .yaml file interface
  withr::with_tempdir({
    path <- file.path(
      system.file("extdata", package = "phiper"),
      "config.yaml"
    )

    ## basic run
    expect_no_error(suppressWarnings(
      phip_convert_legacy(config_yaml = path))
    )
  })

})
