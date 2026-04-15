# Discrete fill scale using the PHIP palette

A thin wrapper around
[`ggplot2::scale_fill_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html)
that applies `phip_palette`. Set as the session default on
[`library(phiper)`](https://github.com/Polymerase3/phiper).

## Usage

``` r
scale_fill_phip(...)
```

## Arguments

- ...:

  Arguments passed on to
  [`discrete_scale`](https://ggplot2.tidyverse.org/reference/discrete_scale.html)

  `limits`

  :   One of:

      - `NULL` to use the default scale values

      - A character vector that defines possible values of the scale and
        their order

      - A function that accepts the existing (automatic) values and
        returns new ones. Also accepts rlang
        [lambda](https://rlang.r-lib.org/reference/as_function.html)
        function notation.

  `drop`

  :   Should unused factor levels be omitted from the scale? The
      default, `TRUE`, uses the levels that appear in the data; `FALSE`
      includes the levels in the factor. Please note that to display
      every level in a legend, the layer should use
      `show.legend = TRUE`.

  `na.translate`

  :   Unlike continuous scales, discrete scales can easily show missing
      values, and do so by default. If you want to remove missing values
      from a discrete scale, specify `na.translate = FALSE`.

  `name`

  :   The name of the scale. Used as the axis or legend title. If
      `waiver()`, the default, the name of the scale is taken from the
      first mapping used for that aesthetic. If `NULL`, the legend title
      will be omitted.

  `minor_breaks`

  :   One of:

      - `NULL` for no minor breaks

      - `waiver()` for the default breaks (none for discrete, one minor
        break between each major break for continuous)

      - A numeric vector of positions

      - A function that given the limits returns a vector of minor
        breaks. Also accepts rlang
        [lambda](https://rlang.r-lib.org/reference/as_function.html)
        function notation. When the function has two arguments, it will
        be given the limits and major break positions.

  `labels`

  :   One of the options below. Please note that when `labels` is a
      vector, it is highly recommended to also set the `breaks` argument
      as a vector to protect against unintended mismatches.

      - `NULL` for no labels

      - `waiver()` for the default labels computed by the transformation
        object

      - A character vector giving labels (must be same length as
        `breaks`)

      - An expression vector (must be the same length as breaks). See
        ?plotmath for details.

      - A function that takes the breaks as input and returns labels as
        output. Also accepts rlang
        [lambda](https://rlang.r-lib.org/reference/as_function.html)
        function notation.

  `guide`

  :   A function used to create a guide or its name. See
      [`guides()`](https://ggplot2.tidyverse.org/reference/guides.html)
      for more information.

  `call`

  :   The `call` used to construct the scale for reporting messages.

  `super`

  :   The super class to use for the constructed scale

## Value

A ggplot2 scale.

## See also

Other phip-ggplot:
[`scale_colour_phip()`](https://polymerase3.github.io/phiper/reference/scale_colour_phip.md),
[`theme_phip()`](https://polymerase3.github.io/phiper/reference/theme_phip.md)

## Examples

``` r
ggplot2::ggplot(iris,
  ggplot2::aes(Species, Sepal.Length, fill = Species)) +
  ggplot2::geom_boxplot() +
  scale_fill_phip()
```
