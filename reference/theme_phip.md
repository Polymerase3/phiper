# Theme `theme_phip`

A clean, publication-ready ggplot2 theme inspired by the provided
`theme_Publication` snippet, tuned for **facetted** plots and consistent
use of the **Montserrat** font (register it with
[`phip_use_montserrat()`](https://polymerase3.github.io/phiper/reference/phip_use_montserrat.md)).

## Usage

``` r
theme_phip(base_size = 14, base_family = "Montserrat")
```

## Arguments

- base_size:

  Base font size.

- base_family:

  Base font family (default `"Montserrat"`). Call
  [`phip_use_montserrat()`](https://polymerase3.github.io/phiper/reference/phip_use_montserrat.md)
  once per session to register and enable rendering.

## Value

A ggplot2 `theme` object.

## See also

Other phip-ggplot:
[`phip_use_montserrat()`](https://polymerase3.github.io/phiper/reference/phip_use_montserrat.md),
[`scale_colour_phip()`](https://polymerase3.github.io/phiper/reference/scale_colour_phip.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Register Montserrat once per session
phip_use_montserrat()

ggplot2::ggplot(iris, ggplot2::aes(Sepal.Length, Sepal.Width, colour = Species)) +
  ggplot2::geom_point() +
  scale_colour_phip() +
  theme_phip()
} # }
```
