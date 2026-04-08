# Theme `theme_phip`

A clean, publication-ready ggplot2 theme tuned for **facetted** plots
with the **Montserrat** font. The font is registered and **showtext**
rendering is enabled automatically when the package loads — no setup
required.

## Usage

``` r
theme_phip(base_size = 14, base_family = "Montserrat")
```

## Arguments

- base_size:

  Base font size.

- base_family:

  Base font family (default `"Montserrat"`).

## Value

A ggplot2 `theme` object.

## See also

Other phip-ggplot:
[`scale_colour_phip()`](https://polymerase3.github.io/phiper/reference/scale_colour_phip.md),
[`scale_fill_phip()`](https://polymerase3.github.io/phiper/reference/scale_fill_phip.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ggplot2::ggplot(iris, ggplot2::aes(Sepal.Length, Sepal.Width, colour = Species)) +
  ggplot2::geom_point() +
  scale_colour_phip() +
  theme_phip()
} # }
```
