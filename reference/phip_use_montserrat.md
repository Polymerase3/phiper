# Register and enable the Montserrat font for plotting (showtext)

Registers the Google font **Montserrat** (via **sysfonts**) and enables
**showtext** so the font renders consistently in all devices (PNG, PDF,
etc.). Safe to call multiple times; silently no-ops if already active.

## Usage

``` r
phip_use_montserrat(family = "Montserrat", enable = TRUE)
```

## Arguments

- family:

  Internal family alias to register. Default `"Montserrat"`.

- enable:

  If `TRUE`, turns on `showtext::showtext_auto(TRUE)`.

## Value

(Invisibly) the family name to use in themes.

## Details

Requires packages **showtext** and **sysfonts**. If the Google download
fails (e.g., offline), the theme will still try to use a locally
installed Montserrat by name; otherwise it falls back to the default
device font.

## See also

Other phip-ggplot:
[`scale_colour_phip()`](https://polymerase3.github.io/phiper/reference/scale_colour_phip.md),
[`theme_phip()`](https://polymerase3.github.io/phiper/reference/theme_phip.md)
