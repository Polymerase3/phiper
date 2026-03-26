# Path to example PhIP-Seq datasets shipped with phiper

Path to example PhIP-Seq datasets shipped with phiper

## Usage

``` r
phip_example_path(name = c("phip_mixture"))
```

## Arguments

- name:

  Character scalar. Name of the example dataset. Currently supported:
  `"phip_mixture"`.

## Value

A character scalar with an absolute path to the file.

## Examples

``` r
sim_path <- phip_example_path("phip_mixture")
# phip_obj <- phip_convert(sim_path)
```
