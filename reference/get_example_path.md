# Path to Example PhIP-Seq Datasets

Return the path to an example dataset shipped with phiperio, suitable
for use with
[`load_example_data`](https://polymerase3.github.io/phiper/reference/load_example_data.md)
or
[`convert_standard`](https://polymerase3.github.io/phiperio/reference/convert_standard.html).

## Usage

``` r
get_example_path(name = c("phip_mixture"))
```

## Arguments

- name:

  Character scalar. Name of the example dataset. Currently supported:
  `"phip_mixture"`.

## Value

A character scalar with an absolute path to the file.

## Examples

``` r
sim_path <- get_example_path("phip_mixture")
# phip_obj <- convert_standard(sim_path)
```
