# Disconnect backend database connections

Closes any open database connections held inside a `phip_data` object,
including the `data_long` connection and the optional `peptide_library`
connection.

Call this when you want to release resources explicitly; it is also safe
to rely on garbage collection to close the connections automatically at
the end of the R session.

## Usage

``` r
disconnect(x)
```

## Arguments

- x:

  A valid `phip_data` object.

## Value

The input `phip_data` object, invisibly, with its connections closed.
