# Anamorph

Creates functions that transform arbitrary distributions into Gaussian
distributions, and vice versa.

## Usage

``` r
anamorph(x, k, plot = FALSE)
```

## Arguments

- x:

  vector of data to transform

- k:

  number of Hermite polynomials to use

- plot:

  logical; plot the transformation?

## Value

Returns two function in a list

- xtog :

  Function to transform arbitrary variable x into a Gaussian
  distribution

- gtox :

  The back transformation

## Details

Increasing k can give a better fit.

## References

Wackernagel, H. (2003) *Multivariate Geostatistics.* 3rd edition,
Springer-Verlag, Berlin.
[doi:10.1007/978-3-662-05294-5](https://doi.org/10.1007/978-3-662-05294-5)

## Author

Richard Telford <Richard.Telford@bio.uib.no>

## Examples

``` r
set.seed(42)
x <- c(rnorm(50, 0, 1), rnorm(50, 6, 1))
hist(x)

ana.fun <- anamorph(x, 30, plot = TRUE)

xg <- ana.fun$xtog(x)
qqnorm(xg)
qqline(xg)

all.equal(x, ana.fun$gtox(xg))
#> [1] TRUE
```
