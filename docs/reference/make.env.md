# Generates simulated environmental variables

Simulates environmental variables used for generating species
abundances. Environmental variables may be correlated, and may follow
different distributions.

## Usage

``` r
make.env(n, elen, emean, edistr, ecor, ndim)
```

## Arguments

- n:

  Number of samples to be generated.

- elen:

  Range of the environmental variables. Single number or vector of
  length `ndim`.

- emean:

  Mean of the environmental variables. Single number or vector of length
  `ndim`.

- edistr:

  Distribution of the environmental variables. Currently 'uniform' and
  'Gaussian' are supported.

- ecor:

  Correlation matrix of the environmental variables. Object can be
  generated with
  [`cor.mat.fun`](https://richardjtelford.github.io/palaeoSig/reference/cor.mat.fun.md).
  If omitted environmental variables are not correlated.

- ndim:

  Number of environmental variables to generate.

## Value

Matrix of environmental variables. `n` rows and `ndim` columns.

## References

Minchin, P.R. (1987) Multidimensional Community Patterns: Towards a
Comprehensive Model. *Vegetatio*, **71**, 145-156.
[doi:10.1007/BF00039167](https://doi.org/10.1007/BF00039167)

## See also

[`cor.mat.fun`](https://richardjtelford.github.io/palaeoSig/reference/cor.mat.fun.md)

## Author

Mathias Trachsel and Richard J. Telford

## Examples

``` r
env.vars <- make.env(100,
  elen = rep(100, 10), emean = rep(50, 10),
  edistr = "uniform", ndim = 10
)
```
