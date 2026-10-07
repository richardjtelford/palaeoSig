# Generate correlation matrix

Generates a correlation matrix for the environmental variables generated
in
[`make.env`](https://richardjtelford.github.io/palaeoSig/reference/make.env.md)
and for correlated species optima in
[`species`](https://richardjtelford.github.io/palaeoSig/reference/species.md).
Only used when correlated environmental variables or optima are
generated.

## Usage

``` r
cor.mat.fun(ndim, cors)
```

## Arguments

- ndim:

  Number of environmental variables that are subsequently generated with
  [`make.env`](https://richardjtelford.github.io/palaeoSig/reference/make.env.md).

- cors:

  List of correlations between environmental variables. Each element of
  the list consists of three numbers, the first two numbers indicate the
  variables that are correlated, the third number is the correlation
  coefficient. If correlations between two variables are omitted the
  correlation remains 0.

## Value

A correlation matrix

## See also

[`make.env`](https://richardjtelford.github.io/palaeoSig/reference/make.env.md),
[`species`](https://richardjtelford.github.io/palaeoSig/reference/species.md)

## Author

Mathias Trachsel

## Examples

``` r
correlations <- list(c(1, 2, 0.5), c(1, 4, 0.1), c(2, 5, 0.6))
cor.mat <- cor.mat.fun(5, correlations)
```
