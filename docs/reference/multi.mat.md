# MAT for multiple variables

MAT for many environmental variables simultaneously. More efficient than
calculating them separately for each variable.

## Usage

``` r
multi.mat(
  training.spp,
  envs,
  core.spp,
  noanalogues = 10,
  method = "sq-chord",
  run = "both"
)
```

## Arguments

- training.spp:

  Community data

- envs:

  Environmental variables - or simulations

- core.spp:

  Optional fossil data to make predictions for

- noanalogues:

  Number of analogues to use

- method:

  distance metric to use

- run:

  Return LOO predictions or predictions for fossil data

## Value

If `run = "both"`, a list with two elements:

- jack:

  Matrix of leave-one-out cross-validation predictions for the
  calibration set

- core:

  Matrix of predictions for the fossil data

Otherwise, one of these matrices is returned.

## References

Telford, R. J. and Birks, H. J. B. (2009) Evaluation of transfer
functions in spatially structured environments. *Quaternary Science
Reviews* **28**: 1309–1316.
[doi:10.1016/j.quascirev.2008.12.020](https://doi.org/10.1016/j.quascirev.2008.12.020)

## Author

Richard Telford <Richard.Telford@bio.uib.no>

## Examples

``` r
data(arctic.env)
data(arctic.pollen)

mMAT <- multi.mat(arctic.pollen, arctic.env[, 9:67], noanalogues = 5)
```
