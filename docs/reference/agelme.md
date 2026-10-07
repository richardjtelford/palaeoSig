# Estimation of the relationship between Calibrated age and depth

Estimates the relationship of Calibrated age and depth for
palaeorecords. The function uses a smooth spline from the mgcv library
by Simon Wood. It produces predicted confidence interval for the
relationship approximating a mixed effect model, as there are two levels
of uncertainty, i.e. within dated object and between dated objects.

## Usage

``` r
agelme(
  depup,
  depdo,
  bpup,
  bpdo,
  use,
  weights = c(1, rep(0, length(depup) - 1)),
  vspan = 1,
  k = length(depup) - 1,
  m = 2,
  diagnostic = FALSE
)
```

## Arguments

- depup:

  The upper depths of the dated slides

- depdo:

  The lower depths of the dated slides

- bpup:

  The younger calibrated ages of the dated slides

- bpdo:

  The older calibrated ages of the dated slides

- use:

  Logical vector of dates to include in the model. Default is to use
  all.

- weights:

  Weights to be used for the estimation, default is fixed top-layer
  followed by inverse variance of within dated object

- vspan:

  The span to be used for the diagnostic plots, default span = 1

- k:

  Number of base function to start the shrinkage in the gam estimation
  procedure

- m:

  The order of penalty for the term, i.e. the degree of continuity at
  the knots (default, m = 2 gives cubic smooth spline)

- diagnostic:

  Logical, should diagnostic plots be made.

## Value

A list of class 'agelme' with the following items:

- tdf :

  Degrees of freedom used by the cubic smooth spline, a vector with
  first value for constant variance and second vector for variance equal
  to mu.

- weights :

  A vector of the weights used by the cubic smooth spline

- Constant:

  A matrix with the numerical results for the dated points using a
  constant variance

- Muvar:

  A matrix with the numerical results for the dated points using
  variance equal to mu

- RES:

  A vector of the Residual sum of squares

- Models:

  A list with the models from the cubic smooth spline, constant and mu
  variance, respectively

- Data:

  A data.frame including the data used for the estimation

## Details

Note that the fixation of the top layer is done by a weight = 1, whereas
the other weights follows inverse variance within object.

The diagnostic plots is used to check the quality of the estimation and
to see if there is a need for an assumption of between object variance
proportional to mean. The latter however is rarely encountered for
palaeodata.

## References

Heegaard, E., Birks, HJB. & Telford, RJ. 2005. Relationships between
calibrated ages and depth in stratigraphical sequences: an estimation
procedure by mixed-effect regression. *The Holocene* **15**: 612-618
[doi:10.1191/0959683605hl836rr](https://doi.org/10.1191/0959683605hl836rr)

## Author

Einar Heegaard

## Examples

``` r
data(STOR)

fit.mod <- with(STOR, agelme(depthup, depthdo, cageup, cagedo))

# Predicting using the constant variance model,
# for each cm between 70 and 400 cm.
fit.pre <- predict(fit.mod, 1, 70:400)
plot(fit.pre)
```
