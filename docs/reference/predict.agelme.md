# Predicts the Calibrated age for agelme models

This function uses the output from 'agelme' to predict the Calibrated
ages for specified depths.

## Usage

``` r
# S3 method for class 'agelme'
predict(object, v = 1, depth, ...)
```

## Arguments

- object:

  An `agelme` model

- v:

  Using constant (1) or mu (2) variance

- depth:

  A vector of the depths to be predicted

- ...:

  Other arguments, currently unused.

## Value

A list with three items

- v:

  Whether constant variance or mu variance used.

- fit:

  A data.frame with the model predictions, containing columns:

  Depth

  : The depths for the predicted ages

  Estage

  : Predicted age

  Lowlim

  : Lower 95 % confidence interval

  Upplim

  : Upper 95 % confidence interval

  Tsd

  : Total standard deviation

&nbsp;

- data:

  A data.frame containing the age and depth information of the
  radiocarbon dates.

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
