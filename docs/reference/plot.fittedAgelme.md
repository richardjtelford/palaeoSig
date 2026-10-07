# Plots fitted agelme model and dates

Plots fitted agelme model and dates

## Usage

``` r
# S3 method for class 'fittedAgelme'
plot(x, main, xlab = "Depth", ylab = "Calibrated Age", ...)
```

## Arguments

- x:

  Fitted agelme model.

- main:

  Title of the plot.

- xlab:

  x axis label of the plot.

- ylab:

  y axis label of the plot.

- ...:

  Other arguments to plot.

## Examples

``` r
data(STOR)
fit.mod <- with(STOR, agelme(depthup, depthdo, cageup, cagedo))
#Predicting using the constant variance model,
#for each cm between 70 and 400 cm.
fit.pre <- predict(fit.mod, 1, 70:400)
plot(fit.pre)
```
