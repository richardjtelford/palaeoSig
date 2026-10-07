# Calculate the effective number of species in the fossil data

Gives a measure of the species diversity in the fossil data.

## Usage

``` r
Hill.N2.core(spp)
```

## Arguments

- spp:

  data.frame of Species data

## Value

Minimum, first quartile and median effective number of species

## Details

Uses [`rioja::Hill.N2()`](https://rdrr.io/pkg/rioja/man/utils.html) from
the rioja package.

## Note

If the effective number of species is small, WA based reconstructions
are unlikely to be significant, and MAT based reconstructions should be
tested instead.

## References

Hill, M. O. (1973) Diversity and evenness: a unifying notation and its
consequences. *Ecology* **54**: 427–432.

## See also

[`rioja::Hill.N2()`](https://rdrr.io/pkg/rioja/man/utils.html)

## Author

Richard Telford

## Examples

``` r
require(rioja)
#> Loading required package: rioja
#> This is rioja  1.0-7
data(RLGH)
Hill.N2.core(RLGH$spec)
#>       0%      25%      50% 
#> 11.79277 13.11629 13.89376 
```
