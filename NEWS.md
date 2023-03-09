# palaeoSig 2.1-3

Vignettes work on a subset of data for speed.

# palaeoSig 2.1-2

Arguments to autoplot() functions have changed so first is now `object` rather than `x` to meet CRAN rules. This probably does not effect anybody unless they had named the first argument.

# palaeoSig 2.1-1

CRAN compliance version re CITATION

# palaeoSig 2.1-0

## Vignettes

New vignette showing how to use a red-noise null model for `randomTF()`.

Changed h-block cross-validation to use `sf` package for geographic data (rather than `sp`).

## Rewrites

`rne()` is now several orders of magnitude faster.

## bug fixes

In `randomTF()`, species and fossil data were automatically transformed to matrices so this was done once, not in every iteration.
This sometimes failed with MAT. 
Now this step is recommended in the help file.

Sometimes great-circle distances from a point to itself, calculated with `fields::rdist.earth()`, are non-zero due to rounding errors. 
This lead to them not being excluded in the h-block resampling in `rne()` when the threshold distance was zero. 
Now the diagonal of the distance matrix is set to zero, removing this problem.


# palaeoSig 2.0-8

## Bug fix

In randomTF(), if argument `col` was missing, all columns of the reconstruction would be used, because lazy evaluation was not done before use in a formula. Code now check `col` is set when needed.

# palaeoSig 2.0-7

Maintenance release to fix CRAN check error.

# palaeoSig 2.0-6

## Bug fix

MAT shortcut now works when autosim is set

# palaeoSig 2.0-1

Lots of changes in this version

## Breaking changes

* In the plotting functions `p.val` becomes `p_val` and is not the p-value rather than the percentile of the null distribution
* The functionality of `modelMaker` and `randomTFmm` (which can be used when a calibration set is used for many cores) is moved within `randomTF`

## Other changes

* Lot more sanity checks on inputs to `randomTF` and `obs.cor`
* `randomTF` now accepts a vector for the environmemtal data
* `randomTF` and `obs.cor` have a `permute` argument that permutes the environmental variable rather than simulating one from a uniform distribution. This is only possible with a single environmental variable.

## New functions

 * `autoplot` functions for the output of `randomTF` and `obs.cor`


## Deleted functions

 * `mat.h` Deleted as `rioja::MAT` can now do _h_-block cross-validation 
 * `modelMaker` and `randomTFmm` - this functionality is now implemented in `randomTF`


# Changes in version 1.1-3

## GENERAL

 * Bug fixes and minor improvement, merged package with age-depth functions and R implementation of Minchin's COMPAS. 
 * Changes to the namespace mean that the rioja package needs to be loaded explicitly before some of the function will work properly.
  
## NEW FUNCTIONS
  
 *  `agelme`, `predict.agelme` and `plot.fittedAgelme` age-depth modelling functions from Heegaard et al (2005)
 *  `coverage.plot` Diagnostic plot for reconstructions
 *  `centipede.plot` Plot WA optima and tolerances
 *  `make.env`, `species`, `make.env` Simulate species-environment relationships based on Minchin's (1983) compas.

## REMOVED FUNCTIONS
  
 *  `rotate` No longer needed as gstat now takes geodesic distance.
 *  `simulate.spatial` Was a wrapper for krige(), now easier to use krige() directly.

## MODIFIED FUNCTIONS
  
 *  `randomTF` and `randomTFmm` bug when partialling reconstructions out fixed


## BUG FIXES
  
 *  `obs.cor` bug when species not in same order in calibration set and fossil data fixed.
  

# Changes in version 1.1-2

## GENERAL
      
 *  Bug fixes and minor improvement.
      
  
## NEW FUNCTIONS
     
 *  `ModelMaker` and `randomTFmm` allows models for randomTF to be fitted once and used for several reconstructions. This can be much faster, but does not work with MAT.


## MODIFIED FUNCTIONS
  
 *  `plot.palaeoSig` p-value highlighted can be selected.
 *  `plot.obscor` p-value highlighted can be selected.
 *  `randomTF` Modified to accept predictions in a vector. 
      

## BUG FIXES
  
 *  `obs.cor` bug when species not in same order in calibration set and fossil data fixed.


# Changes in version 1.1-1

## GENERAL
      
 *  The package \pkg{autocorTF} has been merged with \pkg{palaeoSig} for ease of maintenance.

  
## NEW FUNCTIONS
     
 *  `identify.obscor` allows species names to be added to plots interactively.

## MODIFIED FUNCTIONS
    
 *  `obs.cor` Now includes correlations with several different species weights
 *  `plot.obscor` upgraded to allow the a choice of which abundance weighting is used, and the code for scaling points has been improved.
      

## BUG FIXES
  
 *  `randomTF` bug when partialling out reconstruction other than MAT now fixed.


# Changes in version 1.1

## GENERAL
      
 *  The package \pkg{autocorTF} has been merged with \pkg{palaeoSig} for ease of maintenance.
  
## NEW and REWRITTEN FUNCTIONS

      
 *  `RNE` replaces function `mat.rne`, `mat.he` and `mat.rd`. `RNE` allows any of the transfer function methods in \pkg{rioja} to be used to find the dependence of transfer function performance on spatially close observations.
 *  `plot.RNE` replaces `plot.rne` to work with the output of `RNE`.
 *  `obs.cor` replaces the previous version of `obs.cor`, `sim.cor` and `obscor.sig` to make the function easier to use and more similar to `randomTF`
 *  `plot.obscor` replaces the previous version of `plot.obscor` and `plot.simcor` to use the new output of `obs.cor`. The two plot types can be selected with the argument `which`.
 *  `jointsig`: tests if two environmental variables have joint control on fossil assemblage composition. 
 *  `plot.js`: a `plot` function for `jointsig`
 *  `Hill.N2.core` now calculates the Minimum, first quartile and median effective number of species for all fossil observations. This makes much more sense than the previous version, which conflated diversity in individual levels and turnover between levels.
 *  `plot.palaeoSig` has been improved to give neater figures.
  
## REMOVED FUNCTIONS

 * `wajack` was removed as `WA` is now implemented in \pkg{rioja}. 

## FORTHCOMING CHANGES
      
 * `mat.h` will be removed in the next release, when code for h-block resampling should be implemented in \pkg{rioja}.

