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

