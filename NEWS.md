# palaeoSig 2.0-1

Lots of changes in this version

## Breaking changes

* In the plotting functions `p.val` becomes `p_val` and is not the p-value rather than the 
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
