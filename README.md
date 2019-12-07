# Imperfect instrumental variable estimation

This package was born out of a term paper on the subject completed in early 2018. The package implements the imperfect instrumental variable estimation as seen in the excellent paper "Identification with Imperfect Instruments" by Aviv Nevo and Adam M. Rosen in Review of Economics and Statistics V. 94 2012.

## Content

Note that the Code is far from efficient or perfect. The idea was to have a reusable implementation. 

* The folder iivestimatoR contains the R-Package with functions that implement what is seen in the Paper.
* the example_implementation.Rmd file contains a mark-down file that implements the iiv identification with the widely available card dataset with distance to college as imperfect instrument.

## Download

To download the package, install devtools and then run:

`install_github("chtonios/IIV/iivestimatoR")`
