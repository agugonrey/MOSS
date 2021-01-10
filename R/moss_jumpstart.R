.onAttach <- function(libname, pkgname) {
  packageStartupMessage("
 _____________________________________________________________________
|MOSS: Multi-omic integration via sparse singular value decomposition.|
 _____________________________________________________________________

Agustin Gonzalez-Reymundez, Alexander Grueneberg, and Ana I. Vazquez.

Maintainer:                                   
  Agustin Gonzalez-Reymundez <gonza650@msu.edu>,                    

Documentation:
  Type 'help(moss)' for package overview and examples.

Type 'sim.data <- simulate_data()',
to generate a small simulated data set with three omics and 
one categorical variable. 
")
}
