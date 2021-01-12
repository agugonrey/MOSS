.onAttach <- function(libname, pkgname) {
  packageStartupMessage("
 _____________________________________________________________________
|MOSS: Multi-Omic integration via Sparse Singular value decomposition.|
 _____________________________________________________________________

Agustin Gonzalez-Reymundez, Alexander Grueneberg, and Ana I. Vazquez.

Maintainer: Agustin Gonzalez-Reymundez <agugonrey@gmail.com>
URL: https://github.com/agugonrey/MOSS
BugReports: https://github.com/agugonrey/MOSS/issues

Type 'help(moss)' for package overview and examples.
Type 'sim.data <- simulate_data()' to generate a small simulated data set.
")
}
