.onAttach <- function(libname, pkgname) {
  version <- utils::packageVersion("MOSS")
  pkg_message <- paste0("
 _____________________________________________________________________
|MOSS: Multi-Omic integration via Sparse Singular value decomposition.|
 _____________________________________________________________________

Agustin Gonzalez-Reymundez, Alexander Grueneberg, and Ana I. Vazquez.

Maintainer: Agustin Gonzalez-Reymundez <agugonrey@gmail.com>
URL: https://github.com/agugonrey/MOSS
BugReports: https://github.com/agugonrey/MOSS/issues
  
Version ",version,

".\n\nType 'help(moss)' for package overview and examples.
Type 'sim.data <- simulate_data()' to generate a small simulated data set.

To cite in published research:

 Agustin Gonzalez-Reymundez, Alexander Grueneberg, Guanqi Lu,
 Filipe Couto Alves, Gonzalo Rincon, Ana I Vazquez.
 
 'MOSS: Multi-omic integration with Sparse Value Decomposition'
 Bioinformatics (2022).

  https://doi.org/10.1093/bioinformatics/btac179

")
  packageStartupMessage(pkg_message)
}
