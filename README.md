## MOSS: Multi-omic integration via sparse singular value decomposition.

Agustin Gonzalez-Reymundez, Alexander Grueneberg, and Ana I. Vazquez.

### Installing and loading MOSS

  The code below illustrates how to install and load MOSS.

```{r echo=T}
library(devtools)
install_github("agugonrey/MOSS")
library("MOSS")
```  

### Documentation.

  For a description of the package's main function. 

```{r echo=T}
help(moss)
```

  For more documentation, see the package's [vignette](https://github.com/agugonrey/MOSS/blob/master/vignettes/MOSS_working_example.pdf). An example of using MOSS on a multi-omic "big" pan-cancer data can be found [here](https://github.com/agugonrey/MOSS/blob/master/vignettes/MOSS_pancancer_example.pdf). The data can be found [here](https://data.mendeley.com/datasets/r8p67nfjc8/1).
