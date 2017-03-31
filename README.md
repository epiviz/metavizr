metavizR
========

The `metavizr` package implements two-way communication between the [R/Bioconductor](http://bioconductor.org) environment and the [metaviz](http://metaviz.cbcb.umd.edu) web app for interactive visualization of microbiome sequencing results. The hierarchy of features from a microbiome sequencing result can be visualized with a navigation utility and count values are displayed dynamically updated heatmaps or stacked bar plots. Metavizr uses Websockets for communication between the browser Javascript client and the R environment.

## Installation and requirements
metaivizr is available from github. To install `metavizr`:

```{r}
library(devtools)
install_github("epiviz/metavizr", build_vignettes = TRUE)
```

## Try it out

The easiest way to try `metavizr` is to follow the package vignette:

```{r}
require(metvizr)
browseVignettes("metavizr")
```

## Tutorials

You can get a quick tour of metaviz here: [https://epiviz.github.io/tutorials/metaviz/overview/](https://epiviz.github.io/tutorials/metaviz/overview/)

## Non-blocking

Metavizr supports a non-blocking workflow on both UNIX-like and Windows systems where data is served to the webapp without blocking the R/bioc interactive session. Make sure you are using the latest version of the [httpuv package](http://cran.r-project.org/web/packages/httpuv/index.html) to use this. (Thanks to the[Rstudio](http://rstudio.org) folks for folding our daemonizing code into the main httpuv release).

## More info

Check out the `Metaviz` [http://metaviz.org](http://metaviz.org) for more details and documentation.
