% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/startMetaviz.R
\name{startMetavizStandalone}
\alias{startMetavizStandalone}
\title{Start metaviz app in standalone (locally) and create \code{\link[metavizr]{MetavizApp}} object to manage connection.}
\usage{
startMetavizStandalone(register_function = .register_all_metaviz_things,
  use_viewer_option = FALSE, ...)
}
\arguments{
\item{register_function}{(function) function used to register actions and charts on the metaviz app.}

\item{use_viewer_option}{(function) run application in viewer defined by \code{getOption("viewer")}.
This allows standalone app to run in Rstudio's viewer (FALSE by default)}

\item{...}{additional parameters passed to \code{\link[epivizrStandalone]{startStandalone}}.}
}
\value{
An object of class \code{\link{MetavizApp}}
}
\description{
Start metaviz app in standalone (locally) and create \code{\link[metavizr]{MetavizApp}} object to manage connection.
}
\examples{

#' # see package vignette for example usage
app <- startMetavizStandalone(non_interactive=TRUE)
app$stop_app()


}
