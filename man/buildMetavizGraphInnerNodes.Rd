% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MetavizGraphInnerNodes-class.R
\docType{methods}
\name{buildMetavizGraphInnerNodes}
\alias{buildMetavizGraphInnerNodes}
\alias{buildMetavizGraphInnerNodes,MRexperiment-method}
\title{Build a MetavizTree object from another object}
\usage{
buildMetavizGraphInnerNodes(object, ...)

\S4method{buildMetavizGraphInnerNodes}{MRexperiment}(object, feature_order,
  ...)
}
\arguments{
\item{object}{The object from which taxonomy data is extracted}

\item{...}{Additional arguments}

\item{feature_order}{Ordering of leaves (features) in taxonomy tree}
}
\value{
a \code{\link{MetavizGraphInnerNodes}} object
}
\description{
Build a MetavizTree object from another object
}
\section{Methods (by class)}{
\itemize{
\item \code{MRexperiment}: Build graph from a \code{\link[metagenomeSeq]{MRexperiment-class}} object
}}

