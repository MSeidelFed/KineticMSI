% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AcqTimeRandomSamplingMSI.R
\name{AcqTimeRandomSamplingMSI}
\alias{AcqTimeRandomSamplingMSI}
\title{A function to sort previous to unification independent kMSI samples according to pixel acquisition time}
\usage{
AcqTimeRandomSamplingMSI(
  path,
  PatternEnrichment = "MeanEnrichment.csv",
  RepsOrPath = c("reps", "path"),
  featureChoice = 0
)
}
\arguments{
\item{path}{parent directory were all enrichment files are contained, digs within recursive folders.}

\item{PatternEnrichment}{defaults to "MeanEnrichment". Inherits from parent function.}

\item{RepsOrPath}{defines whether the input is a vector with paths to individual reps or already a path itself.}

\item{featureChoice}{Allows selecting the index of a feature of interest to be sorted.}
}
\description{
This function allows KineticMSI functions to build a matrix with independent samples by unifying them.
}
\examples{

...
}
\keyword{Dynamics}
\keyword{MSI}
\keyword{Subsetting}
\keyword{Tracer}
