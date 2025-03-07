% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{switchReadStrand}
\alias{switchReadStrand}
\title{Switch read strand}
\usage{
switchReadStrand(strand)
}
\arguments{
\item{strand}{Character of read strand, either '+' or '-'.}
}
\value{
Returns character '+' or '-'.
}
\description{
Function to switch read strand when calling the Tagmentation reaction strand.
Reverse strands are switched to forward and vice versa, the output is
directly returned.
}
\details{
This function is used for \code{\link{mapInsertionSites}} where the direction
of one Tagmentation reaction can be determined from the reads of the opposite
reaciton
}
\author{
Koen Rademaker, \email{k.rademaker@nki.nl}
}
