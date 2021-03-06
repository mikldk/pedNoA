% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pedNoA_to_df.R
\name{pedNoA_to_df}
\alias{pedNoA_to_df}
\title{Number of shared alleles in data.frame format for several markers}
\usage{
pedNoA_to_df(pedigree_info, ncores)
}
\arguments{
\item{pedigree_info}{A \code{ped}}

\item{ncores}{Number of cores to use in parallelization.}
}
\description{
Calculates the probabilities of sharing \code{0,1,2,...} alleles between pedigree members defined in \code{pack_pedigree_information} over several markers and possibly the convolution of these.
}
\examples{
# To come...

}
\seealso{
\code{\link{paramlink}}, \code{\link{plot.pedNoA}}
}
\author{
Mads Lindskou, \email{mads@math.aau.dk}
}
