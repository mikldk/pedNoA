% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pedNoA.R
\name{pedNoA}
\alias{pedNoA}
\title{Number of shared alleles}
\usage{
pedNoA(x, ids, afreq, nindep = 0, marker_list = NULL, verbose = TRUE,
  nc = 1)
}
\arguments{
\item{x}{The pedigree - a linkdat object from \code{\link{paramlink}}.}

\item{ids}{A vector of indicies from \code{x}.}

\item{afreq}{A list of allele frequencies.}

\item{nindep}{The number of independent contributors.}

\item{marker_list}{A list of lists with known genotypes along with the respective locus index. See notes.}

\item{verbose}{If TRUE, information about \code{x} will be printet.}

\item{nc}{Number of cores to be used in parallelization,}
}
\description{
Calculates the probabilities of sharing \code{0,1,2,...} alleles between \code{ids} and \code{nindep} for the markers in \code{afreq}
}
\note{
The procedure is much faster when \code{length(ids) == 2}. See paper. The \code{marker_list} argument
is similar to that of \code{...} in the \code{marker} function in \code{paramlink}. If the the genotype
of individual \code{1} is known to be heterozygotic with allele \code{2} and \code{3} at locus \code{2}, and individual \code{2} is heterozygote with \code{1} and \code{2} at locus \code{3}, the argument
is \code{marker_list = list( list(2,3), list(1, c(2,2), 2, c(1,2)))}.
}
\examples{
# Uncle and nephew plus an independent individual:
library(paramlink)
x <- nuclearPed(2)
x <- addSon(x, 4)
plot(x)
p1 <- c(.1, .2, .3, .4)
p2 <- c(.3,.7)
pNx <- pedNoA(x, c(3,6), list(p1, p2), 1)
plot(pNx)

# First order cousins with grandfather being homozygote
# with allele 3 at maker \\code{2}  and grandmother being heterozygotic with allele 1 and 2 at marker \\code{3}
y <- cousinsPed(1)
p3 <- c(.2,.2,.3,.3)
plot(y)
pNy <- pedNoA(y, c(5,8), list(p2, p3), marker_list = list(list(2,3), list(1, c(3,3), 2, c(1,2))))
plot(pNy)

}
\seealso{
\code{\link{paramlink}}, \code{\link{plot.pedNoA}}
}
\author{
Mads Lindskou, \email{mads@math.aau.dk}
}
