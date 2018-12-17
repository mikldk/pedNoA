.nsibs <- function(n, afreq) {
    k1 <- 2^(2-2*n)
    k2 <- (2^(2-n) - 2^(3-2*n))
    k3 <- (1 - 2^(1-n))^2
    m2 <- sum(afreq^2)
    m3 <- sum(afreq^3)
    m4 <- sum(afreq^4)
    p1 <- k1*m2 + k2*m3 + k3*m4
    p2 <- k1 + (3*k2 - k1)*m2 + (4*k3 - 3*k2)*m3 - 7*k3*m4 + 3*k3*m2^2
    p3 <- k2 + (6*k3 - 3*k2)*m2 + (2*k2 - 12*k3)*m3 + 12*k3*m4 - 6*k3*m2^2
    p4 <- k3*(1 - 6*m2 + 8*m3 - 6*m4 + 3*m2^2)
    structure(c(p1,p2,p3,p4), names = 1:4)
}

.uncle_nephew <-function(afreq) {
    m2 <- sum(afreq^2)
    m3 <- sum(afreq^3)
    m4 <- sum(afreq^4)
    p1 <- .5 * (m3 + m4)
    p2 <- .5 * (3*m2 + m3 - 7*m4 + 3*m2^2)
    p3 <- .5 * (1 + 3*m2 - 10*m3 + 12*m4 - 6*m2^2)
    p4 <- .5 * (1 - 6*m2 + 8*m3 - 6*m4 + 3*m2^2)
    structure(c(p1,p2,p3,p4), names = 1:4)
}

.first_cousins <- function(afreq) {
    m2 <- sum(afreq^2)
    m3 <- sum(afreq^3)
    m4 <- sum(afreq^4)
    p1 <- .25 * ( m3 + 3*m4 )
    p2 <- .25 * ( 3*m2 + 9*m3 - 21*m4 + 9*m2^2)
    p3 <- .5 * (.5 + 15/2*m2 - 17*m3 + 18*m4 - 9*m2^2)
    p4 <- 3/4 * (1 - 6*m2 + 8*m3 - 6*m4 + 3*m2^2)
    structure(c(p1,p2,p3,p4), names = 1:4)
}

#' Number of shared alleles
#'
#' Calculates the probabilities of sharing \code{0,1,2,...} alleles between members with the relation \code{relation}
#'
#' @param relation A valid relation; see notes.
#' @param afreq A vector of allele frequencies.
#' @param verbose If TRUE, information about \code{relation} will be printet.
#'
#' @note Valid relations are:
#' "nsibs" with \deqn{n \geq 1} for \code{n} siblings.
#' "uncle_nephew" for uncle and nephew
#' "first_cousins" for first cousins.
#'
#' @author Mads Lindskou, \email{mads@@math.aau.dk}
#' @seealso \code{\link{paramlink}}, \code{\link{plot.pedNoA_}}
#'
#' @examples
#' # Uncle and nephew plus an independent individual:
#' library(paramlink)
#' x <- nuclearPed(2)
#' x <- addSon(x, 4)
#' plot(x)
#' p <- c(.1, .2, .3, .4)
#' pNx <- pedNoA(x, c(3,6), p, 1)
#' plot(pNx)
#'
#' # First order cousins with grandfather being homozygote
#' # with allele 3 and grandmother being heterozygotic with allele 1 and 2
#' y <- cousinsPed(1)
#' plot(y)
#' pNy <- pedNoA(y, c(5,8), p, marker_list = list(1, c(3,3), 2, c(1,2)))
#' plot(pNy)
#'
#' @export
#' @import paramlink
#'

pedNoA_ <- function(relation, afreq, verbose = T) {

    ## if( !isTRUE(all.equal( 1, sum(afreq))) ) stop("Allele frequencies do not sum to one!")

    if( grepl("sibs", relation) ) {
        n <- as.numeric(gsub("sibs", "", relation))
        if( is.na(n) ) stop("Not a valid relation")
        relation <- "nsibs"
        x <- paramlink::nuclearPed(n)
        ids <- 3:(n+2)
    }

    if( relation == "uncle_nephew" ) {
        x <- paramlink::nuclearPed(2)
        x <- paramlink::addSon(x, 4)
        ids <- c(3,6)
    }

    if( relation == "first_cousins" ) {
      x <- paramlink::cousinsPed(1)
      ids <- c(5,8)
    }

    pN <- switch(relation,
                 nsibs = .nsibs(n, afreq),
                 uncle_nephew  = .uncle_nephew(afreq),
                 first_cousins = .first_cousins(afreq),
                 stop("Not a valid relation")
                 )

    out <- list(pN = pN, x = x, ids = ids, afreq = afreq)
    class(out) <- append(class(out),"pedNoA_")
    if(verbose) {
        cat("\n m =", ids, "\n #alleles: ",
            length(afreq), "\n P(N(m)=n):", "\n")
        print(pN)
    }
    invisible(out)
}

#' A function for plotting the pedigree used in \code{pedNoA_}
#'
#' Individuals in \code{x} which are being tested are grayed out in the plot.
#'
#' @param x A \code{pedNoA_} object
#' @param title The title of the plot
#' @param ... Further arguments passed on to \code{plotPedList}
#' @author Mads Lindskou, \email{mads@@math.aau.dk}
#' @seealso \code{\link{pedNoA}}, \code{\link{pedNoA_}}
#'
#' @examples
#' # 5 sibslings:
#' library(paramlink)
#  p <- rchisq(50, 20)
#  p <- p/sum(p)
#' pNx <- pedNoA_("5sibs", p)
#' plot(pNx)
#'
#' # Probability of being het or hom:
#' pN_1 <- pedNoA_("1sibs", p)
#' plot(pN_1)
#'
#' # Uncle and nephew
#' pNy <- pedNoA_("uncle_nephew", p)
#' plot(pNy)
#'
#' @export
#' @import paramlink
#'

plot.pedNoA_ <- function(x, title = "", ...){
    x$x <- setAvailable(x$x, x$ids)
    ped_list <- list(x$x, available='shaded', symbolsize = 1)
    plotPedList(list(ped_list), frames = F, available='shaded', title = title, ...)
}

#' Printing a \code{pedNoA_} object
#'
#' Printing information about the pedigree used in \code{pedNoA_}
#'
#' @param x A pedNoA object from \code{\link{pedNoA}}
#' @param ... Further arguments passed on to \code{print.default}
#' @author Mads Lindskou, \email{mads@@math.aau.dk}
#' @seealso \code{\link{pedNoA}}, \code{\link{pedNoA_}}
#'
#' @examples
#' # 5 sibslings:
#' library(paramlink)
#  p <- rchisq(50, 20)
#  p <- p/sum(p)
#' pNx <-pedNoA_("5sibs", p)
#' plot(pNx)
#'
#' # Probability of being het or hom:
#' pN_1 <- pedNoA_("1sibs", p)
#' plot(pN_1)
#'
#' # Uncle and nephew
#' pNy <- pedNoA_("uncle_nephew", p)
#' plot(pNy)
#'
#' @export
#' @import paramlink
#'

print.pedNoA_ <- function(x, ...) {
    cat("\n m =", x$ids, "\n #alleles: ",
        length(x$afreq), "\n P(N(m)=n):", "\n")
    print(x$pN)
}
