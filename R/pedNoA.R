.char_split <- function(x, pattern) unlist(strsplit(unlist(x), pattern))

.pfounder <- function(gf, af) {
  # Output: Joint probability of founder genotypes under HWE
  out = sapply(gf, function(z) {
    z = as.numeric(.char_split(z, "/"))
    ifelse(z[1] == z[2], af[z[1]]^2, 2*af[z[1]]*af[z[2]])
  })
  prod(out)
}

.conditional_marker <- function(marker_list, alleles, afreq, x){
  marker_data <- marker_list ## eval(substitute(marker_list))
  marker_data$alleles <- alleles
  marker_data$afreq <- afreq
  marker_data$x <- x
  do.call(paramlink::marker, marker_data)
}

.pedNoA2 <- function(x, ids, afreq, nindep = 0, m) {
    ## Computes the number of shared alleles when length(ids) = 2
    nalleles <- length(afreq)
    G <- paramlink::allGenotypes(nalleles)
    nG <- nrow(G)

    ## Picking out the (nG * (nG + 1) / 2) neccesarry joint genotypes
    id1 <- unlist(sapply(1:nG, function(i) i:nG))
    id2 <- unlist(sapply(1:nG, function(i) rep(i, nG-i+1)))
    grid <- matrix(c(id1, id2), ncol = 2)
    ngrid <- nrow(grid)

    ## Joint probability
    px <- paramlink::oneMarkerDistribution(x, ids, m, verbose = F, grid.subset = grid)
    if(nindep) {
        indep <- sapply(dimnames(px)[[1]], .pfounder, afreq)
        for( i in 1:nindep ) px <- outer(px, indep)
    }

    max_alleles <- min(nalleles, 2*(length(x$founders)+nindep), 2*(length(ids) + nindep))

    ## Outputvector
    pN <- structure(vector(length = max_alleles), names = 1:max_alleles)

    ## Summing over the neccesarry genotypes for the two pedigree members
    for( j in 1:ngrid ) {
        ped_index <- sapply(grid[j,], function(i) paste(G[i, ], collapse = "/"))
        ## All non-diagnoal elements must be multiplied by two
        px_coef <- ifelse(ped_index[1] == ped_index[2], 1, 2)

        if(!nindep) {
            N <- length(unique(.char_split(ped_index, "/")))
            pN[N] <- pN[N] + px_coef * px[matrix(ped_index, 1)]
        }

        if(nindep) {
            indep_index <- do.call(expand.grid, replicate(nindep, dimnames(px)[[1]], simplify = FALSE))
            ## Summing over the independent genotypes for each combination of the neccesarry pedigree genotypes
            for( k in 1:nrow(indep_index) ) {
                joint_index <- c(ped_index, c(as.matrix(indep_index[k,])))
                N <- length(unique(.char_split(joint_index, "/")))
                if( N <= max_alleles ) pN[N] <- pN[N] + px_coef * px[matrix(joint_index, 1)]
            }
        }
    }
    return(pN)
}

## .pedNoA_sim <- function(ped, ids, afreq, n_sim, unrelated) {
##     partial <- paramlink::marker(ped, afreq = afreq, alleles = seq_along(afreq))
##     sim <- paramlink::markerSim(ped, N=n_sim, partialmarker=partial, verbose = FALSE)
##     sim <- as.data.frame(sim)
##     sim <- as.matrix(sim[ids, 6:ncol(sim)])

##     if( unrelated ) {
##         for ( k in 1:unrelated ) {
##             new <- as.character(sample(seq_along(afreq),
##                                        2*n_sim,
##                                        prob = afreq,
##                                        replace = TRUE))
##             sim <- rbind(sim, new)
##         }
##     }

##     counts <- purrr::map_dbl( seq(2, ncol(sim), 2), function(x) {
##         q <- sim[, (x-1):x]
##         length(unique(c(q)))
##     })
    
##     table(counts) / sum(table(counts))
## }

.pedNoA_one_marker <- function(x, ids, afreq, nindep = 0, marker_list = NULL) {
    ## Compute the number of shared alleles between ids for a single marker
    nalleles <- length(afreq)
    if(!is.null(marker_list)) m <- .conditional_marker(marker_list, 1:nalleles, afreq, x)
    else m <- paramlink::marker(x, alleles = 1:nalleles, afreq = afreq)

    if( length(ids) == 2 ) return(.pedNoA2(x, ids, afreq, nindep, m))

    time <- Sys.time()
    px <- paramlink::oneMarkerDistribution(x, ids, m, verbose = F)

    if(nindep) {
        indep <- sapply(dimnames(px)[[1]], .pfounder, afreq)
        for( i in 1:nindep ) px <- outer(px, indep)
    }

    max_alleles <- min(length(afreq), 2*(length(x$founders)+nindep), 2*(length(ids) + nindep))
    pN <- structure(vector(length = max_alleles), names = 1:max_alleles)

    joint_index <- expand.grid(dimnames(px))
    for( j in 1:nrow(joint_index)) {
        index <- c(as.matrix((joint_index[j,])))
        N <- length(unique(.char_split(index, "/")))
        if( N <= max_alleles ) pN[N] <- pN[N] + px[matrix(index,1)]
    }
    return(pN)
}


#' Number of shared alleles
#'
#' Calculates the probabilities of sharing \code{0,1,2,...} alleles between \code{ids} and \code{nindep} for the markers in \code{afreq}
#'
#' @param x The pedigree - a linkdat object from \code{\link{paramlink}}.
#' @param ids A vector of indicies from \code{x}.
#' @param afreq A list of allele frequencies.
#' @param nindep The number of independent contributors.
#' @param marker_list A list of lists with known genotypes along with the respective locus index. See notes.
#' @param verbose If TRUE, information about \code{x} will be printet.
#' @param nc Number of cores to be used in parallelization,
#'
#' @note The procedure is much faster when \code{length(ids) == 2}. See paper. The \code{marker_list} argument
#' is similar to that of \code{...} in the \code{marker} function in \code{paramlink}. If the the genotype
#' of individual \code{1} is known to be heterozygotic with allele \code{2} and \code{3} at locus \code{2}, and individual \code{2} is heterozygote with \code{1} and \code{2} at locus \code{3}, the argument
#' is \code{marker_list = list( list(2,3), list(1, c(2,2), 2, c(1,2)))}.
#' @author Mads Lindskou, \email{mads@@math.aau.dk}
#' @seealso \code{\link{paramlink}}, \code{\link{plot.pedNoA}}
#'
#' @examples
#' # Uncle and nephew plus an independent individual:
#' library(paramlink)
#' x <- nuclearPed(2)
#' x <- addSon(x, 4)
#' plot(x)
#' p1 <- c(.1, .2, .3, .4)
#' p2 <- c(.3,.7)
#' pNx <- pedNoA(x, c(3,6), list(p1, p2), 1)
#' plot(pNx)
#'
#' # First order cousins with grandfather being homozygote
#' # with allele 3 at maker \code{2}  and grandmother being heterozygotic with allele 1 and 2 at marker \code{3}
#' y <- cousinsPed(1)
#' p3 <- c(.2,.2,.3,.3)
#' plot(y)
#' pNy <- pedNoA(y, c(5,8), list(p2, p3), marker_list = list(list(2,3), list(1, c(3,3), 2, c(1,2))))
#' plot(pNy)
#'
#' @export
#' @import paramlink
#'

pedNoA <- function(x, ids, afreq, nindep = 0, marker_list = NULL, verbose = TRUE, nc = 1){

    time <- Sys.time()
    x_copy <- x
    markers_index <- marker_list[[1]]

    if(!is.null(markers_index)) stopifnot(length(markers_index) != length(marker_list) - 1)

    pN_per_locus <- parallel::mclapply(mc.cores = nc, X = 1:length(afreq), FUN = function(l) {
        marker_l <- NULL
        x_l <- x_copy
        if( l %in% markers_index ) {
            marker_l <- unlist(marker_list[l+1], recursive = FALSE)
            na <- length(afreq[[l]])
            cm <- .conditional_marker(marker_l, 1:na, afreq[[l]], x)
            x_l <- paramlink::addMarker(x_l, cm)
            x <<- paramlink::addMarker(x, cm)
        }
        .pedNoA_one_marker(x_l, ids, afreq[[l]], nindep, marker_list = marker_l)
    })

    pN_convolution <- list(pN_per_locus[[1]])

    if( length(pN_per_locus) > 1 ) {
        for( locus in 2:length(pN_per_locus) ) {
            pNl <- c()
            n_alleles_total_up_to_locus <- sum(sapply(afreq[1:locus], length))
            n_min <- locus
            n_max <- min(2*(length(ids)+nindep)*length(afreq[1:locus]), n_alleles_total_up_to_locus)
            for( .n in n_min:n_max ) {
                I <- min(2*(length(ids)+nindep)*length(afreq), .n-1)
                n_i <- (.n - 1):1
                i <- 1:I
                n_i_not_allowed <- which(!(n_i %in% names(pN_convolution[[locus-1]])))
                i_not_allowed <- which(!(i %in% names(pN_per_locus[[locus]])))
                not_allowed <- unique(c(n_i_not_allowed, i_not_allowed))
                if( length(not_allowed) > 0 ){
                    n_i <- n_i[-not_allowed]
                    i <- i[-not_allowed]
                }
                ## Actual convolution
                .pNl <- mapply(function(x,y){
                    nx <- which(names(pN_convolution[[locus-1]]) == x) ## Neccesarry?
                    ny <- which(names(pN_per_locus[[locus]]) == y)
                    sum((pN_convolution[[locus-1]])[nx] * (pN_per_locus[[locus]])[ny])
                },
                n_i, i)
                add_pNl <- ifelse(length(.pNl) == 0, 0, sum(.pNl))
                pNl <- c(pNl, add_pNl)
            }
            pN_convolution[[(length(pN_convolution)+1)]] <- structure(pNl, names = n_min:n_max)
        }
    }
    time <- as.numeric(Sys.time() - time)
    out <- list(pN = pN_convolution[[length(pN_convolution)]],
                x = x,
                ids = ids,
                afreq = afreq,
                nindep = nindep,
                time = time)
    class(out) <- append(class(out),"pedNoA")
    if(verbose) {
        cat("\n ----------------------",
            "\n time: ", time,
            "\n m = ", ids,
            "\n #independent: ", nindep,
            "\n #alleles: ", sum(sapply(afreq, length)),
            "\n #markers: ", length(afreq),
            "\n E[N(m)=n]: ", sum(as.numeric(names(out$pN)) * out$pN),
            "\n ----------------------",
            "\n P(N(m)=n):\n")
        print(out$pN)
    }
    invisible(out)
}


#' A function for plotting the pedigree used in \code{pedNoA} along with information about which members and independent individuals are being tested.
#'
#' Individuals in \code{x} which are being tested are grayed out in the plot.
#'
#' @param x A pedNoA object from \code{\link{pedNoA}}.
#' @param title The title of the plot
#' @param ... Further arguments passed on to \code{plotPedList}
#' @author Mads Lindskou, \email{mads@@math.aau.dk}
#' @seealso \code{\link{pedNoA}}, \code{\link{pedNoA_}}
#'
#' @examples
#' # Uncle and nephew plus an independent individual:
#' library(paramlink)
#' x <- nuclearPed(2)
#' x <- addSon(x, 4)
#' plot(x)
#' p <- c(.1, .2, .3, .4)
#' pNx <-pedNoA(x, c(3,6), list(p), 1)
#' plot(pNx)
#'
#' # First order cousins with grandfather being homozygote
#' # with allele 3 and grandmother being heterozygotic with allele 1 and 2
#' y <- cousinsPed(1)
#' plot(y)
#' pNy <- pedNoA(y, c(5,8), list(p), marker_list = list(list(1), list(1, c(3,3), 2, c(1,2))))
#' plot(pNy)
#'
#' @export
#' @import paramlink
#'

plot.pedNoA <- function(x, title = "", ...) {
    ## Update, such that all marker information is displayed
    x$x <- setAvailable(x$x, x$ids)
    if(x$nindep) {
        margs_ped <- c(8,1,1,1)
        ped_list <- list(x$x, margins=margs_ped, available='shaded', symbolsize=1)
        s_list <- list()
        margs_singles <- c(0,0,0,0)
        for( j in 1:x$nindep ){
            id <- max(x$x$orig.ids)+j
            single <- singleton(id)
            single <- setAvailable(single, id)
            s_list[[j]] <- list(single, margins = margs_singles, available='shaded', symbolsize=1)
        }
        plotPedList(s_list, frames = F, available='shaded', dev.height = 5)
        graphics::par(new = TRUE)
        plotPedList(list(ped_list), frames = F, available='shaded', title = title, ...)
    }
    else {
        ped_list <- list(x$x, available='shaded', symbolsize = 1)
        plotPedList(list(ped_list), frames = F, available='shaded', title = title, ...)
    }
}

#' Printing a \code{pedNoA} object
#'
#' Printing information about the pedigree used in \code{pedNoA}
#'
#' @param x A pedNoA object from \code{\link{pedNoA}}.
#' @param ... Further arguments passed on to \code{print.default}
#' @author Mads Lindskou, \email{mads@@math.aau.dk}
#' @seealso \code{\link{pedNoA}}, \code{\link{pedNoA_}}
#'
#' @examples
#' # Uncle and nephew plus an independent individual:
#' library(paramlink)
#' x <- nuclearPed(2)
#' x <- addSon(x, 4)
#' plot(x)
#' p <- c(.1, .2, .3, .4)
#' pNx <-pedNoA(x, c(3,6), list(p), 1)
#' print(pNx)
#'
#' # First order cousins with grandfather being homozygote
#' # with allele 3 and grandmother being heterozygotic with allele 1 and 2
#' y <- cousinsPed(1)
#' plot(y)
#' pNy <- pedNoA(y, c(5,8), list(p), marker_list = list(list(1), list(1, c(3,3), 2, c(1,2))))
#' print(pNy)
#'
#' @export
#' @import paramlink
#'

print.pedNoA <- function(x, ...) {
    cat("\n ----------------------",
        "\n time: ", x$time,
        "\n m = ", x$ids,
        "\n #independent: ", x$nindep,
        "\n #alleles: ", sum(sapply(x$afreq, length)),
        "\n #markers: ", length(x$afreq),
        "\n E[N(m)=n]: ", sum(as.numeric(names(x$pN)) * x$pN),
        "\n ----------------------",
        "\n P(N(m)=n):\n")
    print(x$pN)
}
