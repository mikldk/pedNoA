#' Helper function that packs the information of several pedigrees as input to \code{pedNoA_to_df}
#'
#' Calculates the probabilities of sharing \code{0,1,2,...} alleles between \code{ids} and \code{nindep} for the markers in \code{afreq} over several markers
#'
#' @param pedigrees A list with named pedigrees with class \code{linkdat}; see \code{paramlink}.
#' @param ids A list of vectors of indicies corresponding to \code{pedigrees}.
#' @param is_unrelated A list of the number of unrelated individuals corresponding to \code{pedigrees}.
#' @param markers A named list of marker frequencies.
#' @param convoluted A boolean indicating if markers should be convoluted
#'
#' @examples
#' # To come:
#'
#' @export
#'

pack_pedigree_information <- function(pedigrees, ids, is_unrelated, markers, convoluted = FALSE) {

    if(length(pedigrees) != length(ids)) stop("Unequal lengths")
    ## if(names(pedigrees) == NULL) stop("pedigrees must be named")
    ## if(names(markers) == NULL) stop("markers must be named")
    
    out <- list()
    if( convoluted ) {
        out <- lapply( seq_along(pedigrees), function(x) {
            list(pedigree = pedigrees[[x]],
                 pedigree_name = names(pedigrees)[x],
                 ids = ids[[x]],
                 unrelated = is_unrelated[[x]],
                 markers = markers)
        })
        class(out) <- append(class(out),"pack_pedigree_information_convoluted")
    }
    else {
        cnt <- 1
        for ( k in seq_along(pedigrees) ) {
            for ( j in seq_along(markers) ) {
                new_relation <- list(pedigree = pedigrees[[k]],
                                     pedigree_name = names(pedigrees)[k],
                                     ids = ids[[k]],
                                     unrelated = is_unrelated[[k]],
                                     marker = markers[j],
                                     marker_name = names(markers[j]))
                out[[cnt]] <- new_relation
                cnt <- cnt + 1
            }
        }
        class(out) <- append(class(out),"pack_pedigree_information")
    }
    out
}

#' Number of shared alleles in data.frame format for several markers
#'
#' Calculates the probabilities of sharing \code{0,1,2,...} alleles between pedigree members defined in \code{pack_pedigree_information} over several markers and possibly the convolution of these.
#'
#' @param pedigree_info A \code{ped}
#' @param ncores Number of cores to use in parallelization.

#' @author Mads Lindskou, \email{mads@@math.aau.dk}
#' @seealso \code{\link{paramlink}}, \code{\link{plot.pedNoA}}
#'
#' @examples
#' # To come...
#'
#' @export
#'

pedNoA_to_df <- function(pedigree_info, ncores) {
    
    if( "pack_pedigree_information_convoluted" %in% class(pedigree_info) ) {
        pedlist <- pbapply::pblapply(X = pedigree_info, cl = ncores, FUN = function(x) {
            z <- pedNoA(x$pedigree,
                        x$ids,
                        afreq = x$marker,
                        nindep = x$unrelated,
                        verbose = FALSE) ##, nc = ncores)
            pn <- z$pN
            counts <- as.numeric(names(pn))
            relation <- rep(x$pedigree_name, length(counts))
            data.frame(prob = pn, counts = counts, relation = relation)
        })
    }
    else {
        pedlist <- pbapply::pblapply(X = pedigree_info, cl = ncores, FUN = function(x) {
            z <- pedNoA(x$pedigree,
                        x$ids,
                        afreq = x$marker,
                        nindep = x$unrelated,
                        verbose = FALSE) ##, nc = ncores)
            pn <- z$pN
            counts <- as.numeric(names(pn))
            relation <- rep(x$pedigree_name, length(counts))
            marker <- rep(names(z$afreq), length(counts))
            data.frame(prob = pn, counts = counts, marker = marker, relation = relation)
        })
    }
    do.call(rbind, pedlist)
}
