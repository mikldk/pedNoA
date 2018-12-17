## Remove (located in pedNoA)
.char_split <- function(x, pattern) unlist(strsplit(unlist(x), pattern))

## Remove (located in pedNoA)
.conditional_marker <- function(marker_list, alleles, afreq, x){
  marker_data <- marker_list ## eval(substitute(marker_list))
  marker_data$alleles <- alleles
  marker_data$afreq <- afreq
  marker_data$x <- x
  do.call(paramlink::marker, marker_data)
}

.pedNoA_zscore_sim_onemarker <- function(x,
                                         idsx,
                                         y = x, ## the sample DNA mixture
                                         idsy = idsx,
                                         afreq, ## a SNP vector
                                         nindep = 0,
                                         marker_list = NULL) {
    nalleles <- length(afreq)
    if(!is.null(marker_list)) m <- .conditional_marker(marker_list, 1:nalleles, afreq, x)
    else m <- paramlink::marker(x, alleles = 1:nalleles, afreq = afreq)

    px <- paramlink::oneMarkerDistribution(x, idsx, m, verbose = F)

    if(nindep) {
        indep <- sapply(dimnames(px)[[1]], .pfounder, afreq)
        for( i in 1:nindep ) px <- outer(px, indep)
    }

    joint_index <- expand.grid(dimnames(px))
    n <- nrow(joint_index)

    hom_het_list <- sapply(1:n, function(x) {
        all_alleles <-.char_split(as.character(unlist(joint_index[x,])), "/")
        out <- ""
        if( length(unique(all_alleles)) == 1 ) {
            out <- all_alleles[1]
        }
        out
    })

    p1 <- px[as.matrix(joint_index[which(hom_het_list == "1"),])]
    p2 <- px[as.matrix(joint_index[which(hom_het_list == "2"),])]
    p  <- structure(c(p1, p2, 1-p1-p2), names = 1:3)

    partial <- paramlink::marker(y, afreq = afreq, alleles = seq_along(afreq))
    sim <- paramlink::markerSim(y, N = 1, partialmarker = partial, verbose = FALSE)
    sim <- as.data.frame(sim)
    sim <- as.matrix(sim[idsy, 6:ncol(sim)])

    if( nindep ) {
        for ( k in 1:nindep ) {
            new <- as.character(sample(seq_along(afreq),
                                       2*nsim,
                                       prob = afreq,
                                       replace = TRUE))
            sim <- rbind(sim, new)
        }
    }

    x0 <- c(0,0,0)
    if( length(unique(c(sim))) == 1 && c(sim)[1] == "1"  ) {
        x0[1] <- 1
    } else if ( length(unique(c(sim))) == 1 && c(sim)[1] == "2") {
        x0[2] <- 1
    } else {
        x0[3] <- 1
    }

    loglik <- sum(x0*log(p))
    mu <- sum(p*log(p))
    var <- sum(p*(1-p)*log(p)^2)
    data.frame(L = loglik, mu = mu, var = var)
}

pedNoA_zscore_sim <- function(x,
                              idsx,
                              y = x,
                              idsy = idsx,
                              afreq, ## a list of SNPS
                              nindep = 0,
                              nsim = 100,
                              marker_list = NULL,
                              ncores = 1) {

    zs <- vector(mode = "numeric", length = nsim)
    pbapply::pbsapply(X = zs, cl = ncores, FUN = function(z) {
        dfs <- lapply(afreq, function(af) {
            .pedNoA_zscore_sim_onemarker(x,
                               idsx,
                               y,
                               idsy,
                               unlist(af),
                               nindep,
                               marker_list)
        })
        df <- do.call(rbind, dfs)
        sum(df$L - df$mu) / sqrt(sum(df$var))
    })
}

# load(file = "danish_snps.Rdat") ## danish_snps
# markers <- danish_snps
# names(markers) <- lapply(markers, function(x) unique(names(x)))
# markers <- lapply(markers, unname)
# af <- markers[[1]]
# library(paramlink)
# library(pedNoA)
# x <- nuclearPed(2); idsx <- 3:4
# y <- cousinsPed(1); idsy <- c(3, 8)
#
# out <- pedNoA_zscore_sim(x, idsx, y, idsy, afreq = markers, nsim = 50, ncores = 1)
# plot(out)
#
# ## make a t-test based on the mean
# mean(out)
# plot(out)
