# `pedNoA`: Number of shared alleles

## Description


 Calculates the probabilities of sharing `0,1,2,...` alleles between `ids` and `nindep` for the markers in `afreq` 


## Usage

```r
pedNoA(x, ids, afreq, nindep = 0, marker_list = NULL, verbose = TRUE,
  nc = 1)
```


## Arguments

Argument      |Description
------------- |----------------
```x```     |     The pedigree - a linkdat object from [`paramlink`](paramlink.html) .
```ids```     |     A vector of indicies from `x` .
```afreq```     |     A list of allele frequencies.
```nindep```     |     The number of independent contributors.
```marker_list```     |     A list of lists with known genotypes along with the respective locus index. See notes.
```verbose```     |     If TRUE, information about `x` will be printet.
```nc```     |     Number of cores to be used in parallelization,

## Seealso


 [`paramlink`](paramlink.html) , [`plot.pedNoA`](plot.pedNoA.html) 


## Note


 The procedure is much faster when `length(ids) == 2` . See paper. The `marker_list` argument
 is similar to that of `...` in the `marker` function in `paramlink` . If the the genotype
 of individual `1` is known to be heterozygotic with allele `2` and `3` at locus `2` , and individual `2` is heterozygote with `1` and `2` at locus `3` , the argument
 is `marker_list = list( list(2,3), list(1, c(2,2), 2, c(1,2)))` .


## Author


 Mads Lindskou, mads@math.aau.dk 


## Examples

```r 
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
 # with allele 3 at maker \code{2}  and grandmother being heterozygotic with allele 1 and 2 at marker \code{3}
 y <- cousinsPed(1)
 p3 <- c(.2,.2,.3,.3)
 plot(y)
 pNy <- pedNoA(y, c(5,8), list(p2, p3), marker_list = list(list(2,3), list(1, c(3,3), 2, c(1,2))))
 plot(pNy)
 
 ``` 

