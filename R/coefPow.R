.coefList <- function(x){
  l <- length(x)
  if(l==1) return(list(sign=1,pow=x,tail=FALSE))
  y <- x[-l]
  h <- x[l]
  res <- list(list(sign=1,pow=h,tail=.coefList(y)))
  for(i in 1:(l-1)){
    z <- y
    z[i] <- z[i]+h
    res <- c(res,list(list(sign=-1,pow=NULL,tail=.coefList(z))))
  }
  res
}

#' Coefficients and powers
#'
#' Determines the coefficients and powers needed to express
#' \deqn{\sum p_{i_1}^{\alpha_1} p_{i_2}^{\alpha_2} \cdots p_{i_n}^{\alpha_n}} in terms of \eqn{\sum p^{m}}.
#'
#' @param x A vector of powers
#' @author Poul Svante Eriksen, \email{svante@@math.aau.dk}
#'
#' @examples
#' # Computing \sum p_{i}^2 p_{j}^2 for i \neq j and p = c(.1,.2,.3,.4):
#' coefPow(c(2,2))
#' p = c(.1,.2,.3,.4)
#' -sum(p^4) + sum(p^2)^2
#'
#' @export
#'

coefPow <- function(x){
  # Make a print.coefPow and a sum.coefPow - maybe include m and ms (funny).
  res <- .coefList(x)
  ud <- c()
  oneX <- function(x){
    if(is.logical(x$tail$tail)) return(list(list(sign=x$sign*x$tail$sign,pow=c(x$pow,x$tail$pow))))
    res <- list()
    for(h in x$tail){
      h$sign <- h$sign*x$sign
      h$pow <- c(h$pow,x$pow)
      res <- c(res,oneX(h))
    }
    res
  }
  id <- function(x,y){
    if(length(x)!=length(y)) return(FALSE)
    sum((x-y)^2)==0
  }

  for(x in res) ud <- c(ud,oneX(x))
  ll <- c()
  for(i in 1:length(ud)) {
    ud[[i]]$pow <- sort(ud[[i]]$pow)
    ll <- c(ll,length(ud[[i]]$pow))
  }
  ud <- ud[order(ll)]
  L <- length(ud)
  if(L==1) return(ud)
  powL <- ud[1]
  for(i in 2:L){
    x <- ud[[i]]
    lp <- length(powL)
    for(j in 1:lp){
      if(id(x$pow,powL[[j]]$pow)) {
        powL[[j]]$sign <- powL[[j]]$sign+x$sign
        break
      }
      if(j==lp) powL <- c(powL,list(x))
    }
  }
  for(i in 1:length(powL)) names(powL[[i]])[1] <- "coef"
  powL
}
