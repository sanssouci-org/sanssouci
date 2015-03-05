k.max  <- structure(function(
### k^{th} maximum of a numeric vector
                   x,
### A numeric vector
                   k
### An integer vector of test statistics of length m
                   ) {
  x.len <- length(x)
  if (k > x.len) {
    result <- -Inf;
    warning("no non-missing arguments to max; returning -Inf");
  } else if (1==k) {
    result <- max(x)
  } else {
    x.sort <- sort(x)
    result <- x.sort[x.len-k+1]
  }
  result
### the k^{th} maximum of x
}, ex=function(){
  x <- rnorm(100)
  k.max(x, 10);
  sort(x)[100-10+1]
})

############################################################################
## HISTORY:
## 2012-10-08
## o Created from Romano & Wolf's implementation of kFWER control.
############################################################################

