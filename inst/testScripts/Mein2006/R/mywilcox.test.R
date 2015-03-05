## copied from stats:::wilcox.test.default
## speeded up by removing redundant tests and calls when calling the function multiple times
mywilcox.test <- function (x, y=NULL, alternative=c("two.sided", "less", "greater"), 
                           mu=0, paired=FALSE, exact=NULL, correct=TRUE, conf.int=FALSE, 
                           conf.level=0.95,
                           p.value=TRUE,
                           ## calculate p-value ? (if not, it is faked using -abs(test statistic)...)
                           bypassTiesMethod=FALSE,
                           ## bypass argument 'ties.method' of the 'rank' function?
                           ...) {
  ## stopifnot(!is.null(y))
  ## stopifnot(!conf.int)
  DNAME <- "x and y" ## fake but fast
  alternative <- match.arg(alternative)
  METHOD <- "Wilcoxon rank sum test"

  ## bypass argument 'ties.method' of the 'rank' function, as the default value "average" is forced by wilcox.test
  if (!bypassTiesMethod) {
    r <- rank(c(x - mu, y))
  } else {
    xx <- c(x - mu, y)
    r <- .Internal(rank(xx, 
                        length(xx), "average"))  # forcing ties.method="average"
  }
  
  n.x <- as.double(length(x))
  n.y <- as.double(length(y))
  if (is.null(exact)) 
      exact <- (n.x < 50) && (n.y < 50)
  STATISTIC <- c(W = sum(r[seq_along(x)]) - n.x * (n.x + 
                     1)/2)
  if (p.value) {
    TIES <- (length(r) != length(unique(r)))
    if (exact && !TIES) {
      PVAL <- switch(alternative, two.sided = {
        p <- if (STATISTIC > (n.x * n.y/2)) pwilcox(STATISTIC - 
                                                    1, n.x, n.y, lower.tail = FALSE) else pwilcox(STATISTIC, 
                                                                   n.x, n.y)
      min(2 * p, 1)
      }, greater = {
        pwilcox(STATISTIC - 1, n.x, n.y, lower.tail = FALSE)
      }, less = pwilcox(STATISTIC, n.x, n.y))
      if (conf.int) {
        alpha <- 1 - conf.level
        diffs <- sort(outer(x, y, "-"))
        cint <- switch(alternative, two.sided = {
          qu <- qwilcox(alpha/2, n.x, n.y)
          if (qu == 0) qu <- 1
          ql <- n.x * n.y - qu
          achieved.alpha <- 2 * pwilcox(trunc(qu) - 1, 
                                        n.x, n.y)
          c(diffs[qu], diffs[ql + 1])
        }, greater = {
          qu <- qwilcox(alpha, n.x, n.y)
          if (qu == 0) qu <- 1
          achieved.alpha <- pwilcox(trunc(qu) - 1, n.x, 
                                    n.y)
        c(diffs[qu], +Inf)
        }, less = {
          qu <- qwilcox(alpha, n.x, n.y)
          if (qu == 0) qu <- 1
          ql <- n.x * n.y - qu
          achieved.alpha <- pwilcox(trunc(qu) - 1, n.x, 
                                    n.y)
          c(-Inf, diffs[ql + 1])
        })
        if (achieved.alpha - alpha > alpha/2) {
          warning("Requested conf.level not achievable")
          conf.level <- 1 - achieved.alpha
        }
        attr(cint, "conf.level") <- conf.level
      ESTIMATE <- c(`difference in location` = median(diffs))
      }
    }
    else {
      NTIES <- table(r)
      z <- STATISTIC - n.x * n.y/2
      SIGMA <- sqrt((n.x * n.y/12) * ((n.x + n.y + 1) - 
                                    sum(NTIES^3 - NTIES)/((n.x + n.y) * (n.x + n.y - 
                                                                         1))))
    if (correct) {
      CORRECTION <- switch(alternative, two.sided = sign(z) * 
                           0.5, greater = 0.5, less = -0.5)
      METHOD <- paste(METHOD, "with continuity correction")
    }
      z <- (z - CORRECTION)/SIGMA
      PVAL <- switch(alternative,
                     less = pnorm(z),
                     greater = pnorm(z, lower.tail = FALSE),
                     two.sided = 2 * min(pnorm(z), pnorm(z, lower.tail = FALSE)))
      if (conf.int) {
        alpha <- 1 - conf.level
        mumin <- min(x) - max(y)
        mumax <- max(x) - min(y)
        wdiff <- function(d, zq) {
          dr <- rank(c(x - d, y))
        NTIES.CI <- table(dr)
          dz <- (sum(dr[seq_along(x)]) - n.x * (n.x + 
                                                1)/2 - n.x * n.y/2)
          CORRECTION.CI <- if (correct) {
            switch(alternative, two.sided = sign(dz) * 
                   0.5, greater = 0.5, less = -0.5)
          }
          else 0
          SIGMA.CI <- sqrt((n.x * n.y/12) * ((n.x + n.y + 
                                              1) - sum(NTIES.CI^3 - NTIES.CI)/((n.x + n.y) * 
                                                                               (n.x + n.y - 1))))
          if (SIGMA.CI == 0) 
              stop("cannot compute confidence interval when all observations are tied", 
                   call. = FALSE)
          (dz - CORRECTION.CI)/SIGMA.CI - zq
        }
        root <- function(zq) {
          f.lower <- wdiff(mumin, zq)
          if (f.lower <= 0) 
              return(mumin)
          f.upper <- wdiff(mumax, zq)
          if (f.upper >= 0) 
              return(mumax)
          uniroot(wdiff, c(mumin, mumax), f.lower = f.lower, 
                  f.upper = f.upper, tol = 1e-04, zq = zq)$root
        }
        cint <- switch(alternative, two.sided = {
          l <- root(zq = qnorm(alpha/2, lower.tail = FALSE))
          u <- root(zq = qnorm(alpha/2))
          c(l, u)
        }, greater = {
          l <- root(zq = qnorm(alpha, lower.tail = FALSE))
          c(l, +Inf)
        }, less = {
          u <- root(zq = qnorm(alpha))
          c(-Inf, u)
        })
        attr(cint, "conf.level") <- conf.level
        correct <- FALSE
        ESTIMATE <- c(`difference in location` = uniroot(wdiff, 
                          c(mumin, mumax), tol = 1e-04, zq = 0)$root)
      }
      if (exact && TIES) {
        warning("cannot compute exact p-value with ties")
        if (conf.int) 
            warning("cannot compute exact confidence intervals with ties")
      }
    }
  } else {  ## just symmetrize the test statistic wrt its mean value
    med <- length(x) * length(y)/2 
    PVAL <- -abs(STATISTIC-med)
  }
  names(mu) <- if (paired || !is.null(y)) 
      "location shift"
  else "location"
  RVAL <- list(statistic = STATISTIC, parameter = NULL, p.value = as.numeric(PVAL), 
               null.value = mu, alternative = alternative, method = METHOD, 
               data.name = DNAME)
  if (conf.int) 
      RVAL <- c(RVAL, list(conf.int = cint, estimate = ESTIMATE))
  class(RVAL) <- "htest"
  RVAL
}
