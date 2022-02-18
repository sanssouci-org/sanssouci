# simulateGaussianNullsFromFactorModel
#
# Simulate test statistic null distribution in a Gaussian factor model
#
#
# @param m Number of tests
# @param n Number of replications of the simulation
# @param flavor The type of factor model to be simulated
# @param rho \code{1-rho} is the standard deviation of the noise
# @param cov should the covariance matrix of the model be returned ?
# @return A \code{m x n} \code{Matrix} simulated test statistics.  If
# \code{cov} is set to TRUE, the covariance matrix of the model may be
# accessed by \code{attr(Y, "Sigma")}.
# @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
#' @importFrom Matrix Matrix t
# @examples
#
# m <- 20
#
# ## independent
# Y <- sanssouci:::simulateGaussianNullsFromFactorModel(m, 
#     flavor = "independent")
#   image(Y)
#
# ## equi-correlated
# Y <- sanssouci:::simulateGaussianNullsFromFactorModel(m, 
#     flavor = "equi-correlated", rho = 0.2, cov = TRUE)
# S <- attr(Y, "Sigma")
# image(S)
#
# ## check equi-correlation:
# Y <- sanssouci:::simulateGaussianNullsFromFactorModel(m, 
#     n = 1237, flavor = "equi-correlated", rho = 0.2)
# covmat <- cov(t(Y))
# image(covmat)
# diag(covmat) <- NA
# dim(covmat) <- NULL
# summary(covmat)
#
# ## 3-factor model
# m <- 4*floor(m/4) ## make sure m/4 is an integer
# S3 <- sanssouci:::simulateGaussianNullsFromFactorModel(m, 
#     flavor = "3-factor", rho = 0.5)
# image(S3)
#
simulateGaussianNullsFromFactorModel <- function(
    m, n = 1, flavor = c("independent", "equi-correlated", "3-factor"), 
    rho = 0, cov = FALSE) {
    flavor <- match.arg(flavor)
    
    if (flavor=="independent") {
        h <- numeric(0)
        P <- Matrix(nrow=m, ncol=length(h))
        rho <- 0
        Y <- replicate(n, simulateFactorModelNullsFromSingularValuesAndLoadings(m, h=h, P=P, rho=rho)$Y)
    } else if (flavor=="equi-correlated") {
        h <- 1
        P <- Matrix(1, m, length(h))
        if (FALSE) { ## possibly slow for large m*n
            Y <- replicate(n, simulateFactorModelNullsFromSingularValuesAndLoadings(m, h=h, P=P, rho=rho)$Y)
        } else {
            Y <- simulateGaussianEquiCorrelatedNulls(m, n, rho)
        }
    } else if (flavor=="3-factor") {
        if (m != 4*floor(m/4)) {
            stop("Argument 'm' should be a multiplier of 4 for flavor '3-factor'")
        }
        h <- c(0.5, 0.3, 0.2)*m
        gamma1 <- rep(1,m)
        gamma2 <- rep(c(-1, 1), each=m/2)
        gamma3 <- rep(c(-1, 1), times=m/2)
        P <- Matrix(cbind(gamma1, gamma2, gamma3)/sqrt(m))
        Y <- replicate(n, simulateFactorModelNullsFromSingularValuesAndLoadings(m, h=h, P=P, rho=rho)$Y)
    }
    
    if (cov) {
        Sigma <- getFactorModelCovarianceMatrix(m, h=h, P=P, rho=rho)
        attr(Y, "Sigma") <- Sigma
    }
    Y
}

###########################################################################
# HISTORY:
# 2014-05-23
# o Now returning only 'Y' by default, and optionally 'Sigma' as an
# attribute.
# 2014-04-21
# o SPEEDUP for flavor 'equicorrelated'.
# 2013-04-26
# o Renamed to 'simulateGaussianNullsFromFactorModel'.
# 2013-03-29
# o Created.
###########################################################################

