## general setup
library("sansSouci")
library("howmany.pn")  ## to be able to uncomment the magic line in the implementation of howmany
library("R.utils")

verbose <- TRUE

pkg <- "R.menu"
if (!require(pkg, character.only=TRUE)) {
    install.packages(pkg, repos="http://R-Forge.R-project.org")
}
library(pkg, character.only=TRUE)
## sourceDirectory("R")
dirp <- system.file("testScripts/Mein2006/R", package="sansSouci")
if (!file.exists(dirp)) {
    stop()
}
sourceDirectory(dirp)

## parameters
## m <- textMenu(c(10, 100, 1000), title="Number of hypotheses to test", value=TRUE)
m <- 1000
pi0 <- textMenu(c(1, 0.99, 0.7, 0.6), title="Fraction of true null hypotheses", value=TRUE)
nbSimu <- textMenu(c(10, 100, 1000, 10000), title="Number of simulations", value=TRUE)
B <- textMenu(c(1000, 5e3, 1e4), title="Number of permutations", value=TRUE)
flavorH <-textMenu(c("howmany,2006", "howmany,2013", "sansSouci"), title="Flavor of the analysis", value=TRUE)
if (flavorH=="sansSouci") {
    stepDown <-textMenu(c(TRUE, FALSE), title="step down ?", value=TRUE)
    sv <- packageDescription("sansSouci")$Version
    if (stepDown) {
        ftag <- paste("sansSouci,SD,",  sv, sep="")
    } else {
        ftag <- paste("sansSouci,",  sv, sep="")
    }
} else {
    ftag <- flavorH
}

sort <- TRUE ## compare 'Q' with *sorted* p-values ?  Should be the case, but there seems to be a bug in Meinshausen's code...

p <- 0.5  ## binomial proportion for the two classes
alpha <- 0.05

ns <- selectMenu(c(20, 40, 60, 80, 100), title="Number of observations", selected=5)
ns <- 100

rhos <- selectMenu(c(0, 0.2, 0.4), title="Correlation level", selected=c(1))

SNR <- textMenu(c(0.1, 0.2, 0.5, 1), title="Signal to noise ratio", value=TRUE)

maxcores <- parallel::detectCores()
mc.coress <- sort(unique(c(1:min(maxcores, 4), round(maxcores/(1:4)))))
mc.cores <- textMenu(mc.coress, title=sprintf("Number of cores to be used (out of %s)", maxcores), value=TRUE)


## paths
## study <- ifelse(sort, sprintf("Mein%s,sorted", flavorH), sprintf("Mein%s", flavorH))
## study <- ifelse(sort, sprintf("Mein%s,sorted", flavorH), sprintf("Mein%s", flavorH))
study <- sprintf("Mein2006,%s", ftag)

path <- file.path("resData", study)
path <- Arguments$getWritablePath(path)

figPath <- file.path("fig", study, "fig1")
figPath <- Arguments$getWritablePath(figPath)

##
runForce <- TRUE
