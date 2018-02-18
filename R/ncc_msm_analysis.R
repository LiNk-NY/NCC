# Fitting NCC data using R's msm package with bootstrap
# on correlated data
#
# refer to the toturial for details
#
# Hongbin Zhang
# Marcel Ramos

library(msm)
library(dplyr)
library(abind)
library(readr)

## For installation, run the following:
# source("https://bioconductor.org/biocLite.R")
# BiocInstaller::biocLite("BiocParallel")
library(BiocParallel)

ncc <- read_csv("data/NCCstatus.csv")

txsit <- statetable.msm(STATUS, IDLOC, data=subset(ncc, STATUS!=99))
# save(txsit, file = "data/statetable.Rda")

## Total number of transition and consistent states
sum(statetable.msm(STATUS,IDLOC, data=subset(ncc, STATUS!=99)))

(treated <- statetable.msm(STATUS, IDLOC,
    data=subset(ncc, STATUS!= 99 & drug==1)))
(placebo <- statetable.msm(STATUS, IDLOC,
    data=subset(ncc, STATUS!= 99 & drug==0)))

Q <- rbind(
    c(0, 0.25,    0, 0.25),
    c(0,    0, 0.25, 0.25),
    c(0,    0,    0, 0.25),
    c(0,    0,    0,    0)
)

Q.ini <- crudeinits.msm(STATUS ~ MONTH, IDLOC, data = ncc, qmatrix = Q,
    censor = 99, censor.states=c(2,3,4))

###########################
# hazard ratio estimation #
###########################

fit0 <- msm(STATUS ~ MONTH, subject = IDLOC, data = ncc, censor = 99,
    censor.states = c(2,3,4), qmatrix = Q.ini, covariates = ~ drug)

hazard.msm(fit0)

# note, above fitting assumes independence among patients with multiple
# locations. To correct that, we use the bootstrap method

# get bootstrap based hazard 95% CI

# the bootstrap size (This will take about 4.3 minutes using 6 threads)
B <- 1000

# we get bootstrap from re-sampling on the patients
ids <- unique(ncc$ID)
N   <- length(ids)
split_ncc <- split(ncc, ncc$ID)

set.seed(123)
system.time({
hm.boot0 <- BiocParallel::bplapply(seq_len(B), function(x) {
    # it's important to specify replace=TRUE for bootstrapping
    bt.ids <- sample(ids, N, replace=TRUE)
    ncc.boot <- dplyr::bind_rows({
        idx <- vapply(as.character(bt.ids), function(x) {
            which(names(split_ncc) %in% x)
            }, integer(1L))
        tt <- split_ncc[idx]
        lapply(seq_along(tt), function(i, x) {
            x[[i]]$bootid <- paste(x[[i]]$IDLOC, as.character(i), sep = "-")
            x[[i]]
        }, x = tt)
    })
    fit <- msm(STATUS ~ MONTH, subject = bootid, data = ncc.boot, censor = 99,
        censor.states = c(2, 3, 4), qmatrix = Q.ini,
        covariates = ~ drug + Parenchymal)
    hm <- hazard.msm(fit)
    lapply(hm, function(x) x[, "HR"])
}, BPPARAM = MulticoreParam())
})

HRs <- vector("list", 3L)
for (i in seq_along(hm.boot0)) {
    HRs[[i]] <- do.call(rbind, hm.boot0[[i]])
}

HRarray <- sapply(HRs, function(x) {x}, simplify = "array")
HRresult <- t(apply(HRarray, 1:2, median))

HR.ci <- apply(HRarray, 1:2, function(x)
   quantile(x, probs = c(0.025, 0.975))
)
HR.ci <- aperm(HR.ci, c(2, 1, 3))

modelRes <- abind(HR = t(HRresult), HR.ci, along = 2L)
modelRes <- aperm(modelRes, c(3, 2, 1))

save(modelRes, file = "data/modelRes.Rda")
