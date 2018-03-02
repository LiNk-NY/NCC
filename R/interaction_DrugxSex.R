## Interaction analysis DRUG by SEX
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
library(haven)
## For installation, run the following:
# source("https://bioconductor.org/biocLite.R")
# BiocInstaller::biocLite("BiocParallel")
library(BiocParallel)

ncc <- read_csv("data/NCCstatus.csv")

Q <- rbind(
    c(0, 0.25,    0, 0.25),
    c(0,    0, 0.25, 0.25),
    c(0,    0,    0, 0.25),
    c(0,    0,    0,    0)
)

Q.ini <- crudeinits.msm(STATUS ~ MONTH, IDLOC, data = ncc, qmatrix = Q,
    censor = 99, censor.states=c(2,3,4))

#################################
#                               #
#  subgroup analysis            #
#                               #
#  assuming a group variable    #
#  is coded 0/1                 #
#                               #
#  refer to tutorial for details#
#################################

###########################
# hazard ratio estimation #
###########################

# DRUG X SEX -------------------------------------------------------------

# readin the patient info data
# file.rename("data/Adult & Pediatric Patient Information Form August162004.sav", "data/patinfo.sav")
patinfo <- read_spss("data/patinfo.sav")

age.raw <- as.numeric((patinfo$PIADATE - patinfo$PIA3)/365)

age <- ifelse(age.raw > 40, 1, 0)

male <- ifelse(patinfo$PIA2 == 1, 1, 0);

sub.patinfo <- data.frame(cbind(ID = patinfo$ID, age, male))

ncc <- merge(ncc, sub.patinfo, by="ID")

# note, above fitting assumes independence among patients with multiple
# locations. To correct that, we use the bootstrap method

# get bootstrap based hazard 95% CI

# the bootstrap size (This will take about 4.3 minutes using 6 threads)
B <- 1000

# we get bootstrap from re-sampling on the patients
ids <- unique(ncc$ID)
N   <- length(ids)
split_ncc <- split(ncc, ncc$ID)

set.seed(213)
system.time({
hm.boot <- BiocParallel::bplapply(seq_len(B), function(x) {
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
        censor.states = c(2,3,4), qmatrix = Q.ini,
        covariates = ~ drug * male)

    hm <- hazard.msm(fit)
    beta.mat <- vapply(hm, function(x) log(x[, "HR"]), numeric(5L))

    hazardA <- exp(rowSums(beta.mat))
    hazardB <- exp(0 + beta.mat[, "male"] + 0)

    HR.male1 <- hazardA/hazardB

    hazardC <- exp(beta.mat[, "drug"] + 0 + 0)
    hazardD <- exp(0 + 0 + 0)

    HR.male0 <- hazardC/hazardD

    cbind(HR.male0, HR.male1)

    }, BPPARAM = MulticoreParam(workers = 20, RNGseed = 213))
})

HRarray <- sapply(hm.boot, function(x) {x}, simplify = "array")
HRresult <- t(apply(HRarray, 1:2, median))

HR.ci <- apply(HRarray, 1:2, function(x)
   quantile(x, probs = c(0.025, 0.975))
)
HR.ci <- aperm(HR.ci, c(2, 1, 3))

drugmaleint <- abind(HR = t(HRresult), HR.ci, along = 2L)

save(drugmaleint, file = "data/interactions/drugmaleint.Rda")
