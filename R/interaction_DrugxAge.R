## Interaction analysis DRUG by AGE
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

# Drug X Age --------------------------------------------------------------

# readin the patient info data
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

set.seed(321)
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
        covariates = ~ drug * age)

    hm <- hazard.msm(fit)
    lapply(hm, function(x) x[, "HR"])
}, BPPARAM = MulticoreParam())
})

HRs <- vector("list", 3L)
for (i in seq_along(hm.boot)) {
    HRs[[i]] <- do.call(rbind, hm.boot[[i]])
}

HRarray <- sapply(HRs, function(x) {x}, simplify = "array")
HRresult <- t(apply(HRarray, 1:2, median))

HR.ci <- apply(HRarray, 1:2, function(x)
   quantile(x, probs = c(0.025, 0.975))
)
HR.ci <- aperm(HR.ci, c(2, 1, 3))

drugageint <- abind(HR = t(HRresult), HR.ci, along = 2L)
drugageint <- aperm(drugageint, c(3, 2, 1))

save(drugageint, file = "data/interactions/drugageint.Rda")

# hazard for drug.x = 1, and loc = 1
hazardA <- exp(rowSums(log(drugageint[, "HR", ])))
# note I did not include the baseline beta since it will
#be cancelled anyway

# hazard for drug.x = 0, and loc = 1
hazardB <- exp(0 + log(drugageint[, "HR", "age"]) + 0)

# hazard ratio (treatment vs. placebo ) for loc = 1
HR.age1 <- hazardA/hazardB

# hazard for drug.x = 1, and loc = 0
hazardC <- exp(drugageint[, "HR", "drug"] + 0 + 0)

# hazard for drug.x = 0, and loc = 0
hazardD <- exp(0 + 0 + 0)

# hazard ratio (treatment vs. placebo ) for group = 1
HR.age0 <- hazardC/hazardD

# put results together
(interaction.age <- t(rbind(HR.age0, HR.age1)))
