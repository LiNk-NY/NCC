#
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

## For installation, run the following:
# source("https://bioconductor.org/biocLite.R")
# BiocInstaller::biocLite("BiocParallel")
library(BiocParallel)

fullData <- read.csv("data/NCClong.csv", header=TRUE, stringsAsFactors = FALSE)
ncc <- fullData %>% select(ID:STATUS, drug) %>% arrange(ID, LocCode, MONTH)

cleanLocationDF <- function(datframe) {
    stat4 <- datframe[["STATUS"]] == 4L
    if (any(stat4, na.rm = TRUE)) {
    datframe <- datframe[seq_len(which.max(stat4)), ]
    datframe <- datframe[!is.na(datframe[["STATUS"]]), ]
    }
    nas <- is.na(datframe[["STATUS"]])
    if (any(nas)) {
        restNA <- identical(which(nas), seq(which.max(nas), length(nas)))
        if (restNA) {
            datframe[["STATUS"]][which.max(nas)] <- 99L
        } else {
            if (sum(nas) > 1L) {
            datframe[["STATUS"]][6L - which.max(rev(nas))] <- 99L
            }
        }
    }
    datframe[!is.na(datframe[["STATUS"]]), ]
}

ncc <- dplyr::bind_rows(lapply(split(ncc, ncc[["IDLOC"]]), cleanLocationDF))
ncc$Parenchymal <- dplyr::recode(ncc$Tissue, "Parenchymal" = "Yes", .default = "No")
## Recode non-parenchymal regions as "No"
ncc$Parenchymal[is.na(ncc$Parenchymal)] <- "No"

## minimum number of recorded visits
group_by(ncc, ID) %>% summarise(nMonth = n_distinct(MONTH)) %>% select(nMonth) %>% min
group_by(ncc, ID) %>% summarize(dislv = sum(STATUS == 4), dlogic = sum(dislv != 0L)) %>%
    summarize(ndislv = sum(dlogic))
round(101/117 * 100, 2)
filter(ncc, MONTH == 0L) %>% group_by(ID) %>%
    summarize(activ = sum(STATUS == 1L)) %>%
    summarize(n1cysts = sum(activ >= 1L))
round(62/117 * 100, 2)

txsit <- statetable.msm(STATUS, IDLOC, data=subset(ncc, STATUS!=99))

sum(statetable.msm(STATUS,IDLOC, data=subset(ncc, STATUS!=99)))
(treated <- statetable.msm(STATUS, IDLOC, data=subset(ncc, STATUS!= 99 & drug==1)))
(placebo <- statetable.msm(STATUS, IDLOC, data=subset(ncc, STATUS!= 99 & drug==0)))

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
        censor.states = c(2, 3, 4), qmatrix = Q.ini,
        covariates = ~ drug + Parenchymal)
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

modelRes <- abind(HR = t(HRresult), HR.ci, along = 2L)
modelRes <- aperm(modelRes, c(3, 2, 1))

# save(modelRes, file = "data/modelRes.Rda")
# load("data/modelRes.Rda")

resultTable1 <- do.call(cbind,
lapply(c(Albendazole = 1, Parenchymal = 2), function(g)
    apply(round(modelRes[, , g], 2), 1L, function(x)
        paste(x[1], paste0("(", paste(x[2:3], collapse = " - "), ")")))
    )
)

stateChanges <- gsub("State 1", "Active", rownames(resultTable1), fixed = TRUE) %>%
    gsub("State 2", "Transitional", ., fixed = TRUE) %>%
    gsub("State 3", "Calcified", ., fixed = TRUE) %>%
    gsub("State 4", "Dissolved", ., fixed = TRUE)

rownames(resultTable1) <- stateChanges

countstt <- c(txsit[1, 2], txsit[1, 4], txsit[2, 3], txsit[2, 4], txsit[3, 4])
"n (%)" <- paste(countstt, paste0("(", round(countstt/sum(txsit) * 100, 1), ")"))

resultTable1 <- cbind(`n (%)`, resultTable1)

# write.csv(resultTable1, "data/resultTable1.csv")

#############################
# plot the survival curve   #
#############################

par(mfrow=c(1,2))
plot(fit0, covariate = list(1), show.legend = FALSE,
        xlab="Time after trial start (months)", las=1)
legend(x = 8, y = 1,
    legend = paste("From", c("Active", "Transitional", "Calcified")),
    cex = 0.9, lty = seq(4-1), col = rainbow(3), lwd = 1)
title("Albendazole Group")

plot(fit0, covariate = list(0), show.legend = FALSE,
        xlab="Time after trial start (months)", las=1)
legend(x = 8, y = 1,
    legend = paste("From", c("Active", "Transitional", "Calcified")),
    cex = 0.9, lty = seq(4-1), col = rainbow(3), lwd = 1)
title("Placebo Group")

