# Script for generating tables and figures

# load("data/statetable.Rda")
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

offtx <- txsit
diag(offtx) <- 0L
offdgsum <- sum(offtx)

countstt <- c(txsit[1, 2], txsit[1, 4], txsit[2, 3], txsit[2, 4], txsit[3, 4])
"No. of Events (%)" <- paste(countstt, paste0("(", round(countstt/offdgsum * 100, 1), ")"))

resultTable1 <- cbind(`No. of Events (%)`, resultTable1)

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
