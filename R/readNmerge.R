## Read NCC SAV datasets

## source("https://bioconductor.org/biocLite.R")
## BiocInstaller::biocLite("IRanges")

library(haven)
library(IRanges)
library(stringdist)

dataList <- list.files("data", pattern = "\\.sav$", full.names = TRUE)
names(dataList) <- gsub(".sav", "", basename(dataList), fixed = TRUE)

dataList <- lapply(dataList, read_spss)

for (i in seq_along(dataList)) {
## remove filter_$ variables
    if ("filter_$" %in% names(dataList[[i]])) {
        dataList[[i]] <- dataList[[i]][, -(which(
            names(dataList[[i]]) == "filter_$"))]
    }
## substitute bse to BSE in data
    if (any(grepl("bse", names(dataList[[i]])))) {
        names(dataList[[i]]) <- gsub("bse", "BSE", names(dataList[[i]]),
                                     fixed = TRUE)
    }
}

dataNames <- CharacterList(lapply(dataList, names))

timeStampStart <- LogicalList(lapply(dataNames, function(x) {
    grepl("^S[16Y]E|^BSE", x)
}))

timeVaryingNames <- dataNames[timeStampStart]

timeStampEnd <- dataNames[!timeStampStart]

timePoints <- vapply(timeVaryingNames,
                    function(nameList) {
                        TimeName <- sample(nameList, size = 1L)
                        substr(TimeName, 1, 3)
                    }, character(1L))

stopifnot(identical(names(dataNames), names(timePoints)))
stopifnot(identical(names(timeVaryingNames), names(timePoints)))

## Add TIME column to datasets
for (i in seq_along(dataList)) {
    dataList[[i]][["TIME"]] <- timePoints[[i]]
}

for (i in seq_along(timeVaryingNames)) {
    timeVaryingNames[[i]] <- gsub(timePoints[[i]], "", timeVaryingNames[[i]],
                                 ignore.case = TRUE)
}

endings <- IRanges::CharacterList(
    Baseline = paste(c("YN$", "Ca$"), collapse = "|"),
    M12Scan = paste(c("YNY$", "Y$", "YY$", "CY$", "M12$"), collapse = "|"),
    M1Scan = paste(c("YN1$", "1$", "Y1$", "C1$", "M1$"), collapse = "|"),
    M6Scan = paste(c("YN6$", "6$", "Y6$", "C6$", "M6$"), collapse = "|"))

newNames <- mapply(function(patterns, varnames) {
    gsub(patterns, "", varnames, ignore.case = TRUE)
}, patterns = endings, varnames = timeStampEnd)

outersect <- function(x, y) {
        sort(c(setdiff(x, y), setdiff(y, x)))
}

unMatched0 <- Reduce(outersect, newNames)

## Look at this list (includes coded variables, e.g., '2A2B1')
do.call(cbind, split(unMatched0, c(TRUE, FALSE)))

unmatched <- unMatched0[grepl("^[A-Z]", unMatched0, ignore.case = TRUE)]

lvDist <- stringdistmatrix(tolower(unmatched), method = "lv")
shortDist <- as.matrix(lvDist) == 1L | as.matrix(lvDist) == 2L
validPairs <- reshape2::melt(shortDist)
matIndex <- validPairs[validPairs[[3]], 1:2]

first <- apply(matIndex, 1, function(x) x[1] < x[2])
closeIdx <- matIndex[first,]
closeNames <- cbind.data.frame(first = unmatched[closeIdx[[1]]],
                                 second = unmatched[closeIdx[[2]]],
                                 stringsAsFactors = FALSE)
closeNames[, "ncharfirst"] <- nchar(closeNames[[1]])
closeNames[, "ncharsecond"] <- nchar(closeNames[[2]])
## Create a column indicating longer word
closeNames[, "longer"] <- apply(X = closeNames[, c("ncharfirst", "ncharsecond")],
                                  MARGIN = 1,
                                  FUN = function(x) {
                                      if (x[1] == x[2])
                                          NA_integer_
                                      else
                                          which.max(x)
                                  })
## Calculate distances for pairs
closeNames[, "levDist"] <- stringdist(closeNames[["first"]],
                                        closeNames[["second"]], method = "lv")

## Make NA words with a high distance
closeNames[which(closeNames[, "levDist"] == 3), "longer"] <- NA_integer_

closeNames <- closeNames[complete.cases(closeNames), ]
