## Read NCC SAV datasets

## source("https://bioconductor.org/biocLite.R")
## BiocInstaller::biocLite("IRanges")

library(haven)
library(IRanges)
library(stringdist)

dataList <- list.files("data", pattern = "\\.sav$", full.names = TRUE)
names(dataList) <- gsub(".sav", "", basename(dataList), fixed = TRUE)

dataList <- lapply(dataList, read_spss)

timeVaryingNames <- lapply(dataList,
                           function(x)
                               grep("^S[16Y]E|^BSE", x = names(x), value = TRUE))

TimePoint <- vapply(timeVaryingNames,
                    function(nameList) {
                        TimeName <- sample(nameList, size = 1L)
                        substr(TimeName, 1, 3)
                    }, character(1L))

stopifnot(identical(names(dataList), names(TimePoint)))

for (i in seq_along(dataList)) {
    dataList[[i]][["TIME"]] <- TimePoint[[i]]
    if ("filter_$" %in% names(dataList[[i]])) {
        dataList[[i]] <- dataList[[i]][, -(which(
            names(dataList[[i]]) == "filter_$"))]
    }
    names(dataList[[i]]) <- gsub(TimePoint[[i]], "", names(dataList[[i]]),
                                 ignore.case = TRUE)
}

endings <- IRanges::CharacterList(
    Baseline = paste(c("YN$", "Y$", "YNY$"), collapse = "|"),
    M12Scan = paste(c("YNY$", "Y$", "YY$", "CY$", "M12$"), collapse = "|"),
    M1Scan = paste(c("YN1$", "1$", "Y1$", "C1$", "M1$"), collapse = "|"),
    M6Scan = paste(c("YN6$", "6$", "Y6$", "C6$", "M6$"), collapse = "|"))

newNames <- mapply(function(patterns, varnames) {
    gsub(patterns, "", varnames, ignore.case = TRUE)
}, patterns = endings, varnames = lapply(dataList, names))

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
SimilarNames <- cbind.data.frame(first = unmatched[closeIdx[[1]]],
                                 second = unmatched[closeIdx[[2]]],
                                 stringsAsFactors = FALSE)
SimilarNames[, "ncharfirst"] <- nchar(SimilarNames[[1]])
SimilarNames[, "ncharsecond"] <- nchar(SimilarNames[[2]])

## Decision data.frame from word pairs
cbind.data.frame(SimilarNames, standard =
ifelse(SimilarNames[["ncharfirst"]] > SimilarNames[["ncharsecond"]] &
           nchar(SimilarNames[["first"]]) != nchar(SimilarNames[["second"]]),
       SimilarNames[["first"]],
       ifelse(nchar(SimilarNames[["first"]]) == nchar(SimilarNames[["second"]]),
           NA_character_, SimilarNames[["second"]])),
stringsAsFactors = FALSE)
