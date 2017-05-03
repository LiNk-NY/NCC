## Read NCC SAV datasets

## source("https://bioconductor.org/biocLite.R")
## BiocInstaller::biocLite("IRanges")

suppressPackageStartupMessages({
    library(haven)
    library(dplyr)
    library(IRanges)
    library(stringdist)
    library(purrr)
})

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

drugVar <- dataList[[1]][, c("ID", "drug", "drug_gro")]

## Remove extra variables
dataList[[1]] <- dataList[[1]][, -which(names(dataList[[1]]) %in% c("drug", "drug_gro"))]

## Convenience function for getting the outersect
outersect <- function(x, y) {
    sort(c(setdiff(x, y), setdiff(y, x)))
}

dataNames <- CharacterList(lapply(dataList, names))

timeStampStart <- LogicalList(lapply(dataNames, function(x) {
    grepl("^S[16Y]E|^BSE", x)
}))

timeVaryingNames <- dataNames[timeStampStart]

timeStampEnd <- dataNames[!timeStampStart]

## Sample time stamped names and take the first 3 characters (time indicator)
timePoints <- vapply(timeVaryingNames,
                    function(nameList) {
                        TimeName <- sample(nameList, size = 1L)
                        substr(TimeName, 1, 3)
                    }, character(1L))

stopifnot(identical(names(dataNames), names(timePoints)))
stopifnot(identical(names(timeVaryingNames), names(timePoints)))

## Remove time indicators from variable names
for (i in seq_along(timeVaryingNames)) {
    timeVaryingNames[[i]] <- gsub(timePoints[[i]], "", timeVaryingNames[[i]],
                                 ignore.case = TRUE)
}

## Check all vars match across datasets
stopifnot(!length(Reduce(outersect, timeVaryingNames)))

endings <- IRanges::CharacterList(
    Baseline = paste(c("YN$", "y$", "Ca$"), collapse = "|"),
    M12Scan = paste(c("YNY$", "Y$", "YY$", "CY$", "M12$"), collapse = "|"),
    M1Scan = paste(c("YN1$", "1$", "Y1$", "C1$", "M1$"), collapse = "|"),
    M6Scan = paste(c("YN6$", "6$", "Y6$", "C6$", "M6$"), collapse = "|"))

newNames <- S4Vectors::mendoapply(function(patterns, varnames) {
    gsub(patterns, "", varnames, ignore.case = FALSE) ## Cases correct
}, patterns = endings, varnames = timeStampEnd)

unmatched <- Reduce(outersect, newNames)

## Look at this list (includes coded variables, e.g., '2A2B1')
do.call(cbind, split(unmatched, c(TRUE, FALSE)))

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

shorter <- closeNames[, "longer"] %% 2
shorter[shorter == 0L] <- -1
closeNames[, "shorter"] <- closeNames[["longer"]] + shorter

## Fix rownames
closeNames <- data.frame(closeNames, row.names = seq_len(nrow(closeNames)),
                         stringsAsFactors = FALSE)

## Create a corrections data.frame for mapping short to long variable names
corrections <- cbind.data.frame(
    short = vapply(seq_len(nrow(closeNames)),
                   function(i) {
                       closeNames[i, closeNames[["shorter"]][i]]
                   }, character(1L)),
    long = vapply(seq_len(nrow(closeNames)),
                  function(i) {
                      closeNames[i, closeNames[["longer"]][i]]
                  }, character(1L)),
    stringsAsFactors = FALSE)

## Manual entry
corrections <- rbind(corrections,
                     data.frame(short = "totcyst", long = "TotCyst"))

LongNewNames <- cbind.data.frame(varName = unlist(newNames),
                                 group = rep(names(newNames),
                                             lengths(newNames)),
                                 row.names = NULL,
                                 stringsAsFactors = FALSE)

isShort <- vapply(corrections$short, `==`,
                  logical(nrow(LongNewNames)),
                  LongNewNames$varName)
stopifnot(!any(rowSums(isShort) > 1))

replaceIdx <- apply(isShort, 1, which.max)
replaceIdx[replaceIdx == 1L] <- 0L
truePos <- which(isShort[, 1])
replaceIdx[truePos] <- 1L

NonZero <- replaceIdx != 0L
replaceVec <- corrections$long[replaceIdx[NonZero]]

stopifnot(identical(length(replaceVec), sum(NonZero)))

LongNewNames[NonZero, 1] <- corrections$long[replaceIdx[NonZero]]
reNewNames <- IRanges::splitAsList(LongNewNames[[1]], LongNewNames[[2]])

## Check all vars match across datasets
stopifnot(!length(Reduce(outersect, reNewNames)))

dataNames[!timeStampStart] <- reNewNames
dataNames[timeStampStart] <- timeVaryingNames

stopifnot(!length(Reduce(outersect, dataNames)))

newDataList <- lapply(seq_along(dataList), function(i, dataset) {
    names(dataset[[i]]) <- dataNames[[i]]
    dataset[[i]]
}, dataset = dataList)

names(newDataList) <- names(dataList)

timeNumeric <- c(BSE = 0L, SYE = 12L, S1E = 1L, S6E = 6L)

## Check any variable with the name month already in data
stopifnot(!any(tolower(names(newDataList[[1]])) == "month"))

## Add TIME column to datasets
for (i in seq_along(newDataList)) {
    newDataList[[i]] <- tibble::add_column(newDataList[[i]],
                                           MONTH = rep(timeNumeric[[i]],
                                                       nrow(newDataList[[i]])),
                                           .after = 1)
}

NCCdata <- dplyr::bind_rows(newDataList)

NCCdata <- dplyr::arrange(NCCdata, ID, MONTH)

## Find which IDs have less than 4 measurements
group_by(NCCdata, ID) %>% summarize(N = n())

# Impute missing timepoints
fullNCC <- split(NCCdata, NCCdata[["ID"]]) %>%
    map(~ suppressMessages(right_join(.x,
                     tibble(ID = rep(unique(.x[["ID"]]), length(timeNumeric)),
                            MONTH = sort(timeNumeric)))))

fullNCC <- dplyr::bind_rows(fullNCC)

## TODO: Use proper format and remove `sapply`
validC <- split(aa, aa$ID)[sapply(lapply(split(aa, aa$ID), function(x) {apply(x, 1, .checkSum)}), function(x) !any(x, na.rm = TRUE))]
