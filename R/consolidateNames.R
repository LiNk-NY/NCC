## Script unused to consolidate names in all datasets rather than
## take only intersected names

commentIdx <- unique(vapply(cleanNames,
    function(x) { which(x == "COMME") }, numeric(1L)))

aftComm <- endoapply(timeNames, function(x) x[commentIdx:length(x)])
timeNames <- timeNames[!timeNames %in% aftComm]

stopifnot(identical(names(endings), names(notTimeNames)))
similarNames <- S4Vectors::mendoapply(function(patterns, varnames) {
    gsub(patterns, "", varnames, ignore.case = FALSE) ## Cases correct
}, patterns = endings, varnames = notTimeNames)

keepName <- Reduce(intersect, similarNames)
simNames <- Filter(length, similarNames[!similarNames %in% keepName])

inAllData <- Reduce(intersect, simNames)
differs <- simNames[!simNames %in% inAllData]

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

dataNames[!haveTimeStamp] <- reNewNames
