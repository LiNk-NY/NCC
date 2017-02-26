## Read NCC SAV datasets

library(haven)

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
