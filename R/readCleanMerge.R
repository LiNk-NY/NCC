## Cleaning and merging datasets
## NOTE. run R/download.R to download data

## Load dependencies
source("R/loadPackages.R")

## Load helper dataset
source("R/helperData.R")

## Read NCC SAV datasets
dataList <- list.files("data", pattern = "^[BM].*\\.sav$", full.names = TRUE)
names(dataList) <- gsub(".sav", "", basename(dataList), fixed = TRUE)

dataList <- lapply(dataList, read_spss)


# Cleaning column names ---------------------------------------------------

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
# readr::write_csv(drugVar, "data/drugVars.csv")

## Remove extra variables
dataList[[1]] <-
    dataList[[1]][, -which(names(dataList[[1]]) %in% c("drug", "drug_gro"))]

## Convenience function for getting the outersect
outersect <- function(x, y) {
    sort(c(setdiff(x, y), setdiff(y, x)))
}

dataNames <- CharacterList(lapply(dataList, names))

haveTimeStamp <- LogicalList(lapply(dataNames, function(x) {
    grepl("^S[16Y]E|^BSE|^SLE", x) & !grepl("COMME", x, fixed = TRUE)
}))

timeNames <- dataNames[haveTimeStamp]

notTimeNames <- dataNames[!haveTimeStamp]

## Sample time stamped names and take the first 3 characters (time indicator)
timePoints <- vapply(timeNames,
                    function(nameList) {
                        TimeName <- sample(nameList, size = 1L)
                        substr(TimeName, 1, 3)
                    }, character(1L))

stopifnot(identical(names(dataNames), names(timePoints)))
stopifnot(identical(names(timeNames), names(timePoints)))

cleanNames <- CharacterList(vector("list", 5L))
names(cleanNames) <- names(timeNames)

## Remove time indicators from variable names
for (i in seq_along(timeNames)) {
    cleanNames[[i]] <- gsub(timePoints[[i]], "", timeNames[[i]],
                                 ignore.case = TRUE)
}

dataNames[haveTimeStamp] <- cleanNames

intNames <- Reduce(intersect, cleanNames)

## Take only names in all datasets
logicalSub <- dataNames %in% intNames

# Loading data ------------------------------------------------------------

## First, rename datasets with dataNames elements

newDataList <- lapply(seq_along(dataList), function(i, dataset) {
    names(dataset[[i]]) <- dataNames[[i]]
    dataset[[i]]
}, dataset = dataList)

names(newDataList) <- names(dataList)

## Second, Subset datasets with logicalSub
newDataList <- Map(function(x, y) { x[, y] }, x = newDataList, y = logicalSub)

timeNumeric <- c(BSE = 0L, SYE = 12L, S1E = 1L, S6E = 6L, SLE = 24L)

## Check any variable with the name month already in data
stopifnot(!any(tolower(names(newDataList[[1]])) == "month"))

## Add TIME column to datasets
for (i in seq_along(newDataList)) {
    newDataList[[i]] <- tibble::add_column(newDataList[[i]],
                                           MONTH = rep(timeNumeric[[i]],
                                                       nrow(newDataList[[i]])),
                                           .after = 1)
}
## Bind all time points
NCCdata <- dplyr::bind_rows(newDataList)

## Full dataset with missing observations (time points)
NCCdata <- dplyr::arrange(NCCdata, ID, MONTH)

## Tally by ID entries (interactive)
# group_by(NCCdata, ID) %>% summarize(N = n())

# Impute missing timepoints
fullNCC <- split(NCCdata, NCCdata[["ID"]]) %>%
    map(~ suppressMessages(right_join(.x,
                     tibble(ID = rep(unique(.x[["ID"]]), length(timeNumeric)),
                            MONTH = sort(timeNumeric)))))

NCCdf <- dplyr::bind_rows(fullNCC)

## Type convert characters to integers where possible
charVars <- vapply(NCCdf, is.character, logical(1L))
NCCdf[, charVars] <- readr::type_convert(NCCdf[, charVars])

## Serialize and store data
# readr::write_rds(NCCdf, path = "data/NCCmerged.rds")
