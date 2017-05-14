## Load packages
source("R/loadPackages.R")

## Load helper dataset
source("R/helperData.R")

## Load clean data
NCCdf <- readr::read_rds("data/NCCmerged.rds")

## Find variable names for right and left hemispheres
findHemi <- grep("^[23]", names(NCCdf), value = TRUE)

## Obtain hit matrix for matching patterns in variable names
namesHitMat <- vapply(infoFrame[["pattern"]], function(x) grepl(x, findHemi),
                      logical(length(findHemi)))
rownames(namesHitMat) <- findHemi

## Matches in infoFrame for each variable (see interpretation column)
decodedVariable <- apply(namesHitMat, 1, function(g) infoFrame[g, ])

## Take interpretation column and create new data from variables
vars <- t(vapply(decodedVariable, function(x) x[["interpretation"]], character(4L)))
vars <- as.data.frame(vars, stringsAsFactors = FALSE)
colnames(vars) <- unique(infoFrame[["variableName"]])
vars$Code <- rownames(vars)
vars$Status <- factor(vars$Status,
                      levels = c("Active", "Transitional", "Inactive"),
                      ordered = TRUE)

regions <- arrange(vars, Tissue, Status) %>%
    split(., list(.$Hemisphere, .$Area, .$Tissue)) %>%
    map(function(x) x$Code)

## Example data chunk by region
regionID <- lapply(regions, function(reg) split(NCCdf[, reg], NCCdf$ID))

## Checker functions
## No sum of rows greater than 1
.sumRowsGT1 <- function(x) {
    !any(apply(x, 1, function(g) { sum(g) > 1 }), na.rm = TRUE)
}

## Total matrix sum is not zero
.totalSumZero <- function(x) {
    sum(x, na.rm = TRUE) != 0L
}

## Create checker function for data chunks (region by ID)
.validCystMatrix <- function(x) {
    stopifnot(nrow(x) == 4L)
    all(.sumRowsGT1(x), .totalSumZero(x))
}

.validCystMatrix(regionID[[1L]][[1]])

## Check how many code chunks per region and ID are valid
lapply(regionID, function(reg) { sum(vapply(reg, function(x) {.validCystMatrix(x)}, logical(1L)))})

## Breakdown valid chunks by region and ID
lapply(regionID, function(reg) { vapply(reg, function(x) {.validCystMatrix(x)}, logical(1L))})

## Get IDs that have a valid matrix
validIDs <- lapply(regionID, function(reg) {
    names(Filter(isTRUE, vapply(reg,
        function(x) {
            .validCystMatrix(x)
            }, logical(1L))))
    })

IDbyRegion <- Filter(length, validIDs)

regionReduce <- regionID[names(IDbyRegion)]

stopifnot(identical(names(IDbyRegion), names(regionReduce)))

validMats <- lapply(seq_along(regionReduce), function(matList, i, vecList) {
    matList[[i]][vecList[[i]]]
}, matList = regionReduce, vecList = IDbyRegion)

names(validMats) <- names(IDbyRegion)

## Checking example region and ID
regionID$RightHemi.Temporal.SubarachNumber$`1109`
