## Load packages
source("R/loadPackages.R")

## Load helper dataset
source("R/helperData.R")

## Load validation functions
source("R/validation.R")

## Load clean data
NCCdf <- readr::read_rds("data/NCCmerged.rds")

## Find variable names for right and left hemispheres
findHemi <- grep("^[23]", names(NCCdf), value = TRUE)

## Shorten location names
locNames <- gsub("([23][A-Z])(.)([A-Z]*)", "\\1.\\3", findHemi)
names(locNames) <- findHemi

## Check there are no matches in 3rd info spot in coded locations (STATUS)
any(vapply(infoFrame[["pattern"]], function(x) { grepl(x, locNames) },
    logical(length(locNames)))[, infoFrame[infoFrame$variableName=="Status",
                                           "pattern"]])

## Obtain hit matrix for matching patterns in variable names
namesHitMat <- vapply(infoFrame[["pattern"]], function(x) grepl(x, findHemi),
                      logical(length(findHemi)))
rownames(namesHitMat) <- findHemi

## Matches in infoFrame for each variable (see interpretation column)
decodedVariable <- apply(namesHitMat, 1, function(g) infoFrame[g, ])

## Take interpretation column and create new data from variables
locationDat <- t(vapply(decodedVariable, function(x) x[["interpretation"]], character(4L)))
locationDat <- as.data.frame(locationDat, stringsAsFactors = FALSE)
colnames(locationDat) <- unique(infoFrame[["variableName"]])
locationDat$Code <- rownames(locationDat)
locationDat$Status <- factor(locationDat$Status,
                      levels = c("Active", "Transitional", "Inactive"),
                      ordered = TRUE)

regions <- arrange(locationDat, Tissue, Status) %>%
    split(., list(.$Hemisphere, .$Area, .$Tissue)) %>%
    map(function(x) x$Code)

colsInterest <- vapply(regions, function(x) c("ID", "MONTH", x), character(5L))
test <- colsInterest[, 1]

## Long and skinny format
unite(NCCdf[, test], IDMONTH, c(ID, MONTH)) %>%
    gather(LOCATION, COUNT, -IDMONTH) %>%
    separate(IDMONTH, c("ID", "MONTH")) %>%
    mutate(STAGE = factor(gsub("([23][A-Z])(.)([A-Z]*)", "\\2", LOCATION),
                          levels = 1:3,
                          labels = c("Active", "Transitional", "Inactive"))) %>%
    arrange(ID, MONTH, STAGE)

## Example data chunk by region
regionID <- lapply(regions, function(reg) split(NCCdf[, reg], NCCdf$ID))

## Function from validation file
validCystMatrix(regionID[[1L]][[1]])

## Check how many code chunks per region and ID are valid
lapply(regionID, function(reg) { sum(vapply(reg, function(x) {validCystMatrix(x)}, logical(1L)))})

## Breakdown valid chunks by region and ID
lapply(regionID, function(reg) { vapply(reg, function(x) {validCystMatrix(x)}, logical(1L))})

## Get IDs that have a valid matrix
validIDs <- lapply(regionID, function(reg) {
    names(Filter(isTRUE, vapply(reg,
        function(x) {
            validCystMatrix(x)
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
