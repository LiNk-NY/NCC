## Load packages
source("R/loadPackages.R")

## Load helper dataset
source("R/helperData.R")

## Load validation functions
source("R/validation.R")

## Load clean data
NCCdf <- readr::read_rds("data/NCCmerged.rds")

colNames <- names(NCCdf)

hitVec <- apply(vapply(infoFrame[["pattern"]], function(x) {
        grepl(x, colNames)
    }, logical(length(colNames))), 1, any)
names(hitVec) <- colNames

## Check pattern is picking up correct targets
(pickedUp <- Filter(isTRUE, hitVec))

## Check left out
setdiff(names(NCCdf), names(pickedUp))

## Obtain hit matrix for matching patterns in variable names
namesHitMat <- vapply(infoFrame[["pattern"]], function(x) grepl(x, colNames),
                      logical(length(colNames)))
rownames(namesHitMat) <- colNames
namesHitMat <- namesHitMat[rowSums(namesHitMat) > 0L, ]

## Matches in infoFrame for each variable (see interpretation column)
decodedVariable <- apply(namesHitMat, 1, function(g) infoFrame[g, ])

decodedVariable <- lapply(decodedVariable,
   function(z) {
       suppressMessages(
           right_join(z, data.frame(variableName =
                        c("Region", "Area", "Status", "Tissue"),
                        stringsAsFactors = FALSE)
                   ))
       })

## Take interpretation column and create new data from variables
locationDat <- t(vapply(decodedVariable, function(x) x[["interpretation"]], character(4L)))
locationDat <- as.data.frame(locationDat, stringsAsFactors = FALSE)
colnames(locationDat) <- unique(infoFrame[["variableName"]])
locationDat$Code <- rownames(locationDat)
locationDat$Status <- factor(locationDat$Status,
                      levels = c("Active", "Transitional", "Inactive"),
                      ordered = TRUE)
locationDat$CombCode <-  gsub("(^[235]*.)([1-3])([A-Z]*.)", "\\1X\\3",
                            locationDat$Code) %>%
                                gsub("(^[67][A-F])([1-3])", "\\1X", .)
locationDat$Area[locationDat$Region == "PosteriorFossa"] <- "BS/Cerebellum"
locationDat <- locationDat[!is.na(locationDat$Status) &
                               grepl("X", locationDat$CombCode), ]

regions <- arrange(locationDat, Tissue, Status) %>% select(-Tissue) %>%
    split(., list(.$Region, .$Area, .$CombCode)) %>%
        Filter(function(x) nrow(x) != 0L, .) %>% map(function(x) x$Code)

colsInterest <- lapply(regions, function(x) c("ID", "MONTH", x))

## Long and skinny format
## FIX HERE
dataByLOC <- lapply(colsInterest, function(region) {
unite(NCCdf[, region], IDMONTH, c(ID, MONTH)) %>%
    gather(LOCATION, COUNT, -IDMONTH) %>%
    separate(IDMONTH, c("ID", "MONTH")) %>%
    arrange(ID, MONTH, LOCATION)
})

dataByLOC <- dplyr::bind_rows(dataByLOC)
NCCFULL <- left_join(dataByLOC, locationDat, by = c("LOCATION" = "Code"))

NCCwide <- spread(NCCFULL, key = Status, value = COUNT)

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

validLocs <- vapply(strsplit(names(IDbyRegion), split = "\\."),
                    function(vec) vec[[3]], character(1L))
names(validLocs) <- names(IDbyRegion)

dataByCode <- dplyr::filter(NCCwide, NCCwide$CombCode %in% validLocs) %>%
    split(., .$CombCode)

for (i in seq_along(dataByCode)) {
    dataByCode[[i]] <- dataByCode[[i]][dataByCode[[i]]$ID %in% IDbyRegion[[i]], ]
}

dataByCode <- dplyr::bind_rows(dataByCode)
