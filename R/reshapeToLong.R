## Load packages
source("R/loadPackages.R")

## Load helper dataset
source("R/helperData.R")

## Load validation functions
source("R/validation.R")

## Load clean data
NCCdf <- readr::read_rds("data/NCCwide.rds")

colNames <- names(NCCdf)

hitVec <- apply(vapply(infoFrame[["pattern"]], function(x) {
        grepl(x, colNames)
    }, logical(length(colNames))), 1, any)
names(hitVec) <- colNames

## Other variables
restVars <- NCCdf[, colNames[!hitVec]]

## Check pattern is picking up correct targets
(pickedUp <- Filter(isTRUE, hitVec))

## Check left out
identical(setdiff(names(NCCdf), names(pickedUp)), colNames[!hitVec])
colNames[!hitVec]

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
                      levels = c("ACTIVE", "TRANSITIONAL", "INACTIVE"),
                      ordered = TRUE)
locationDat$LocCode <-  gsub("(^Q[235]*.)([1-3])([A-Z]*.)", "\\1X\\3",
                            locationDat$Code) %>%
                                gsub("(^Q[67][A-F])([1-3])", "\\1X", .)
locationDat$Area[locationDat$Region == "PosteriorFossa"] <- "BS/Cerebellum"
locationDat <- locationDat[!is.na(locationDat$Status) &
                               grepl("X", locationDat$LocCode), ]

regions <- arrange(locationDat, Tissue, Status) %>% select(-Tissue) %>%
    split(., list(.$Region, .$Area, .$LocCode)) %>%
        Filter(function(x) nrow(x) != 0L, .) %>% map(function(x) x$Code)

colsInterest <- lapply(regions, function(x) c("ID", "MONTH", x))

## Long and skinny format

dataByLOC <- lapply(colsInterest, function(region) {
unite(NCCdf[, region], IDMONTH, c(ID, MONTH)) %>%
    gather(LOCATION, COUNT, -IDMONTH) %>%
    separate(IDMONTH, c("ID", "MONTH")) %>%
    arrange(ID, MONTH, LOCATION)
})

dataByLOC <- dplyr::bind_rows(dataByLOC)

NCCFULL <- left_join(dataByLOC, locationDat, by = c("LOCATION" = "Code"))

NCCwide <- NCCFULL %>% unite(IDLOC, c(ID, LocCode)) %>% select(-LOCATION) %>%
    spread(key = Status, value = COUNT) %>%
    mutate(MONTH = type.convert(MONTH)) %>%
    arrange(IDLOC, MONTH) %>% separate(IDLOC, c("ID", "LocCode"))

## Example data chunk by region
regionID <- lapply(regions, function(reg) split(NCCdf[, reg], NCCdf$ID))

## Function from validation file
validCystMatrix(regionID[[1L]][[1]])

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

dataByCode <- dplyr::filter(NCCwide, NCCwide$LocCode %in% validLocs) %>%
    split(., .$LocCode)

for (i in seq_along(dataByCode)) {
    dataByCode[[i]] <- dataByCode[[i]][dataByCode[[i]]$ID %in%
                                           IDbyRegion[[i]], ]
}

dataByCode <- dplyr::bind_rows(dataByCode)

dataByCode <- dataByCode %>% mutate(IDVAR = paste0(ID, "_", MONTH))

restVars <- restVars %>% mutate(IDVAR = paste0(ID, "_", MONTH))

FullNCC <- left_join(dataByCode, restVars %>% select(-c(MONTH, ID)),
                     by = "IDVAR") %>% select(-IDVAR)

FullNCC <- FullNCC %>%  mutate(ID = type.convert(ID))  %>%
    add_column(STATUS = apply(
        FullNCC[, c("ACTIVE", "TRANSITIONAL", "INACTIVE")], 1L,
        function(x) {.stageRecode(x)}), .after = "Tissue")

drug <- readr::read_csv("data/drugVars.csv")

NCClong <- left_join(FullNCC, drug, by = "ID")

NCClong <- NCClong %>% mutate(IDLOC = paste0(ID, "_", LocCode)) %>%
    select(ID, LocCode, IDLOC, everything()) %>%
    mutate_at(vars(drug), function(x) case_when(x == 2 ~ 0, TRUE ~ 1))

## Replace wide dataset
# readr::write_csv(NCClong, "data/NCClong.csv")
# haven::write_sav(NCClong, "data/NCClong.sav")
