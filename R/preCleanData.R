## Select variables for analysis and recode STATUS variable
fullData <- read.csv("data/NCClong.csv", header=TRUE, stringsAsFactors = FALSE)
ncc <- fullData %>% select(ID:STATUS, drug) %>% arrange(ID, LocCode, MONTH)

cleanLocationDF <- function(datframe) {
    stat4 <- datframe[["STATUS"]] == 4L
    if (any(stat4, na.rm = TRUE)) {
    datframe <- datframe[seq_len(which.max(stat4)), ]
    datframe <- datframe[!is.na(datframe[["STATUS"]]), ]
    }
    nas <- is.na(datframe[["STATUS"]])
    if (any(nas)) {
        restNA <- identical(which(nas), seq(which.max(nas), length(nas)))
        if (restNA) {
            datframe[["STATUS"]][which.max(nas)] <- 99L
        } else {
            if (sum(nas) > 1L) {
            datframe[["STATUS"]][6L - which.max(rev(nas))] <- 99L
            }
        }
    }
    datframe[!is.na(datframe[["STATUS"]]), ]
}

ncc <- dplyr::bind_rows(lapply(split(ncc, ncc[["IDLOC"]]), cleanLocationDF))

## Recode cyst location to dichotomous (Parenchymal/Not)
ncc$Parenchymal <- dplyr::recode(ncc$Tissue, "Parenchymal" = "Yes", .default = "No")

## Recode non-parenchymal regions as "No"
ncc$Parenchymal[is.na(ncc$Parenchymal)] <- "No"

# readr::write_csv(ncc, "data/NCCstatus.csv")

# DESCRIPTIVES ------------------------------------------------------------

## minimum number of recorded visits
group_by(ncc, ID) %>%
    summarize(nMonth = n_distinct(MONTH)) %>%
        select(nMonth) %>%
            min

## number of IDs with at least one dissolved cyst
group_by(ncc, ID) %>%
    summarize(dislv = sum(STATUS == 4), dlogic = (dislv != 0L)) %>%
        summarize(ndislv = sum(dlogic))

## Percent of patients that experience a dissolved cyst
round(101/117 * 100, 2)

## Number of single active cysts
filter(ncc, MONTH == 0L) %>% group_by(ID) %>%
    summarize(activ = sum(STATUS == 1L)) %>%
        summarize(n1cysts = sum(activ >= 1L))

## Percent of patients with single active cysts
round(62/117 * 100, 2)
