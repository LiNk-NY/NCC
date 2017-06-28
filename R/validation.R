## Data validation functions

## Checker functions
## No sum of rows greater than 1
.sumRowsGT1 <- function(x) {
    !any(apply(x, 1, function(g) { sum(g) > 1 }), na.rm = TRUE)
}

## Total matrix sum is not zero
.totalSumNotZero <- function(x) {
    sum(x, na.rm = TRUE) != 0L
}

.notAllNA <- function(x) {
    !all(is.na(x))
}

.valueCheck <- function(x) {
    naRows <- rowSums(is.na(x)) == ncol(x)
    if (sum(naRows) >= (nrow(x) - 1L))
        return(FALSE)
    x <- x[!naRows, ]
    ##check that row sums are conserved
    tryCatch({xdiff <- cbind(0, (x[-1, ] - x[-nrow(x), ]))}, error = function(e) browser())
    checklog <- vector("logical")
    for (i in 2:ncol(xdiff)) {
        ##the change in one column plus the sum of changes of columns to the left
        ##must be non-increasing - ie, if it increases, then the decreases to its
        ##left must be equal or greater
        checklog <- c(checklog, all(xdiff[, i] + rowSums(xdiff[, 1:(i-1), drop=FALSE]) <= 0))
    }
    all(checklog)
}

## Create checker function for data chunks (region by ID)
validCystMatrix <- function(x) {
    stopifnot(nrow(x) == 4L)
    if (.notAllNA(x)) {
        res <- all(.sumRowsGT1(x), .totalSumNotZero(x), .valueCheck(x))
    } else {
        res <- FALSE
    }
    return(res)
}

.stageRecode <- function(x) {
    if (sum(x, na.rm = TRUE) == 1L)
        res <- which(x == 1L)
    else if (any(is.na(x)))
        res <- NA
    else if (all(x == 0L))
        res <- 4L
    res
}
