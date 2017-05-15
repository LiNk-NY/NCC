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

.valueCheck <- function(x) {
    ##check that row sums are conserved
    xdiff <- cbind(0, (x[-1, ] - x[-nrow(x), ]))
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
    all(.sumRowsGT1(x), .totalSumNotZero(x), .valueCheck(x))
}

