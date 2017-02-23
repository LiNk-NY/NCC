## Download data from Dropbox

## Install dropbox API package
## devtools::install_github("karthik/rdrop2")

library(rdrop2)
library(dplyr)

## Authenticate to Dropbox API at FIRST RUN
# drop_auth()
# drop_acc() %>% select(uid, display_name, email_verified, quota_info.quota)
dataFiles <- rdrop2::drop_dir("NCC_PROJECT/DATA")[["path"]]
savData <- grepl("^[^A].+\\.sav$", basename(dataFiles))
dataFiles <- dataFiles[savData]

invisible(lapply(dataFiles, function(file) {
    drop_get(file, local_file = file.path("data", basename(file)))
})
)

## Rename file
file.rename("data/Baseline SCan & Drug.sav", "data/BaselineScanDrug.sav")
