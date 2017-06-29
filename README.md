# NCC

Neurocysticercosis Project

# Instructions for Reproducibility

1. Please run the `R/loadPackages.R` script interactively. Please make sure
that you have all packages installed on your system.

2. The `R/helperData.R` script creates data objects for cleaning and parsing
the column names in the datasets. This is not to be modified unless you are
familiar with the datasets.

3. To download the resources from `Dropbox`, please run the `R/download.R`
script interactively. This requieres a GitHub package download `karthik/rdrop2`
to connect to the `Dropbox` API. Please note that folder access is required
in order to be able to download the data.

4. If you do not have the serialized version of the dataset, please run the
`R/readCleanMerge.R` script with the `SAV` files downloaded from `Dropbox`.
The SPSS (`SAV`) files can be obtained from `Dropbox` with folder permissions.
This script will then save the datasets in a serialized `rds`, `sav`, or `csv`
format.

5. Run the `R/validation.R` script for obtaining some helper functions to
validate and recode the data.

6. Open `R/analyzeNCC.R` and run the analysis interactively. Please note, this
script will source the previous scripts.

# Issues

Please submit any issues at our GitHub repository issue tracker:
https://github.com/LiNk-NY/NCC/issues
