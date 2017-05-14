## Static datasets for data cleaning

## A patterns CharacterList for cleaning names in data
endings <- IRanges::CharacterList(
    Baseline = paste(c("YN$", "y$", "Ca$"), collapse = "|"),
    M12Scan = paste(c("YNY$", "Y$", "YY$", "CY$", "M12$"), collapse = "|"),
    M1Scan = paste(c("YN1$", "1$", "Y1$", "C1$", "M1$"), collapse = "|"),
    M6Scan = paste(c("YN6$", "6$", "Y6$", "C6$", "M6$"), collapse = "|"))

## Create codebook data frame for extracting data from variable names
infoFrame <- data.frame(
    pattern = c("^2[A-Z]*", "^3[A-Z]*", "^[23]A.", "^[23]B.",
                "^[23]C.", "^[23]D.", "^[23]E.", "^[23][A-E]1.",
                "^[23][A-E]2.", "^[23][A-E]3.", "^[23][A-E][1-3]A",
                "^[23][A-E][1-3]B1", "^[23][A-E][1-3]B2", "^[23][A-E][1-3]C"),
    interpretation = c("LeftHemi", "RightHemi", "Frontal", "Temporal",
                       "Parietal", "Occipital", "Basal", "Active",
                       "Transitional", "Inactive", "Parenchymal",
                       "SubarachNumber", "SubarachClust", "Indeterminate"),
    level = c(rep("First", 2), rep("Second", 5), rep("Third", 3),
              rep ("Fourth", 4)),
    variableName = c(rep("Hemisphere", 2), rep("Area", 5),
                     rep("Status", 3), rep("Tissue", 4)),
    stringsAsFactors = FALSE
)
