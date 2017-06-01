## Static datasets for data cleaning

## A patterns CharacterList for cleaning names in data
endings <- IRanges::CharacterList(
    Baseline = paste(c("YN$", "y$", "Ca$"), collapse = "|"),
    M12Scan = paste(c("YNY$", "Y$", "YY$", "CY$", "M12$"), collapse = "|"),
    M1Scan = paste(c("YN1$", "1$", "Y1$", "C1$", "M1$"), collapse = "|"),
    M6Scan = paste(c("YN6$", "6$", "Y6$", "C6$", "M6$"), collapse = "|"))

## Create codebook data frame for extracting data from variable names
infoFrame <- data.frame(
    pattern = c("^2[A-Z]*", "^3[A-Z]*", "^5[1-3]*", "^6[A-D]*", "^7[A-F]*",
        "^[23]A.", "^[23]B.", "^[23]C.", "^[23]D.",
        "^[23]E.", "^6A", "^6B", "^6C", "^6D",
        "^6[A-D]1", "^6[A-D]2", "^6[A-D]3",
        "^7A[1-4]", "^7B[1-4]", "^7C[1-4]", "^7D[1-4]", "^7E[1-4]", "^7F[1-4]",
        "^7[A-F]1", "^7[A-F]2", "^7[A-F]3", "^7[A-F]4",
        "^[23][A-E]1.", "^[23][A-E]2.", "^[23][A-E]3.",
        "^51[A-C]*", "^52[A-C]*", "^53[A-C]*",
        "^[23][A-E][1-3]A", "^[23][A-E][1-3]B1",
        "^[23][A-E][1-3]B2", "^[23][A-E][1-3]C",
        "^5[1-3]A", "^5[1-3]B1", "^5[1-3]B2", "^5[1-3]C"),
    interpretation = c("LeftHemi", "RightHemi", "PosteriorFossa",
                       "Intraventricular", "Cisternal",
                       "Frontal", "Temporal", "Parietal", "Occipital", "Basal",
                       "LateralVent", "IIIVent", "IVVent", "Doubtful",
                       "Active", "Transitional", "Inactive",
                       "RSylvian", "LSylvian", "Suprasellar", "PeriMesence",
                       "Prepontine", "Other",
                       "Active", "Transitional", "Inactive", "SubarachClust",
                       rep(c("Active", "Transitional", "Inactive"), 2),
                       rep(c("Parenchymal", "SubarachNumber",
                             "SubarachClust", "Indeterminate"), 2)),
    level = c(rep("First", 5), rep("Second", 9), rep("Third", 3),
              rep("Second", 6), rep("Third", 4), rep("Third", 3),
              rep("Second", 3), rep ("Fourth", 4),
              rep("Third", 4)),
    variableName = c(rep("Region", 5), rep("Area", 9), rep("Status", 3),
                     rep("Area", 6), rep("Status", 3), "Tissue",
                     rep("Status", 6), rep("Tissue", 8)),
    stringsAsFactors = FALSE
)
