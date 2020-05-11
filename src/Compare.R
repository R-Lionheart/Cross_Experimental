library(tidyverse)
source("src/Functions.R")


# Enter user data for file.pattern
file.pattern <- "combined"

# Import required files
filename <- RemoveCsv(list.files(path = "data_processed", pattern = regex(file.pattern, ignore_case = TRUE)))

filepath <- file.path("data_processed", paste(filename, ".csv", sep = ""))
assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE))

full.file <- MSDial_B12_Transect_combined_2020.05.04 %>%
  select(Metabolite.name, Replicate.Name, Column, Dataset, Area.Value:SN.Value)