library(RCurl)
library(tidyverse)
source("src/Functions.R")

# Import processed B12 Incubation files
filenames <- RemoveCsv(list.files(path = "data_processed", pattern = "combined"))

for (i in filenames) {
  filepath <- file.path("data_processed", paste(i, ".csv", sep = ""))
  assign(make.names(i), read.csv(filepath, stringsAsFactors = FALSE))
}

B12_Only <- MSDial_B12_Transect_combined_2020.06.16 %>%
  mutate(Metabolite.name = recode(Metabolite.name, 
                                  "2-Hydroxy-4-(methylthio)butyric" = "2-Hydroxy-4-(methylthio)butyric_acid")) %>%
  filter(Dataset == "B12_Incubation") %>%
  select(Metabolite.name, Replicate.Name, Column, Area.Value) 

# Import most recent standards sheet from Github
Raw.Standards <- getURL("https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards_NEW.csv")
Ingalls.Standards <- read.csv(text = Raw.Standards, stringsAsFactors = FALSE) %>%
  select(Compound.Name, Compound.Name_old) %>%
  rename(Metabolite.name = Compound.Name_old) %>%
  unique()

Name.Replacement <- B12_Only %>%
  left_join(Ingalls.Standards)

Treatment <- B12_Only %>%
  ungroup() %>%
  select(Replicate.Name) %>%
  separate(Replicate.Name, into = c("Date", "runtype", "Supergroup", "replicate"), remove = FALSE) %>%
  mutate(Control.Status = ifelse(str_detect(Supergroup, "IT0"),
                                 "Incubation", ifelse(str_detect(Supergroup, "DSW"), "DeepSeaWater", 
                                                      ifelse(str_detect(Supergroup, "Control"), "Control", "Treatments")))) %>%
  mutate(Treatment.Status = ifelse(Control.Status == "Control", "Control",
                                   ifelse(Control.Status == "DeepSeaWater", "DeepSeaWater",
                                          ifelse(Control.Status == "Incubation", "TimeZero",
                                                 ifelse(str_detect(Supergroup, "DMBnoBT"), "DMBnoB12",
                                                        ifelse(str_detect(Supergroup, "WBT"), "B12",
                                                               ifelse(str_detect(Supergroup, "DMB"), "DMB", "noB12"))))))) %>%
  select(Replicate.Name, Control.Status, Treatment.Status, Supergroup)

AllB12 <- B12_Only %>%
  left_join(Treatment) %>%
  unique()