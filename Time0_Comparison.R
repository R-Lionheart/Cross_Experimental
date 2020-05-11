# Comparing Time 0 on B12 Incubations with Eddy Transect to determine similarity

library(tidyverse)

all.datasets <- read.csv("data_processed/MSDial_B12_Transect_combined_2020-05-04.csv", stringsAsFactors = FALSE) %>%
  select(Metabolite.name, Replicate.Name, Column, Dataset, Area.Value) %>%
  filter(!str_detect(Replicate.Name,
                     "ProcessBlk|Sept29QC|TruePooWeek1|TruePooWeek2|TruePooWeek3|TruePooWeek4|DSW700m|Std|Blk")) %>%
  filter(!str_detect(Replicate.Name, regex("poo", ignore_case = TRUE))) %>%
  filter(!str_detect(Metabolite.name, ","))

compounds <- c("Glutamic acid", "Glutamine", "Ketoglutaric Acid", "N(e)-Acetyl-Lysine")

Time0 <- all.datasets %>%
  filter(Dataset == "B12_Incubation") %>%
  filter(str_detect(Replicate.Name, "IT0"))

Eddy <- all.datasets %>%
  filter(Dataset == "Eddy_Transect")


Eddy_Time0 <- Time0 %>%
  rbind(Eddy) %>%
  group_by(Metabolite.name) %>%
  filter(!length(unique(Dataset)) < 2) %>%
  filter(Metabolite.name %in% compounds)
  

ggplot(data = Eddy_Time0, aes(x = Metabolite.name, y = Area.Value, fill = Dataset)) +
  geom_bar(stat = "identity")