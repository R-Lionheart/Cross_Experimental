# Comparing Time 0 on B12 Incubations with Eddy Transect to determine similarity

library(tidyverse)
library(vegan)

all.datasets <- read.csv("data_processed/MSDial_B12_Transect_combined_2020-06-16.csv", stringsAsFactors = FALSE) %>%
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
  filter(Dataset == "Eddy_Transect") %>%
  filter(str_detect(Replicate.Name, "MS6|MS12")) %>%
  filter(str_detect(Replicate.Name, "MS6C315m|MS12C115m"))


Eddy_Time0 <- Time0 %>%
  rbind(Eddy) %>%
  separate(Replicate.Name, into = c("one", "two", "SampID", "four")) %>%
  group_by(Metabolite.name, SampID) %>%
  mutate(Average.Area = mean(Area.Value, na.rm = TRUE)) %>%
  mutate(Cyclone = ifelse(str_detect(SampID, "MS12|IL1"), 
                          "Cyclonic", "Anticyclonic")) %>%
  mutate(Size.Fraction = ifelse(Dataset == "B12_Incubation" & str_detect(SampID, "5um"), "5um_SizeFraction", 
                                ifelse(Dataset == "B12_Incubation" & !str_detect(SampID, "5um"), "0.2_SizeFraction", "Transect"))) %>%
  group_by(Metabolite.name) %>%
  filter(!length(unique(Dataset)) < 2) %>%
  filter(Metabolite.name %in% compounds) %>%
  select(Metabolite.name, SampID, Dataset, Cyclone, Size.Fraction, Average.Area) %>%
  unique()

ggplot(data = Eddy_Time0, aes(x = Size.Fraction, y = Average.Area, fill = Metabolite.name)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Cyclone) +
  ggtitle("B12 Incubation T0 and Eddy Transect Comparison")

# With error bars
Eddy_Time0_Anticyclonic <- Eddy_Time0 %>%
  filter(Cyclone = "Anticyclonic")
Eddy_Time0_Cyclonic <- Eddy_Time0 %>%
  filter(Cyclone = "Cyclonic")

ggplot(Eddy_Time0_Anticyclonic, aes(Size.Fraction, Average.Area, fill= Metabolite.name)) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge") +
  ggtitle("B12 Incubation T0 and Eddy Transect Comparison")



## Normalization
Eddy_Time0_normalized <- Eddy_Time0 %>%
  select(-Dataset, -Cyclone) %>%
  pivot_wider(names_from = Metabolite.name, values_from = Average.Area)
EddytoJoin <- Eddy_Time0_normalized 
Eddy_Time0_normalized <- Eddy_Time0_normalized %>% select(-1) %>% as.matrix()
Eddy_Time0_normalized <- decostand(Eddy_Time0_normalized, method = "standardize", 1) %>%
  as.data.frame()
rownames(Eddy_Time0_normalized) <- EddytoJoin$SampID

Eddy_Time0_normalized <- Eddy_Time0_normalized %>% 
  rownames_to_column("SampID") %>%
  pivot_longer(-SampID, names_to = "Metabolite.name", values_to = "Average.Area.Norm") %>%
  mutate(Cyclone = ifelse(str_detect(SampID, "MS12|IL1"), 
                          "Cyclonic", "Anticyclonic")) %>%
  mutate(Dataset = ifelse(str_detect(SampID,"IL1|IL2"), "B12_Incubation", "Eddy_Transect")) 

ggplot() + 
  geom_tile(data = Eddy_Time0_normalized, 
            aes(x = SampID, y = Metabolite.name, fill = Average.Area.Norm)) +
  theme(#axis.text.y  = element_blank(),
    axis.text.x = element_text(angle = 270, vjust = 0, hjust = 0),
    strip.text = element_blank()) +
  ggtitle("Normalized B12 Incubation T0 and Eddy Transect Comparison")
  