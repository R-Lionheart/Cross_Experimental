library(tidyverse)
source("src/Functions.R")
# Import all files
filenames <- RemoveCsv(list.files(path = "data_raw", pattern = "*.csv"))

for (i in filenames) {
  filepath <- file.path("data_raw", paste(i, ".csv", sep = ""))
  assign(make.names(i), read.csv(filepath, stringsAsFactors = FALSE))
}

matching.variable = "HILIC"
#####

columns.to.drop <- c('Average.Rt.min.', 'Formula', 'Ontology', 'INCHIKEY', 
                     'SMILES', 'Isotope.tracking.parent.ID', 'Isotope.tracking.weight.number',
                     'MS1.isotopic.spectrum', 'MS.MS.spectrum', 'Average.Mz', 'Post.curation.result', 
                     'Fill..', 'Annotation.tag..VS1.0.', 'RT.matched',
                     'm.z.matched', 'MS.MS.matched', 'Manually.modified', 'Total.score', 
                     'RT.similarity', 'Dot.product', 'Reverse.dot.product', 'Fragment.presence..')

# Set header, filter unknowns ---------------------------------------
runs <- grep(matching.variable, names(.GlobalEnv), value = TRUE, ignore.case = TRUE)
runlist <- do.call("list", mget(runs))

headers.set <- lapply(names(runlist), function(x) SetHeader(runlist[[x]]))
names(headers.set) <- runs

for (df in seq_along(headers.set)) { 
  headers.set[[df]] <- headers.set[[df]] %>% filter(!Metabolite.name == "Unknown")
  headers.set[[df]] <- headers.set[[df]] %>% select(-one_of(columns.to.drop))
}

# Change variable classes -------------------------------------------------
classes.changed <- lapply(names(headers.set), function(x) ChangeClasses(headers.set[[x]]))
names(classes.changed) <- runs

list2env(classes.changed, globalenv())


# Rearrange datasets ------------------------------------------------------

# Positive B12
B12.Area.pos <- RearrangeDatasets(Area_HILIC.POS_B12.Inc, parameter = "Area.Value")
B12.Mz.pos   <- RearrangeDatasets(Mz_HILIC.POS_B12.Inc, parameter = "Mz.Value")
B12.RT.pos   <- RearrangeDatasets(RT_HILIC.POS_B12.Inc, parameter = "RT.Value")
B12.SN.pos   <- RearrangeDatasets(SN_HILIC.POS_B12.Inc, parameter = "SN.Value")

# Positive Eddy Transect
Transect.Area.pos <- RearrangeDatasets(Area_HILICPos_EddyTransect, parameter = "Area.Value")
Transect.Mz.pos   <- RearrangeDatasets(Mz_HILICPos_EddyTransect, parameter = "Mz.Value")
Transect.RT.pos   <- RearrangeDatasets(RT_HILICPos_EddyTransect, parameter = "RT.Value")
Transect.SN.pos   <- RearrangeDatasets(SN_HILICPos_EddyTransect, parameter = "SN.Value")

# Negative B12
B12.Area.neg <- RearrangeDatasets(Area_HILIC.NEG_B12.Inc, parameter = "Area.Value")
B12.Mz.neg   <- RearrangeDatasets(Mz_HILIC.NEG_B12.Inc, parameter = "Mz.Value")
B12.RT.neg   <- RearrangeDatasets(RT_HILIC.NEG_B12.Inc, parameter = "RT.Value")
B12.SN.neg   <- RearrangeDatasets(SN_HILIC.NEG_B12.Inc, parameter = "SN.Value")

# Negative Eddy Transect
Transect.Area.neg <- RearrangeDatasets(Area_HILICNeg_EddyTransect, parameter = "Area.Value")
Transect.Mz.neg   <- RearrangeDatasets(Mz_HILICNeg_EddyTransect, parameter = "Mz.Value")
Transect.RT.neg   <- RearrangeDatasets(RT_HILICNeg_EddyTransect, parameter = "RT.Value")
Transect.SN.neg   <- RearrangeDatasets(SN_HILICNeg_EddyTransect, parameter = "SN.Value")

# Combine to one dataset --------------------------------------------------
combined.B12.pos <- B12.Area.pos %>%
  left_join(B12.Mz.pos) %>%
  left_join(B12.SN.pos) %>%
  left_join(B12.RT.pos) %>%
  mutate(Column = "HILICPos") %>%
  mutate(Dataset = "B12_Incubation") %>%
  select(Replicate.Name, Column, Dataset, Area.Value, Mz.Value, RT.Value, SN.Value, everything())

combined.transect.pos <- Transect.Area.pos %>%
  left_join(Transect.Mz.pos) %>%
  left_join(Transect.SN.pos) %>%
  left_join(Transect.RT.pos) %>%
  mutate(Column = "HILICPos") %>%
  mutate(Dataset = "Eddy_Transect") %>%
  select(Replicate.Name, Column, Dataset, Area.Value, Mz.Value, RT.Value, SN.Value, everything())


combined.B12.neg <- B12.Area.neg %>%
  left_join(B12.Mz.neg) %>%
  left_join(B12.SN.neg) %>%
  left_join(B12.RT.neg) %>%
  mutate(Column = "HILICNeg") %>%
  mutate(Dataset = "B12_Incubation") %>%
  select(Replicate.Name, Column, Dataset, Area.Value, Mz.Value, RT.Value, SN.Value, everything())

combined.transect.neg <- Transect.Area.neg %>%
  left_join(Transect.Mz.neg) %>%
  left_join(Transect.SN.neg) %>%
  left_join(Transect.RT.neg) %>%
  mutate(Column = "HILICNeg") %>%
  mutate(Dataset = "Eddy_Transect") %>%
  select(Replicate.Name, Column, Dataset, Area.Value, Mz.Value, RT.Value, SN.Value, everything())
  
  
combined.B12 <- rbind(combined.B12.pos, combined.B12.neg) %>%
  mutate(Metabolite.name = ifelse(str_detect(Metabolite.name, "Ingalls_"), sapply(strsplit(Metabolite.name, "_"), `[`, 2), Metabolite.name)) 

combined.transect <- rbind(combined.transect.pos, combined.transect.neg) %>%
  mutate(Metabolite.name = ifelse(str_detect(Metabolite.name, "Ingalls_"), sapply(strsplit(Metabolite.name, "_"), `[`, 2), Metabolite.name)) 

combined <- rbind(combined.transect, combined.B12) %>%
  filter(!str_detect(Metabolite.name, ","))
combined$Replicate.Name <- gsub("^.{0,1}", "", combined$Replicate.Name)

# Fix Time 0 Labeling issues -----------------------------------------------
combined.final <- combined %>%
  mutate(Replicate.Name = recode(Replicate.Name, 
                                 "171002_Smp_IT0_1" ="171002_Smp_IL1IT0_1", 
                                 "171002_Smp_IT0_2" = "171002_Smp_IL1IT0_2",
                                 "171002_Smp_IT0_3" = "171002_Smp_IL1IT0_3",
                                 "171009_Smp_IT05um_1" = "171009_Smp_IL1IT05um_1",
                                 "171009_Smp_IT05um_2" = "171009_Smp_IL1IT05um_2",
                                 "171009_Smp_IT05um_3" = "171009_Smp_IL1IT05um_3",
                                 "171016_Smp_IT0_1" = "171016_Smp_IL2IT0_1",
                                 "171016_Smp_IT0_2" = "171016_Smp_IL2IT0_2",
                                 "171016_Smp_IT0_3" = "171016_Smp_IL2IT0_3",
                                 "171023_Smp_IT05um_1" = "171023_Smp_IL2IT05um_1",
                                 "171023_Smp_IT05um_2" = "171023_Smp_IL2IT05um_2",
                                 "171023_Smp_IT05um_3" = "171023_Smp_IL2IT05um_3")) 
currentDate <- Sys.Date()
csvFileName <- paste("data_processed/MSDial_B12_Transect_combined_", currentDate, ".csv", sep = "")

write.csv(combined.final, csvFileName, row.names = FALSE)
rm(list = ls())



