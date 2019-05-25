### Thioro and Zack are starting the file to analyse the cover crop data.
### This script requires set working directory to source file location.

library(tidyverse)
library(agricolae)

## Import and Manage Data

seeding_rate_table <- data.frame(CropTrt = c("SH", "SS", "VB",
                                             "SHSS", "SHVB", "SSVB"),
                                 seeding_rate = c(1, 1, 1, 0.5, 0.5, 0.5))

biomass_data <- read_csv("data/harvest-012019-TREC.csv") %>% 
  separate(PlotName, c("Location_CropTrt", "Row"), sep = "-") %>%
  separate(Location_CropTrt, c("Location", "CropTrt"), sep = "_", remove=FALSE) %>%
  mutate(LeavesStems_tha = LeavesStems_g_1m2*0.01) %>%
  left_join(seeding_rate_table)

total_biomass <- biomass_data %>%
  group_by(PlotID, CropSp, CropTrt, Location_CropTrt, seeding_rate) %>%
  summarize(avg_LeavesStems_tha = mean(LeavesStems_tha),
            avg_StemCount = mean(StemCounts))

total_biomass_summary <- total_biomass %>%
  group_by(Location_CropTrt, CropSp) %>%
  summarize(Avg_LeavesStems_tha = mean(avg_LeavesStems_tha),
            SD_LeavesStems_tha = sd(avg_LeavesStems_tha),
            Avg_LeavesStems_seed = mean(avg_LeavesStems_tha/seeding_rate),
            SD_LeavesStems_seed = sd(avg_LeavesStems_tha/seeding_rate),
            Avg_StemCount = mean(avg_StemCount),
            SD_StemCount = sd(avg_StemCount))

moisture_content <- read_csv("data/moisture-content.csv")
moisture_species <- moisture_content %>%
  group_by(CropSp) %>%
  summarize(avg_MC = mean(AGB_MC))


species <- c("SH", "SS", "VB")

## ANOVA for TotalBiomass
sink("output/Biomass_anova.txt")
for(s in species){
  moisture_conversion <- moisture_species %>%
    filter(CropSp == s)
  moisture_factor <- (100-moisture_conversion$avg_MC)/100
  
  data_for_analysis <- total_biomass %>%
    select(PlotID, CropSp, CropTrt, Location_CropTrt, 
           avg_LeavesStems_tha, avg_StemCount, seeding_rate) %>%
    mutate(avg_drymass = avg_LeavesStems_tha*moisture_factor) %>%
    filter(CropSp == s)
  
  test <- aov(avg_LeavesStems_tha ~  CropTrt + Location_CropTrt, data=data_for_analysis)
  print(paste(s, "- Fresh Mass"))
  print(summary(test))
  print(HSD.test(test, "CropTrt")$groups)
  print(HSD.test(test, "Location_CropTrt")$groups)
  
  test <- aov(avg_drymass ~  CropTrt + Location_CropTrt, data=data_for_analysis)
  print(paste(s, "- Dry Mass"))
  print(summary(test))
  print(HSD.test(test, "CropTrt")$groups)
  print(HSD.test(test, "Location_CropTrt")$groups)
  
  test <- aov(avg_LeavesStems_tha/seeding_rate ~  CropTrt + Location_CropTrt, data=data_for_analysis)
  print(paste(s, "- Fresh Mass / Seeding Rate"))
  print(summary(test))
  print(HSD.test(test, "CropTrt")$groups)
  print(HSD.test(test, "Location_CropTrt")$groups)
  
  test <- aov(avg_LeavesStems_tha/avg_StemCount ~  CropTrt + Location_CropTrt, data=data_for_analysis)
  print(paste(s, "- Fresh Mass / Stem Count"))
  print(summary(test))
  print(HSD.test(test, "CropTrt")$groups)
  print(HSD.test(test, "Location_CropTrt")$groups)
}
sink()

## Land Equivalent Ratio
data_for_LER <- total_biomass %>%
  select(Location, CropTrt, CropSp, avg_LeavesStems_g_1m2) %>%
  group_by(Location, CropTrt, CropSp) %>%
  summarize(site_avg_biomass = mean(avg_LeavesStems_g_1m2),
            site_sd_biomass = sd(avg_LeavesStems_g_1m2))

spp_combos <- data.frame(sp_1 = c("SH", "SH", "SS"), sp_2 = c("SS", "VB", "VB"))
LER_output <- data.frame(location = c(), sp_1 = c(), sp_2 = c(), LER = c())
for (l in unique(data_for_LER$Location)){
  location_data <- data_for_LER %>%
    filter(Location == l)
  for (i in 1:nrow(spp_combos)){
    mixed_1 <- location_data %>%
      filter(CropTrt == paste(spp_combos$sp_1[i], spp_combos$sp_2[i], sep="") &
               CropSp == spp_combos$sp_1[i])
    mono_1 <- location_data %>%
      filter(CropTrt == spp_combos$sp_1[i] & CropSp == spp_combos$sp_1[i])
    mixed_2 <- location_data %>%
      filter(CropTrt == paste(spp_combos$sp_1[i], spp_combos$sp_2[i], sep="") &
               CropSp == spp_combos$sp_2[i])
    mono_2 <- location_data %>%
      filter(CropTrt == spp_combos$sp_2[i] & CropSp == spp_combos$sp_2[i])
    LER <- round(mixed_1$site_avg_biomass / mono_1$site_avg_biomass + 
      mixed_2$site_avg_biomass / mono_2$site_avg_biomass, 3)
    LER_row <- data.frame(l, c(spp_combos[i, ]), LER)
    LER_output <- bind_rows(LER_output, LER_row)
  }
}

write_csv(LER_output, "output/LER.csv")