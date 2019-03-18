### Thioro and Zack are starting the file to analyse the cover crop data.
### This script requires set working directory to source file location.

library(tidyverse)
library(agricolae)

## Import and Manage Data
biomass_data <- read_csv("data/harvest-092018-TREC.csv")

biomass_data <- biomass_data %>%
  mutate(LeavesStems_tha = LeavesStems_kg_1m2*10, 
         Biomass_tha = round(Biomass_kg_9m2*10/9, 2))

total_biomass <- biomass_data %>%
  group_by(PlotID, Location, RowNum, ColNum, CropTrt) %>%
  summarize(LeavesStems_kg_1m2 = sum(LeavesStems_kg_1m2),
            Biomass_kg_9m2 = sum(Biomass_kg_9m2))

moisture_content <- read_csv("data/moisture-content.csv")
moisture_species <- moisture_content %>%
  group_by(CropSp) %>%
  summarize(avg_MC = mean(AGB_MC))

species <- c("SH", "SS", "VB")

## ANOVA for TotalBiomass
sink("output/Biomass_9m2_anova.txt")
for(s in species){
  moisture_conversion <- moisture_species %>%
    filter(CropSp == s)
  moisture_factor <- (100-moisture_conversion$avg_MC)/100
  
  data_for_analysis <- biomass_data %>%
    select(Location, CropTrt, CropSp, Biomass_tha) %>%
    mutate(Drymass = Biomass_tha*moisture_factor) %>%
    filter(CropSp == s)
  
  test <- aov(Biomass_tha ~  Location + CropTrt + Location*CropTrt, data=data_for_analysis)
  print(s)
  print(summary(test))
  print(HSD.test(test, "CropTrt")$groups)
  print(HSD.test(test, "Location")$groups)
  
  test <- aov(Drymass ~  Location + CropTrt + Location*CropTrt, data=data_for_analysis)
  print(s)
  print(summary(test))
  print(HSD.test(test, "CropTrt")$groups)
  print(HSD.test(test, "Location")$groups)
}
sink()

ggplot(data_for_analysis, aes(Biomass_kg_9m2)) +
  geom_histogram(binwidth=5)

## Land Equivalent Ratio
data_for_LER <- biomass_data %>%
  select(Location, CropTrt, CropSp, LeavesStems_kg_1m2, Biomass_kg_9m2) %>%
  group_by(Location, CropTrt, CropSp) %>%
  summarize(avg_biomass_1m2 = mean(LeavesStems_kg_1m2),
            sd_biomass_1m2 = sd(LeavesStems_kg_1m2),
            avg_biomass_9m2 = mean(Biomass_kg_9m2),
            sd_biomass_9m2 = sd(Biomass_kg_9m2))

spp_combos <- data.frame(sp_1 = c("SH", "SH", "SS"), sp_2 = c("SS", "VB", "VB"))
LER_output <- data.frame(location = c(), sp_1 = c(), sp_2 = c(), 
                         LER_1m2 = c(), LER_9m2 = c())
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
    LER_1m2 <- round(mixed_1$avg_biomass_1m2 / mono_1$avg_biomass_1m2 + 
      mixed_2$avg_biomass_1m2 / mono_2$avg_biomass_1m2, 3)
    LER_9m2 <- round(mixed_1$avg_biomass_9m2 / mono_1$avg_biomass_9m2 + 
      mixed_2$avg_biomass_9m2 / mono_2$avg_biomass_9m2, 3)
    LER_row <- data.frame(l, c(spp_combos[i, ]), LER_1m2, LER_9m2)
    LER_output <- bind_rows(LER_output, LER_row)
  }
}
write_csv(LER_output, "output/LER-2018.csv")



## Regression for 1-9 m2

ggplot(biomass_data, aes(x=LeavesStems_tha, y=Biomass_tha,
                                    color=CropSp, shape=Location)) +
  geom_point(size=3)

ggplot(biomass_data, aes(x=LeavesStems_kg_1m2*9, y=Biomass_kg_9m2,
                                    color=CropTrt, shape=Location)) +
  geom_point(size=3) +
  geom_abline(slope=1)

