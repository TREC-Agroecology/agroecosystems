### This file collects data products for each cover crop summarized 
### per experimental unit with treatments described.

library(tidyverse)

## Summer Data

biomass_summer <- read_csv("data/harvest-092018-TREC.csv") %>%
  mutate(LeavesStems_tha = round(Biomass_kg_9m2*10/9, 2)) %>%
  select(PlotID, Location, CropTrt, CropSp, Row=RowNum, LeavesStems_tha)

height_summer <- read_csv("data/height-summer-2018-TREC.csv") %>%
  group_by(SamplingWeek, PlotID, SiteName, CropTrt, CropName, RowNum) %>% 
  summarize(Height_mm = round(mean(as.numeric(Height_mm), na.rm=TRUE), 0)) %>%
  rename(Location = SiteName, CropSp = CropName)

soil_OM_summer <- read_csv("data/soil-SOM-TP-2018.csv") %>%
  filter(Month == "May") %>%
  select(PlotID, SOM_percent, TP_mgkg)

biomass_soil_summer_SH <- left_join(biomass_summer, soil_OM_summer) %>%
  filter(CropSp == "SH") %>%
  group_by(Location, CropTrt) %>%
  summarize(avg_LeavesStems_tha = mean(LeavesStems_tha), 
            sd_LeavesStems_tha_sd = sd(LeavesStems_tha),
            avg_SOM_percent = mean(SOM_percent),
            sd_SOM_percent = sd(SOM_percent),
            avg_TP_mgkg = mean(TP_mgkg),
            sd_TP_mgkg = sd(TP_mgkg))

height_summer_SH <- height_summer %>%
  filter(CropSp == "SH") %>% 
  group_by(SamplingWeek, Location, CropTrt) %>% 
  summarize(avg_Height_mm = round(mean(Height_mm), 2), 
            sd_Height_mm = sd(Height_mm))

## Winter Data

biomass_winter <- read_csv("data/harvest-012019-TREC.csv") %>%
  separate(PlotName, c("Location", "CropTrt"), sep = "_") %>%
  separate(CropTrt, c("CropTrt", "Row"), sep = "-") %>%
  mutate(LeavesStems_tha = LeavesStems_g_1m2*0.01) %>%
  group_by(PlotID, Location, CropTrt, CropSp, Row) %>%
  summarize(LeavesStems_tha = round(mean(LeavesStems_tha),2))

height_winter <- read_csv("data/height-winter-2018-TREC.csv") %>%
  separate(PlotName, c("SiteName", "CropTrt"), sep = "_") %>%
  separate(CropTrt, c("CropTrt", "Row"), sep = "-") %>%
  group_by(SamplingWeek, PlotID, SiteName, CropTrt, CropName, Row) %>% 
  summarize(Height_mm = round(mean(as.numeric(Height_mm), na.rm=TRUE), 2), 
            Chlorophyll_CI = round(mean(Chlorophyll_CI), 2))

soil_BD_winter <- read_csv("data/soil-BD-092018.csv") %>%
  separate(PlotName, c("SiteName", "CropTrt"), sep = "_") %>%
  separate(CropTrt, c("CropTrt", "Row"), sep = "-")

soil_N_winter <- left_join(
  read_csv("data/soil-N-092018.csv") %>%
    select(PlotID, NitrogenForms, Concentration_mgL) %>%
    spread(NitrogenForms, Concentration_mgL) %>%
    rename(NH4_mgL = `NH4-N`, NO3_mgL = `NO3-N`),
  
  read_csv("data/soil-N-092018.csv") %>%
    select(PlotID, NitrogenForms, Concentration_mgKg) %>%
    spread(NitrogenForms, Concentration_mgKg) %>%
    rename(NH4_mgKg = `NH4-N`, NO3_mgKg = `NO3-N`)
)



