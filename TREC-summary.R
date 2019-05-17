### This file collects data products for each cover crop summarized 
### per experimental unit with treatments described.

library(tidyverse)

crops <- c("SH", "SS", "VB")

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

for (crop in crops){
biomass_soil_summer_crop <- left_join(biomass_summer, soil_OM_summer) %>%
  filter(CropSp == crop) %>%
  group_by(Location, CropTrt) %>%
  summarize(avg_LeavesStems_tha = round(mean(LeavesStems_tha), 2), 
            sd_LeavesStems_tha = round(sd(LeavesStems_tha), 2),
            avg_SOM_percent = round(mean(SOM_percent), 2),
            sd_SOM_percent = round(sd(SOM_percent), 2),
            avg_TP_mgkg = round(mean(TP_mgkg), 2),
            sd_TP_mgkg = round(sd(TP_mgkg), 2))
write_csv(biomass_soil_summer_crop,
          paste0("data/data-products/biomass_soil_summer_", crop, ".csv"))

height_summer_crop <- height_summer %>%
  filter(CropSp == crop) %>% 
  group_by(SamplingWeek, Location, CropTrt) %>% 
  summarize(avg_Height_mm = round(mean(Height_mm), 2), 
            sd_Height_mm = round(sd(Height_mm), 2))
write_csv(height_summer_crop,
          paste0("data/data-products/height_summer_", crop, ".csv"))
}

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
  group_by(SamplingWeek, PlotID, Location=SiteName, CropTrt, CropSp=CropName, Row) %>% 
  summarize(Height_mm = round(mean(as.numeric(Height_mm), na.rm=TRUE), 2), 
            Chlorophyll_CI = round(mean(Chlorophyll_CI), 2))

soil_OM_winter <- read_csv("data/soil-SOM-TP-2018.csv") %>%
  filter(Month == "Nov") %>%
  select(PlotID, SOM_percent)

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

for (crop in crops){
  biomass_soil_winter_crop <- left_join(biomass_winter, soil_OM_winter) %>%
    left_join(soil_BD_winter) %>%
    left_join(soil_N_winter) %>%
    filter(CropSp == crop) %>%
    group_by(Location, CropTrt) %>%
    summarize(avg_LeavesStems_tha = round(mean(LeavesStems_tha), 2), 
              sd_LeavesStems_tha = round(sd(LeavesStems_tha), 2),
              avg_SOM_percent = round(mean(SOM_percent), 2),
              sd_SOM_percent = round(sd(SOM_percent), 2),
              avg_GravelDryWgt_g = round(mean(GravelDryWgt_g), 2),
              sd_GravelDryWgt_g = round(sd(GravelDryWgt_g), 2),
              avg_GravelVol_mL = round(mean(GravelVol_mL), 2),
              sd_GravelVol_mL = round(sd(GravelVol_mL), 2),
              avg_SoilDryWgt_g = round(mean(SoilDryWgt_g), 2),
              sd_SoilDryWgt_g = round(sd(SoilDryWgt_g), 2),
              avg_SoilVol_cm3 = round(mean(SoilVol_cm3), 2),
              sd_SoilVol_cm3 = round(sd(SoilVol_cm3), 2),
              avg_SoilBD_gcm3 = round(mean(SoilBD_gcm3), 2),
              sd_SoilBD_gcm3 = round(sd(SoilBD_gcm3), 2),
              avg_SoilPorosity_percent = round(mean(SoilPorosity_percent), 2),
              sd_SoilPorosity_percent = round(sd(SoilPorosity_percent), 2),
              avg_NH4_mgL = round(mean(NH4_mgL), 2),
              sd_NH4_mgL = round(sd(NH4_mgL), 2),
              avg_NO3_mgL = round(mean(NO3_mgL), 2),
              sd_NO3_mgL = round(sd(NO3_mgL), 2),
              avg_NH4_mgKg = round(mean(NH4_mgKg), 2),
              sd_NH4_mgKg = round(sd(NH4_mgKg), 2),
              avg_NO3_mgKg = round(mean(NO3_mgKg), 2),
              sd_NO3_mgKg = round(sd(NO3_mgKg), 2))
  write_csv(biomass_soil_winter_crop,
            paste0("data/data-products/biomass_soil_winter_", crop, ".csv"))
  
  height_winter_crop <- height_winter %>%
    filter(CropSp == crop) %>% 
    group_by(SamplingWeek, Location, CropTrt) %>% 
    summarize(avg_Height_mm = round(mean(Height_mm), 2), 
              sd_Height_mm = round(sd(Height_mm), 2),
              avg_Chlorophyll_CI = round(mean(Chlorophyll_CI), 2),
              sd_Chlorophyll_CI = round(sd(Chlorophyll_CI), 2))
  write_csv(height_winter_crop,
            paste0("data/data-products/height_winter_", crop, ".csv"))
}


