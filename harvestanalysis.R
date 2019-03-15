### Thioro and Zack are starting the file to analyse the cover crop data.
### This script requires set working directory to source file location.

### To Do:
### - Calculate LER

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


species <- c("SH", "SS", "VB")

## ANOVA for TotalBiomass
sink("output/Biomass_9m2_anova.txt")
for(s in species){
  data_for_analysis <- biomass_data %>%
    select(Location, CropTrt, CropSp, Biomass_kg_9m2) %>%
    filter(CropSp == s)
  
  test <- aov(log(Biomass_kg_9m2) ~  Location + CropTrt + Location*CropTrt, data=data_for_analysis)
  print(s)
  print(summary(test))
  print(HSD.test(test, "CropTrt")$groups)
  print(HSD.test(test, "Location")$groups)
}
sink()

ggplot(data_for_analysis, aes(Biomass_kg_9m2)) +
  geom_histogram(binwidth=5)

## Regression for 1-9 m2

ggplot(biomass_data, aes(x=LeavesStems_tha, y=Biomass_tha,
                                    color=CropSp, shape=Location)) +
  geom_point(size=3)

ggplot(total_biomass, aes(x=LeavesStems_kg_1m2, y=Biomass_kg_9m2,
                                    color=CropTrt)) +
  geom_point(size=3) +
  facet_grid(. ~ Location)

