### Data Check and Analysis for the height and chlorophyll measurements on the
### cover crop experiment at TREC.

library(tidyverse)
library(agricolae)

## Import and Manage Data
height_week4 <- read_csv("data/height-06062018-TREC.csv") %>%
  separate(Height_mm, c("Height_mm", "undermeasure"), sep = "\\+") %>%
  mutate(Height_mm = as.numeric(Height_mm), CI = as.numeric(Chlorophyll))

avg_height_week4 <- height_week4 %>%
  filter(!is.na(Height_mm) | !is.na(CI)) %>% 
  group_by(SiteName, CropTrt, CropName, Position) %>% 
  summarize(avg_height = mean(Height_mm), 
            avg_CI = mean(CI))

height_winter2018 <- read_csv("data/height-winter-2018-TREC.csv") %>%
  separate(PlotName, c("SiteName", "CropTrt"), sep = "_") %>%
  separate(CropTrt, c("CropTrt", "Row"), sep = "-") %>%
  mutate(Height_mm = as.numeric(Height_mm), CI = Chlorophyll_CI)

avg_height_winter2018 <- height_winter2018 %>%
  filter(!is.na(Height_mm) | !is.na(CI)) %>% 
  group_by(SiteName, CropTrt, CropName, Row, SamplingWeeks) %>% 
  summarize(avg_height = mean(Height_mm), 
            avg_CI = mean(CI))

## Table Check

table_check_week4 <- height_week4 %>% 
  filter(!is.na(Height_mm) & !is.na(CI)) %>% 
  group_by(SiteName, CropTrt, CropName, Position) %>% 
  summarize(count = n())

table_check_week4_avg <- avg_height_week4 %>% 
  group_by(SiteName, CropTrt, CropName) %>% 
  summarize(count = n())

table_check_winter2018 <- height_winter2018 %>% 
  filter(!is.na(Height_mm) & !is.na(CI)) %>% 
  group_by(SiteName, CropTrt, CropName, Row, SamplingWeeks) %>% 
  summarize(count = n())

table_check_winter2018_avg <- avg_height_winter2018 %>%
  filter(!is.na(avg_height) & !is.na(avg_CI)) %>%
  group_by(SiteName, CropTrt, CropName) %>% 
  summarize(count = n())

## ANOVA for Height and CI
species <- c("SH", "SS", "VB")
sink("output/height-CI-anova-wk4.txt")
for(s in species){
  size_data <- filter(avg_height_week4, CropName == s)
  
  test <- aov(avg_height ~  SiteName + CropTrt + SiteName*CropTrt, 
              data=size_data)
  print(s)
  print(summary(test))
  print(HSD.test(test, "CropTrt")$groups)
  print(HSD.test(test, "SiteName")$groups)
  
  test <- aov(avg_CI ~  SiteName + CropTrt + SiteName*CropTrt, 
              data=size_data)
  print(s)
  print(summary(test))
  print(HSD.test(test, "CropTrt")$groups)
  print(HSD.test(test, "SiteName")$groups)
}
sink()

## Data Visualization

ggplot(avg_height_week4, aes(avg_height, fill = SiteName)) +
  geom_histogram(bins=20) +
  facet_grid(CropName~CropTrt)

ggplot(avg_height_week4, aes(x = avg_CI, fill = SiteName)) +
  geom_histogram(bins=20) +
  facet_grid(CropName~CropTrt)

species <- c("SH", "SS", "VB")
for (s in species){
  size_data <- filter(avg_height_winter2018, CropName == s)
  print(ggplot(size_data, aes(x = avg_height, fill = SiteName)) +
          geom_histogram(bins=20) +
          facet_grid(SamplingWeeks~CropTrt) + 
          labs(title = s))
}

species <- c("SH", "SS", "VB")
for (s in species){
  size_data <- filter(avg_height_winter2018, CropName == s)
  print(ggplot(size_data, aes(x = avg_CI, fill = SiteName)) +
    geom_histogram(bins=20) +
    facet_grid(SamplingWeeks~CropTrt) + 
    labs(title = s))
}



