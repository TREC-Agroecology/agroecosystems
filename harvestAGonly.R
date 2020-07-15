setwd("/Users/sreader/Documents/GitHub/agrocosystems")
library(tidyverse)
library(agricolae)
library(RColorBrewer)
library(dplyr)

biomass_s2018<-read_csv("data/harvest-092018-ECHO.csv")
biomass_s2018["season"] <- "summer"

##combining biomass data from 4 seasons (excluding summer 2018 which skews data)
biomass_s2019<- read_csv("data/harvest-072019-ECHO.csv")
biomass_s2019["season"] <- "summer"
biomass_w2019<-read_csv("data/harvest-012020-ECHO.csv")
biomass_w2019["season"]<- "winter"
biomass_w2018<-read_csv("data/harvest-012019-ECHO.csv")
biomass_w2018["season"]<-"winter"
biomass_s2020<- read_csv("data/harvest-072020-ECHO.csv")
biomass_s2020["season"]<-"summer"

biomass_data <- merge(biomass_s2019, biomass_s2020, all=TRUE)
biomass_data <- merge(biomass_data, biomass_w2018, all=TRUE)
biomass_data <- merge(biomass_data, biomass_w2019, all=TRUE)

seeding_rate_table <- data.frame(CropTrt = c("SH", "SS", "VB",
                                             "SHSS", "SHVB", "SSVB"),
                                 seeding_rate = c(1, 1, 1, 0.5, 0.5, 0.5))

biomass_data <- biomass_data %>%
  mutate(SiteName = str_sub(str_extract(biomass_data$PlotName, ".._"), end=2),
         CropTrt = str_sub(str_extract(biomass_data$PlotName, "_.*-"), 2,-2),
         LeavesStems_tha = LeavesStems_g_1m2*0.01) %>% 
  filter(CropSp == 'SS' | CropSp == 'SH' | CropSp == 'VB', SiteName== 'AG') %>% 
  left_join(seeding_rate_table)

total_biomass <- biomass_data %>%
  group_by(PlotID, season, CropTrt, CropSp, seeding_rate) %>%
  summarize(sum_LeavesStems_tha = sum(LeavesStems_tha),
            avg_LeavesStems_tha = mean(LeavesStems_tha),
            avg_StemCount = mean(StemCounts))

moisture_content <- read_csv("data/moisture-content.csv")
moisture_species <- moisture_content %>%
  group_by(CropSp) %>%
  summarize(avg_MC = mean(AGB_MC))

species <- c("SH", "SS", "VB")
treatments <- data.frame(treatments = c("SH", "SS", "VB", "SHSS", "SHVB", "SSVB"),
                         color_set = c("#A93226", "#D4AC0D", "#85C1E9", "#AF601A", "#8E44AD", "#196F3D"))

##summary table for biomass production
biomass_summary<- total_biomass %>% 
  group_by(season, CropTrt) %>% 
  summarize(total_biomass = sum(sum_LeavesStems_tha),
            mean_biomass = mean(avg_LeavesStems_tha))

##LER for CropTrt by season
data_for_LER <- total_biomass %>%
  select(CropTrt, CropSp, avg_LeavesStems_tha, season) %>%
  group_by(CropTrt, CropSp, season) %>%
  summarize(site_avg_biomass = mean(avg_LeavesStems_tha),
            sd_biomass = sd(avg_LeavesStems_tha))

spp_combos <- data.frame(sp_1 = c("SH", "SH", "SS"), sp_2 = c("SS", "VB", "VB"))
LER_output <- data.frame(season = c(), sp_1 = c(), sp_2 = c(), LER = c())
for (j in unique(data_for_LER$season)){
  season_data <- data_for_LER %>%
    filter(season == j)
  for (i in 1:nrow(spp_combos)){
    mixed_1 <- season_data %>%
      filter(CropTrt == paste(spp_combos$sp_1[i], spp_combos$sp_2[i], sep="") &
               CropSp == spp_combos$sp_1[i])
    mono_1 <- season_data %>%
      filter(CropTrt == spp_combos$sp_1[i] & CropSp == spp_combos$sp_1[i])
    mixed_2 <- season_data %>%
      filter(CropTrt == paste(spp_combos$sp_1[i], spp_combos$sp_2[i], sep="") &
               CropSp == spp_combos$sp_2[i])
    mono_2 <- season_data %>%
      filter(CropTrt == spp_combos$sp_2[i] & CropSp == spp_combos$sp_2[i])
    LER <- round(mixed_1$site_avg_biomass / mono_1$site_avg_biomass + 
                   mixed_2$site_avg_biomass / mono_2$site_avg_biomass, 3)
    LER_row <- data.frame(j, c(spp_combos[i, ]), LER)
    LER_output <- bind_rows(LER_output, LER_row)
  }
}

write_csv(LER_output, "data/output/LER-ECHO.csv")

##Calculating LER by Row
data_for_LER_row <- biomass_data %>%
  group_by(PlotID, CropSp, CropTrt, RowNum, season) %>%
  summarize(avg_LeavesStems_tha = mean(LeavesStems_tha),
            avg_StemCount = mean(StemCounts)) %>%
  group_by(CropTrt, CropSp, RowNum, season) %>%
  summarize(site_avg_biomass = mean(avg_LeavesStems_tha))

LER_output_row <- data.frame(season = c(), col = c(), sp_1 = c(), sp_2 = c(), LER_1 = c(), LER_2 = c(), LER_3 = c())
for (j in unique(data_for_LER$season)){
  season_data <- data_for_LER_row %>%
    filter(season == j)
  for(r in 1:6){
    season_data %>%
      filter(RowNum == r)
    for (i in 1:nrow(spp_combos)){
      mixed_1 <- season_data %>%
        filter(CropTrt == paste(spp_combos$sp_1[i], spp_combos$sp_2[i], sep="") &
                 CropSp == spp_combos$sp_1[i])
      mono_1 <- season_data %>%
        filter(CropTrt == spp_combos$sp_1[i] & CropSp == spp_combos$sp_1[i])
      mixed_2 <- season_data %>%
        filter(CropTrt == paste(spp_combos$sp_1[i], spp_combos$sp_2[i], sep="") &
                 CropSp == spp_combos$sp_2[i])
      mono_2 <- season_data %>%
        filter(CropTrt == spp_combos$sp_2[i] & CropSp == spp_combos$sp_2[i])
      LER_1 <- mixed_1$site_avg_biomass / mono_1$site_avg_biomass
      LER_2 <- mixed_2$site_avg_biomass / mono_2$site_avg_biomass
      LER_3 <- round(mixed_1$site_avg_biomass / mono_1$site_avg_biomass + 
                       mixed_2$site_avg_biomass / mono_2$site_avg_biomass, 3)
      LER_record <- data.frame(j, r, c(spp_combos[i, ]), LER_1, LER_2, LER_3)
      LER_output_row <- bind_rows(LER_output_row, LER_record)
    }
  }
}
LER_row_summary <- LER_output_row %>% 
  mutate(j = factor(j, levels = c("summer", "winter"))) %>% 
  group_by(j, sp_1, sp_2) %>% 
  summarize(mean_LER = mean(LER_3), ci_LER = 2*sd(LER_3)/sqrt(6), mean_LERfract_sp1 = mean(LER_1), ci_LERfract_sp1 = 2*sd(LER_1)/sqrt(6), mean_LERfract_sp2 = mean(LER_2), ci_LERfract_sp2 = 2*sd(LER_2)/sqrt(6)) %>% 
  mutate(spp = paste0(sp_1, "-", sp_2))

ggplot(LER_row_summary, aes(x=spp, y=mean_LER)) +
  geom_bar(stat="identity", aes(fill=spp)) +
  geom_errorbar(aes(ymin=mean_LER-ci_LER, ymax=mean_LER+ci_LER), width=0.2) +
  facet_grid(.~j) +
  geom_hline(yintercept = 1, linetype="dashed") +
  scale_fill_manual(values=c("#AF601A", "#8E44AD", "#196F3D")) +
  labs(x="Crop Mix", y="LER [+/- 95% CI]") +
  theme_bw(base_size = 24, base_family = "Helvetica") +
  theme(axis.text.x = element_text(size=14, angle=30),
        legend.position = "none")
ggsave("output/LER.png")

#visualizing CI for LER Factors **need to reorganize table to fract all in one column
LER_output_row_fract <- data.frame(season = c(), row = c(), sp = c(), croptrt = c(), LER_fract = c())
for (j in unique(data_for_LER$season)){
  season_data <- data_for_LER_row %>%
    filter(season == j)
  for(r in 1:6){
    season_data %>%
      filter(RowNum == r)
    for (i in 1:nrow(spp_combos)){
      mixed_1 <- season_data %>%
        filter(CropTrt == paste(spp_combos$sp_1[i], spp_combos$sp_2[i], sep="") &
                 CropSp == spp_combos$sp_1[i])
      mono_1 <- season_data %>%
        filter(CropTrt == spp_combos$sp_1[i] & CropSp == spp_combos$sp_1[i])
      mixed_2 <- season_data %>%
        filter(CropTrt == paste(spp_combos$sp_1[i], spp_combos$sp_2[i], sep="") &
                 CropSp == spp_combos$sp_2[i])
      mono_2 <- season_data %>%
        filter(CropTrt == spp_combos$sp_2[i] & CropSp == spp_combos$sp_2[i])
      LER_fract_sp1 <- mixed_1$site_avg_biomass / mono_1$site_avg_biomass
      LER_fract_sp2 <- mixed_2$site_avg_biomass / mono_2$site_avg_biomass
      LER_fract <-cbind(LER_fract_sp_1, LER_fract_sp2)
      LER_record_fract <- data.frame(j, r, species, spp_combos[i, ], LER_fract)
      LER_output_row_fract <- bind_rows(LER_output_row_fract, LER_record_fract)
    }
  }
}
LER_row_summary_fract <- LER_output_row_fract %>% 
  mutate(j = factor(j, levels = c("summer", "winter"))) %>% 
  group_by(j, CropSp) %>% 
  summarize(mean_LER = mean(LER_3), ci_LER = 2*sd(LER_3)/sqrt(6), mean_LERfract_sp1 = mean(LER_1), ci_LERfract_sp1 = 2*sd(LER_1)/sqrt(6), mean_LERfract_sp2 = mean(LER_2), ci_LERfract_sp2 = 2*sd(LER_2)/sqrt(6)) %>% 
  mutate(spp = paste0(sp_1, "-", sp_2))

ggplot(LER_row_summary, aes(x=CropTrt, y=mean_LERfract, group = CropSp)) +
  geom_bar(stat="identity", aes(fill=spp)) +
  geom_errorbar(aes(ymin=mean_LERfract_sp1-ci_LERfract_sp1, ymax=mean_LERfract_sp1+ci_LERfract_sp1), width=0.2) +
  facet_grid(CropSp~j) +
  geom_hline(yintercept = 0.5, linetype="dashed") +
  scale_fill_manual(values=c("#AF601A", "#8E44AD", "#196F3D")) +
  labs(x="Crop Mix", y="LERfraction [+/- 95% CI]") +
  theme_bw(base_size = 24, base_family = "Helvetica") +
  theme(axis.text.x = element_text(size=14, angle=30),
        legend.position = "none")

##Calculating LER by Column
data_for_LER_column <- biomass_data %>%
  group_by(PlotID, CropSp, CropTrt, ColNum, season) %>%
  summarize(avg_LeavesStems_tha = mean(LeavesStems_tha),
            avg_StemCount = mean(StemCounts)) %>%
  group_by(CropTrt, CropSp, ColNum, season) %>%
  summarize(site_avg_biomass = mean(avg_LeavesStems_tha))

LER_output_column <- data.frame(season = c(), col = c(), sp_1 = c(), sp_2 = c(), LER_1 = c(), LER_2 = c(), LER_3 = c())
for (j in unique(data_for_LER$season)){
  season_data <- data_for_LER_column %>%
    filter(season == j)
  for(z in 1:6){
    season_data %>%
      filter(ColNum == z)
    for (i in 1:nrow(spp_combos)){
      mixed_1 <- season_data %>%
        filter(CropTrt == paste(spp_combos$sp_1[i], spp_combos$sp_2[i], sep="") &
                 CropSp == spp_combos$sp_1[i])
      mono_1 <- season_data %>%
        filter(CropTrt == spp_combos$sp_1[i] & CropSp == spp_combos$sp_1[i])
      mixed_2 <- season_data %>%
        filter(CropTrt == paste(spp_combos$sp_1[i], spp_combos$sp_2[i], sep="") &
                 CropSp == spp_combos$sp_2[i])
      mono_2 <- season_data %>%
        filter(CropTrt == spp_combos$sp_2[i] & CropSp == spp_combos$sp_2[i])
      LER_1 <- mixed_1$site_avg_biomass / mono_1$site_avg_biomass
      LER_2 <- mixed_2$site_avg_biomass / mono_2$site_avg_biomass
      LER_3 <- round(mixed_1$site_avg_biomass / mono_1$site_avg_biomass + 
                     mixed_2$site_avg_biomass / mono_2$site_avg_biomass, 3)
      LER_record <- data.frame(j, z, c(spp_combos[i, ]), LER_1, LER_2, LER_3)
      LER_output_column <- bind_rows(LER_output_column, LER_record)
    }
  }
}
LER_column_summary <- LER_output_column %>% 
  mutate(j = factor(j, levels = c("summer", "winter"))) %>% 
  group_by(j, sp_1, sp_2) %>% 
  summarize(mean_LER = mean(LER_3), ci_LER = 2*sd(LER_3)/sqrt(6), mean_LERfract_sp1 = mean(LER_1), ci_LERfract_sp1 = 2*sd(LER_1)/sqrt(6), mean_LERfract_sp2 = mean(LER_2), ci_LERfract_sp2 = 2*sd(LER_2)/sqrt(6)) %>% 
  mutate(spp = paste0(sp_1, "-", sp_2))

write_csv(LER_column_summary, "data/output/LER-ECHO-column-summary.csv")

ggplot(LER_column_summary, aes(x=spp, y=mean_LER)) +
  geom_bar(stat="identity", aes(fill=spp)) +
  geom_errorbar(aes(ymin=mean_LER-ci_LER, ymax=mean_LER+ci_LER), width=0.2) +
  facet_grid(.~j) +
  geom_hline(yintercept = 1, linetype="dashed") +
  scale_fill_manual(values=c("#AF601A", "#8E44AD", "#196F3D")) +
  labs(x="Crop Mix", y="LER [+/- 95% CI]") +
  theme_bw(base_size = 24, base_family = "Helvetica") +
  theme(axis.text.x = element_text(size=14, angle=30),
        legend.position = "none")
ggsave("output/LER.png")

##confidence interval values for graphs (row and column)
LER_row_ci <- LER_output_row %>% 
  mutate(j = factor(j, levels = c(1,2,3,4))) %>% 
  group_by(j, sp_1, sp_2) %>% 
  summarize(mean_LER = mean(LER), ci_LER = 2*sd(LER)/sqrt(6))

LER_column_ci <- LER_output_column %>% 
  mutate(j = factor(j, levels = c(1,2,3,4))) %>% 
  group_by(j, sp_1, sp_2) %>% 
  summarize(mean_LER = mean(LER), ci_LER = 2*sd(LER)/sqrt(6))

## Biomass visualization per species
for(s in species){
  data_for_figure <- total_biomass %>% 
    ungroup() %>%
    filter(CropSp == s) %>%
    mutate(CropTrt = factor(CropTrt, levels = str_subset(treatments$treatments, s))) %>%
    group_by(season, CropTrt) %>%
    summarize(avg_stem_mass = mean(avg_LeavesStems_tha/(0.01*avg_StemCount), na.rm=TRUE),
              ci_stem_mass = 2*sd(avg_LeavesStems_tha/(0.01*avg_StemCount), na.rm=TRUE)/sqrt(n())) %>% 
    left_join(treatments, by = c("CropTrt" = "treatments")) %>% 
    mutate(CropTrt = factor(CropTrt, levels = str_subset(treatments$treatments, s)))
  
  ggplot(data_for_figure, aes(x=CropTrt, y=avg_stem_mass)) +
    geom_bar(stat="identity", aes(fill=CropTrt)) +
    geom_errorbar(aes(ymin = avg_stem_mass-ci_stem_mass, ymax = avg_stem_mass+ci_stem_mass), width=0.2) +
    facet_grid(.~season) +
    labs(x="Crop Mix", y="Fresh Mass / Individual [g +/- 95% CI]", title = paste0(s)) +
    scale_fill_manual(values=as.character(unique(data_for_figure$color_set))) +
    theme_bw(base_size = 20, base_family = "Helvetica") +
    theme(legend.position="none")
}

## ANOVA for total biomass per treatment

sink("output/CombinedBiomass_anova.txt")

test <- aov(sum_LeavesStems_tha ~  season + CropTrt + season*CropTrt, data=combined_biomass)
print(summary(test))
print(HSD.test(test, "season")$groups)
print(HSD.test(test, "CropTrt")$groups)

test <- aov(sum_LeavesStems_tha ~  CropTrt + Location_CropTrt, data=combined_biomass)
print(summary(test))
print(HSD.test(test, "CropTrt")$groups)
print(HSD.test(test, "Location_CropTrt")$groups)

sink()

## Raster for total biomass per treatment

combined_biomass <- total_biomass %>%
  group_by(PlotID, season, CropTrt) %>%
  summarize(sum_LeavesStems_tha = sum(avg_LeavesStems_tha))

global_avg_combined_biomass <- mean(log(combined_biomass$sum_LeavesStems_tha))
global_sd_combined_biomass <- sd(log(combined_biomass$sum_LeavesStems_tha))

combined_biomass_std <- combined_biomass %>%
  mutate(sum_biomass_std = ((log(sum_LeavesStems_tha)-global_avg_combined_biomass))/
           global_sd_combined_biomass)

png(file="output/combined_std.png", width= 1700, height=1100)
matrix_plot(combined_biomass_std, "+/- Global SD")
dev.off()

trt_avg_combined_biomass <- combined_biomass %>% 
  group_by(CropTrt) %>% 
  summarize(mean_biomass = mean(log(sum_LeavesStems_tha)),
            sd_biomass = sd(log(sum_LeavesStems_tha)))

combined_biomass_std_trt <- combined_biomass %>% 
  mutate(sum_biomass_std = treatment_standardization(CropTrt, trt_avg_combined_biomass, sum_LeavesStems_tha))

png(file="output/combined_std_trt.png", width= 1700, height=1100)
matrix_plot(combined_biomass_std_trt, "+/- Treatment SD")
dev.off()

## ANOVA for average stem biomass per species.
sink("output/Biomass_anova.txt")

for(s in species){
  moisture_conversion <- moisture_species %>%
    filter(CropSp == s)
  moisture_factor <- (100-moisture_conversion$avg_MC)/100
  
  data_for_analysis <- total_biomass %>%
    select(PlotID, CropSp, CropTrt, season, 
           avg_LeavesStems_tha, avg_StemCount, seeding_rate) %>%
    mutate(avg_drymass = avg_LeavesStems_tha*moisture_factor) %>%
    filter(CropSp == s)
  
  test <- aov(avg_LeavesStems_tha ~  CropTrt + season, data=data_for_analysis)
  print(paste(s, "- Fresh Mass [t/ha]"))
  print(summary(test))
  print(HSD.test(test, "CropTrt")$groups)
  print(HSD.test(test, "season")$groups)
  
  test <- aov(avg_drymass ~  CropTrt + season, data=data_for_analysis)
  print(paste(s, "- Dry Mass"))
  print(summary(test))
  print(HSD.test(test, "CropTrt")$groups)
  print(HSD.test(test, "season")$groups)
  
  test <- aov(avg_LeavesStems_tha/seeding_rate ~  CropTrt + season, data=data_for_analysis)
  print(paste(s, "- Fresh Mass / Seeding Rate"))
  print(summary(test))
  print(HSD.test(test, "CropTrt")$groups)
  print(HSD.test(test, "season")$groups)
  
  test <- aov(avg_LeavesStems_tha/(0.01*avg_StemCount) ~  CropTrt + season, data=data_for_analysis)
  print(paste(s, "- Fresh Mass / Stem Count [g/ind]"))
  print(summary(test))
  print(HSD.test(test, "CropTrt")$groups)
  print(HSD.test(test, "season")$groups)
  
  test <- aov(avg_LeavesStems_tha/(0.01*avg_StemCount) ~  season + CropTrt + season*CropTrt, data=data_for_analysis)
  print(paste(s, "- Fresh Mass / Stem Count [g/ind]"))
  print(summary(test))
  print(HSD.test(test, "season")$groups)
  print(HSD.test(test, "CropTrt")$groups)
  
}
sink()

for(s in species){
  data_for_figure <- total_biomass %>% 
    ungroup() %>%
    filter(CropSp == s) %>%
    mutate(CropTrt = factor(CropTrt, levels = str_subset(treatments$treatments, s))) %>%
    group_by(season, CropTrt) %>%
    summarize(avg_stem_mass = mean(avg_LeavesStems_tha/(0.01*avg_StemCount), na.rm=TRUE),
              ci_stem_mass = 2*sd(avg_LeavesStems_tha/(0.01*avg_StemCount), na.rm=TRUE)/sqrt(n())) %>% 
    left_join(treatments, by = c("CropTrt" = "treatments")) %>% 
    mutate(CropTrt = factor(CropTrt, levels = str_subset(treatments$treatments, s)))
  
  ggplot(data_for_figure, aes(x=CropTrt, y=avg_stem_mass)) +
    geom_bar(stat="identity", aes(fill=CropTrt)) +
    geom_errorbar(aes(ymin = avg_stem_mass-ci_stem_mass, ymax = avg_stem_mass+ci_stem_mass), width=0.2) +
    facet_grid(.~Location) +
    labs(x="Crop Mix", y="Fresh Mass / Individual [g +/- 95% CI]", title = paste0(s)) +
    scale_fill_manual(values=as.character(unique(data_for_figure$color_set))) +
    theme_bw(base_size = 20, base_family = "Helvetica") +
    theme(legend.position="none")
  ggsave(paste0("output/Stem_mass_", s, ".png"))
}

##T-tests against the value 1
ttest_data<- LER_output_row %>% 
  mutate(spp = paste0(sp_1, "-", sp_2))

ttest_data_SHSS<- ttest_data %>%
  filter(spp == "SH-SS")

t.test(ttest_data_SHSS$LER_3,mu=1, alternative='greater')

ttest_data_SHVB<- ttest_data %>%
  filter(spp == "SH-VB")

t.test(ttest_data_SHVB$LER_3,mu=1, alternative='greater')

ttest_data_SSVB<- ttest_data %>%
  filter(spp == "SS-VB")

t.test(ttest_data_SSVB$LER_3,mu=1, alternative='greater')


##T-tests against the value 1 (from summary)
ttest_data<- LER_row_summary

ttest_data_SHSS<- ttest_data %>%
  filter(spp == "SH-SS")

t.test(ttest_data_SHSS$mean_LER,mu=1, alternative='greater')

ttest_data_SHVB<- ttest_data %>%
  filter(spp == "SH-VB")

t.test(ttest_data_SHVB$mean_LER,mu=1, alternative='greater')


ttest_data_SSVB<- ttest_data %>%
  filter(spp == "SS-VB")

t.test(ttest_data_SSVB$mean_LER,mu=1, alternative='greater')
