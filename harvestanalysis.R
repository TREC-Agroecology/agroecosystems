### Thioro and Zack are starting the file to analyse the cover crop data.
### This script requires set working directory to source file location.

library(tidyverse)
library(agricolae)
library(raster)
library(RColorBrewer)

matrix_plot <- function(data_for_matrix, title){
  plotid_matrix <- matrix(data = c(1:144), nrow=24, ncol=6, byrow=TRUE)
  biomass_std_matrix <- matrix(nrow=24, ncol=6)
  biomass_std_names <- matrix(nrow=24, ncol=6)
  
  for (i in 1:24){
    for (j in 1:6){
      plot_id <- plotid_matrix[i,j]
      biomass_std_matrix[i,j] <- data_for_matrix$sum_biomass_std[plot_id]
      biomass_std_names[i,j] <- data_for_matrix$CropTrt[plot_id]
    }
  }
  
  par(mfrow=c(2,2), mar=c(1,0,1,0))
  plot(raster(biomass_std_matrix[1:6,]), axes=FALSE, box=FALSE, legend=FALSE,
       zlim=c(min(data_for_matrix$sum_biomass_std), max(data_for_matrix$sum_biomass_std)),
       breaks= c(3, 2, 1, 0, -1, -2, -3), col= brewer.pal(7, "RdBu"))
  text(x=rep(c(0.1, 0.3, 0.5, 0.65, 0.8, 1),6), 
       y=rep(c(1, 5/6, 4/6, 3/6, 2/6, 1/6),each=6),
       label=data_for_matrix$CropTrt[1:36], adj=c(1,1))
  plot(raster(biomass_std_matrix[7:12,]), axes=FALSE, box=FALSE, legend=FALSE,
       zlim=c(min(data_for_matrix$sum_biomass_std), max(data_for_matrix$sum_biomass_std)),
       breaks= c(3, 2, 1, 0, -1, -2, -3), col= brewer.pal(7, "RdBu"))
  text(x=rep(c(0.1, 0.3, 0.5, 0.65, 0.8, 1),6), 
       y=rep(c(1, 5/6, 4/6, 3/6, 2/6, 1/6),each=6),
       label=data_for_matrix$CropTrt[36:72], adj=c(1,1))
  plot(raster(biomass_std_matrix[13:18,]), axes=FALSE, box=FALSE, legend=FALSE, 
       zlim=c(min(data_for_matrix$sum_biomass_std), max(data_for_matrix$sum_biomass_std)),
       breaks= c(3, 2, 1, 0, -1, -2, -3), col= brewer.pal(7, "RdBu"))
  text(x=rep(c(0.1, 0.3, 0.5, 0.65, 0.8, 1),6), 
       y=rep(c(1, 5/6, 4/6, 3/6, 2/6, 1/6),each=6),
       label=data_for_matrix$CropTrt[73:108], adj=c(1,1))
  plot(raster(biomass_std_matrix[19:24,]), axes=FALSE,  box=FALSE,
       zlim=c(min(data_for_matrix$sum_biomass_std), max(data_for_matrix$sum_biomass_std)),
       breaks= c(3, 2, 1, 0, -1, -2, -3), col= brewer.pal(7, "RdBu"),
       legend.width=4,  legend.args=list(text=title, side=2, cex=2))
  text(x=rep(c(0.1, 0.3, 0.5, 0.65, 0.8, 1),6), 
       y=rep(c(1, 5/6, 4/6, 3/6, 2/6, 1/6),each=6),
       label=data_for_matrix$CropTrt[109:144], adj=c(1,1))
}

treatment_standardization <- function(trt, std_table, mass_value) {
  trt_record <- filter(std_table, CropTrt == trt)
  std_value <- (log(mass_value) - trt_record$mean_biomass)/ trt_record$sd_biomass
  return(std_value)
}

## Import and Manage Data

seeding_rate_table <- data.frame(CropTrt = c("SH", "SS", "VB",
                                             "SHSS", "SHVB", "SSVB"),
                                 seeding_rate = c(1, 1, 1, 0.5, 0.5, 0.5))

biomass_data <- read_csv("data/harvest-012019-TREC.csv") %>% 
  separate(PlotName, c("Location_CropTrt", "Row"), sep = "-") %>%
  separate(Location_CropTrt, c("Location", "CropTrt"), sep = "_", remove=FALSE) %>%
  mutate(Location = factor(Location, levels = c("VF", "FO", "OC", "HH")),
         CropTrt = factor(CropTrt, levels = c("SH", "SS", "VB", "SHSS", "SHVB", "SSVB")),
         LeavesStems_tha = LeavesStems_g_1m2*0.01) %>%
  left_join(seeding_rate_table)

total_biomass <- biomass_data %>%
  group_by(PlotID, CropSp, Location, CropTrt, Location_CropTrt, seeding_rate) %>%
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
treatments <- c("SH", "SS", "VB", "SHSS", "SHVB", "SSVB")


## ANOVA for total biomass per treatment

sink("output/CombinedBiomass_anova.txt")
  
test <- aov(sum_LeavesStems_tha ~  Location + CropTrt + Location*CropTrt, data=combined_biomass)
print(summary(test))
print(HSD.test(test, "Location")$groups)
print(HSD.test(test, "CropTrt")$groups)

test <- aov(sum_LeavesStems_tha ~  CropTrt + Location_CropTrt, data=combined_biomass)
print(summary(test))
print(HSD.test(test, "CropTrt")$groups)
print(HSD.test(test, "Location_CropTrt")$groups)

sink()

## Raster for total biomass per treatment

combined_biomass <- total_biomass %>%
  group_by(PlotID, Location, CropTrt, Location_CropTrt) %>%
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
    select(PlotID, CropSp, Location, CropTrt, Location_CropTrt, 
           avg_LeavesStems_tha, avg_StemCount, seeding_rate) %>%
    mutate(avg_drymass = avg_LeavesStems_tha*moisture_factor) %>%
    filter(CropSp == s)
  
  test <- aov(avg_LeavesStems_tha ~  CropTrt + Location_CropTrt, data=data_for_analysis)
  print(paste(s, "- Fresh Mass [t/ha]"))
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
  
  test <- aov(avg_LeavesStems_tha/(0.01*avg_StemCount) ~  CropTrt + Location_CropTrt, data=data_for_analysis)
  print(paste(s, "- Fresh Mass / Stem Count [g/ind]"))
  print(summary(test))
  print(HSD.test(test, "CropTrt")$groups)
  print(HSD.test(test, "Location_CropTrt")$groups)
  
  test <- aov(avg_LeavesStems_tha/(0.01*avg_StemCount) ~  Location + CropTrt + Location*CropTrt, data=data_for_analysis)
  print(paste(s, "- Fresh Mass / Stem Count [g/ind]"))
  print(summary(test))
  print(HSD.test(test, "Location")$groups)
  print(HSD.test(test, "CropTrt")$groups)
  
}
sink()

for(s in species){
  data_for_figure <- total_biomass %>% 
    ungroup() %>%
    filter(CropSp == s) %>%
    mutate(CropTrt = factor(CropTrt, levels = str_subset(treatments, s))) %>% 
    group_by(Location, CropTrt) %>%
    summarize(avg_stem_mass = mean(avg_LeavesStems_tha/(0.01*avg_StemCount), na.rm=TRUE),
              se_stem_mass = sd(avg_LeavesStems_tha/(0.01*avg_StemCount*sqrt(n())), na.rm=TRUE))

  ggplot(data_for_figure, aes(x=factor(CropTrt), y=avg_stem_mass)) +
    geom_bar(stat="identity") +
    geom_errorbar(aes(ymin = avg_stem_mass-se_stem_mass, ymax = avg_stem_mass+se_stem_mass), width=0.2) +
    facet_grid(.~Location) +
    labs(x="Crop Mix", y="Fresh Mass / Individual [g]", title = paste0(s)) +
    theme_bw(base_size = 24, base_family = "Helvetica")
  ggsave(paste0("output/Stem_mass_", s, ".png"))
}

## Land Equivalent Ratio
data_for_LER <- total_biomass %>%
  select(Location, CropTrt, CropSp, avg_LeavesStems_tha) %>%
  group_by(Location, CropTrt, CropSp) %>%
  summarize(site_avg_biomass = mean(avg_LeavesStems_tha),
            site_sd_biomass = sd(avg_LeavesStems_tha))

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