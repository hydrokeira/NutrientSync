library(dplyr)
library(purrr)
library(tidyr)

setwd("/Users/keirajohnson/Library/CloudStorage/Box-Box/Keira_Johnson/SiSyn/NutrientSynchrony")

stoich_class<-read.csv("Annual_Seasonality_Differences_Summary.csv")

sd_thresh <- seq(0, 0.35, by = 0.005)
stot_thresh <- seq(0.3, 1, by = 0.01)

threshold_grid <- expand.grid(
  sd_sync_thresh = sd_thresh,
  s_tot_thresh = stot_thresh
)

calc_classes <- function(dat, sd_cut, stot_cut){
  
  dat %>%
    mutate(
      change_class = case_when(
        sd_sync < sd_cut & s_tot > stot_cut ~ "stable",
        sd_sync < sd_cut & s_tot <= stot_cut ~ "magnitude stable",
        sd_sync >= sd_cut & s_tot > stot_cut ~ "timing stable",
        TRUE ~ "unstable"
      )
    ) %>%
    count(change_class) %>%
    mutate(
      prop = n/sum(n),
      sd_sync_thresh = sd_cut,
      s_tot_thresh = stot_cut
    )
}

results <- stoich_class %>%
  group_by(Solute1, Solute2) %>%
  group_split() %>%
  map_dfr(function(df){
    
    map2_dfr(
      threshold_grid$sd_sync_thresh,
      threshold_grid$s_tot_thresh,
      ~calc_classes(df, .x, .y)
    ) %>%
      mutate(
        Solute1 = unique(df$Solute1),
        Solute2 = unique(df$Solute2)
      )
    
  })

p1<-results %>%
  filter(sd_sync_thresh==0.1) %>%
  mutate(solute_pair = paste(Solute1, Solute2, sep = " - ")) %>%
  mutate(solute_pair=factor(solute_pair, levels = c("Si - N", "Si - P", "N - P"))) %>%
  mutate(change_class_2=case_when(
    change_class=="unstable"~"jointly non-stationary",
    change_class=="magnitude stable"~"timing non-stationary",
    change_class=="timing stable"~"magnitude non-stationary",
    .default = "stationary"
  )) %>%
  mutate(change_class_2=factor(change_class_2, levels = c("jointly non-stationary", "timing non-stationary", 
                                                      "magnitude non-stationary", "stationary"))) %>%
  ggplot(aes(s_tot_thresh, prop))+
  geom_vline(xintercept = 0.8)+
  geom_line(aes(col=change_class_2), lwd=1)+
  theme_classic()+theme(text = element_text(size = 15))+
  scale_color_manual(values = carto_pal(4,"Antique"))+facet_wrap(~solute_pair)+
  labs(x="Modal Offset Stability Threshold" , y="Proportion of Sites", col="Stationary Class")+
  guides(color = guide_legend(reverse = TRUE))
  

p2<-results %>%
  filter(s_tot_thresh==0.8) %>%
  mutate(solute_pair = paste(Solute1, Solute2, sep = " - ")) %>%
  mutate(solute_pair=factor(solute_pair, levels = c("Si - N", "Si - P", "N - P"))) %>%
  mutate(change_class_2=case_when(
    change_class=="unstable"~"jointly non-stationary",
    change_class=="magnitude stable"~"timing non-stationary",
    change_class=="timing stable"~"magnitude non-stationary",
    .default = "stationary"
  )) %>%
  mutate(change_class_2=factor(change_class_2, levels = c("jointly non-stationary", "timing non-stationary", 
                                                      "magnitude non-stationary", "stationary"))) %>%
  ggplot(aes(sd_sync_thresh, prop))+
  geom_vline(xintercept = 0.1)+
  geom_line(aes(col=change_class_2), lwd=1)+
  theme_classic()+theme(text = element_text(size = 15))+
  scale_color_manual(values = carto_pal(4,"Antique"))+facet_wrap(~solute_pair)+
  labs(x="SE(Shape Similarity) Threshold", y="Proportion of Sites", col="Stationary Class")+
  guides(color = guide_legend(reverse = TRUE))

pdf("Threshold_Sensitivity_SD.pdf", width = 12, height = 6)

ggarrange(p1, p2, nrow = 2)

dev.off()

