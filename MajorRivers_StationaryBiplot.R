require(ggplot2)
require(trend)

stoich_class<-read.csv("Annual_Seasonality_Differences_Summary.csv")

ggplot(stoich_class, aes(se_sync, s_tot))+geom_point()+theme_classic()+
  geom_hline(yintercept = 0.8)+geom_vline(xintercept = 0.025)

spatial<-read.csv("Spatial_data_average_v3.csv")

spatial_kept<-spatial %>%
  dplyr::select(Research_Network, Stream_Name, watershed_slope, watershed_elevation_median, drainage_area,temperature, precipitation,
         evapotranspiration, npp, greenup, permafrost, snow_prop_area, Max_Daylength, major_rock,
         land_Forest, land_Impervious, land_Cropland, land_Wetland_Marsh)

pdf("Stationarity_Quadrant_V2.pdf", width = 12, height = 4)

stoich_class %>%
  left_join(spatial_kept) %>%
  mutate(solute_pair=paste(Solute1, "-", Solute2)) %>%
  mutate(solute_pair=factor(solute_pair, levels = c("Si - N", "Si - P", "N - P"))) %>%
  mutate(change_class=factor(change_class, levels = c("unstable", "magnitude stable", "timing stable", "stable"))) %>%
  ggplot(aes(CV_sync, s_tot, col=change_class))+geom_point(size=2)+theme_classic()+
  geom_hline(yintercept = 0.8)+geom_vline(xintercept = 0.2)+
  scale_color_manual(values = carto_pal(4,"Antique"))+facet_wrap(~solute_pair)+
  labs(x="SE in Shape Similarity", y="Modal Offset Stability")+
  theme(text = element_text(size = 20), legend.position = "null")

dev.off()

Seine_Rivers<-c("CARRIERES-SOUS-POISSY", "BENNECOURT", "AMFREVILLE-SOUS-LES-MONTS", "IVRY-SUR-SEINE",
                "PARIS-12E--ARRONDISSEMENT", "POSES 3", "TRIEL-SUR-SEINE")

big_rivers<-stoich_class %>%
  left_join(spatial_kept) %>%
  filter(Research_Network %in% c("UMR") | Stream_Name %in% Seine_Rivers) %>%
  mutate(LTER = case_when(
    #Research_Network=="MD"~"Murray Darling",
    Research_Network=="UMR"~"Upper Mississippi River",
    .default = "Seine"
  )) %>%
  mutate(solute_pair=paste(Solute1, "-", Solute2)) %>%
  mutate(solute_pair=factor(solute_pair, levels = c("Si - N", "Si - P", "N - P")))
  
stoich_class<-stoich_class %>%
  mutate(solute_pair=paste(Solute1, "-", Solute2)) %>%
  mutate(solute_pair=factor(solute_pair, levels = c("Si - N", "Si - P", "N - P")))

pdf("MajorRiver_NonStationarity.pdf", width = 13, height = 4)

p1<-ggplot()+geom_point(stoich_class, mapping=aes(CV_sync, s_tot), col="gray", size=3)+
  geom_point(big_rivers, mapping = aes(CV_sync, s_tot, col=LTER), size=3)+
  theme_classic()+
  geom_hline(yintercept = 0.8)+geom_vline(xintercept = 0.2)+
  facet_wrap(~solute_pair)+
  labs(x="SE in Shape Similarity", y="Modal Offset Stability", col="River")+
  theme(text = element_text(size = 20), legend.position = "bottom")

p1

dev.off()  

stoich_class %>%
  left_join(spatial_kept) %>%
  filter(Research_Network=="UK") %>%
  filter(!duplicated(Stream_Name)) %>%
  nrow()

river_names<-unique(big_rivers$Stream_Name)

setwd("/Users/keirajohnson/Library/CloudStorage/Box-Box/Keira_Johnson/SiSyn/NutrientSynchrony")

chem<-read.csv("WRTDS_Outputs_Clean_02022026.csv")

big_chem<-chem %>%
  filter(LTER=="UK")

chem_trends <- big_chem %>%
  group_by(LTER, Stream_Name, chemical) %>%
  arrange(Year, .by_group = TRUE) %>%
  summarise(
    
    # Mann-Kendall test
    mk = list(MannKendall(FNConc_mgL)),
    tau = mk[[1]]$tau,
    p_value = mk[[1]]$sl,
    
    # Sen's slope
    sen = list(sens.slope(FNConc_mgL)),
    sen_slope = sen[[1]]$estimates,
    
    # Direction based on Sen's slope and MK significance
    change_direction = case_when(
      p_value < 0.05 & sen_slope < 0 ~ "decreasing",
      p_value < 0.05 & sen_slope > 0 ~ "increasing",
      TRUE ~ "no significant trend"
    ),
    
    .groups = "drop"
  )

p2<-ggplot(chem_trends, aes(LTER, tau))+geom_boxplot(outliers = F)+facet_wrap(~chemical)+
  theme_classic()+geom_jitter(aes(shape = change_direction, fill = LTER))+geom_hline(yintercept = 0)+
  scale_shape_manual(values = c(25,24,1))+
  theme(text = element_text(size = 20), legend.position="null")+labs(x="", y="Mann Kendall tau", col="", fill="", 
                                             shape="Concentration Change")
p2

pdf("Figure_S8.pdf", width = 10, height = 10)

ggarrange(p3, p2, p1, nrow = 3, heights = c(0.5, 0.5, 0.7))

dev.off()

