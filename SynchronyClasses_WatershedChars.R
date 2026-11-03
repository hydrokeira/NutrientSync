setwd("/Users/keirajohnson/Library/CloudStorage/Box-Box/Keira_Johnson/SiSyn/NutrientSynchrony")

spatial<-read.csv("Spatial_data_average_v3.csv")

spatial_kept<-spatial %>%
  dplyr::select(Research_Network, Stream_Name, watershed_slope, watershed_elevation_median, drainage_area,temperature, precipitation,
         evapotranspiration, npp, greenup, permafrost, snow_prop_area, Max_Daylength, major_rock,
         land_Forest, land_Impervious, land_Cropland, land_Wetland_Marsh)

stoich_class<-read.csv("Annual_Seasonality_Differences_Summary.csv")

stoich_class <- stoich_class %>%
  mutate(solute_pair=paste(Solute1, "-", Solute2))

stoich_spatial <- stoich_class[,c(2,14, 17)] %>%
  left_join(spatial_kept) %>%
  mutate(change_class=factor(change_class, 
                             levels = c("stable", "timing stable", "magnitude stable", "unstable"))) %>%
  mutate(solute_pair=factor(solute_pair, levels = c("Si - N", "Si - P", "N - P")))

unique_sites<-stoich_spatial %>%
  filter(!duplicated(Stream_Name))

x1<-ggplot(unique_sites, aes(x=log(drainage_area)))+
  geom_histogram(fill="black")+
  theme_classic()+
  labs(x="Log(Drainage Area (km2))", y="Number of Sites")+
  theme(text = element_text(size = 22))+
  ylim(0,100)

x2<-ggplot(unique_sites, aes(x=land_Cropland/100))+
  geom_histogram(fill="black")+
  theme_classic()+
  labs(x="Proportional Cropland", y="")+
  theme(text = element_text(size = 22))+
  ylim(0,100)

x3<-ggplot(unique_sites, aes(x=land_Impervious/100))+
  geom_histogram(fill="black")+
  theme_classic()+
  labs(x="Proportional Urban", y="")+
  theme(text = element_text(size = 22))+
  ylim(0,100)


x3<-ggplot(unique_sites, aes(x=snow_prop_area))+
  geom_histogram(fill="black")+
  theme_classic()+
  labs(x="Proportional Snow Covered Area", y="")+
  theme(text = element_text(size = 22))+
  ylim(0,100)

pdf("WatershedChars_Distribution.pdf", width = 15, height = 4.5)

ggarrange(x1, x2, x3, nrow = 1)

dev.off()

p1<-ggplot(stoich_spatial, aes(change_class, log(drainage_area)))+
  geom_boxplot(aes(group = change_class), outliers = F)+
  facet_wrap(~solute_pair, nrow=3)+
  geom_jitter(alpha=0.2, size=0.5)+theme_classic()+
  labs(x="",y="Log(Drainage Area (km2))")+
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust=1))

pdf("Stationarity_Cropland_Prop.pdf", width=14, height = 5)

p2<-ggplot(stoich_spatial, aes(change_class, land_Cropland/100))+
  geom_boxplot(aes(group = change_class), outliers = F)+
  facet_wrap(~solute_pair, nrow=1)+
  geom_jitter(alpha=0.2, size=0.5)+theme_classic()+
  labs(x="", y="Proportional Cropland")+
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust=1))

p2

dev.off()

p3<-ggplot(stoich_spatial, aes(change_class, land_Impervious/100))+
  geom_boxplot(aes(group = change_class), outliers = F)+
  facet_wrap(~solute_pair, nrow=3)+
  geom_jitter(alpha=0.2, size=0.5)+theme_classic()+
  labs(x="", y="Proportional Cropland")+
  theme(text = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust=1))


p3<-ggplot(stoich_spatial, aes(change_class, snow_prop_area))+
  geom_boxplot(aes(group = change_class), outliers = F)+
  facet_wrap(~solute_pair, nrow=3)+
  geom_jitter(alpha=0.2, size=0.5)+theme_classic()+
  labs(x="", y="Proportional Snow Covered Area")+
  theme(text = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust=1))

p3


#calculate change in Q
discharge<-read.csv("WRTDS-input_discharge.csv")
discharge$Date<-as.Date(discharge$Date)

discharge<-discharge %>%
  mutate(Stream_Name=sub(".*__", "", Stream_ID)) %>%
  filter(Stream_Name %in% stoich_spatial$Stream_Name)

discharge_avg<-discharge %>%
  group_by(Stream_Name, year(Date)) %>%
  summarise(mean_Q=mean(Q))

Q_trend<-discharge_avg %>%
  group_by(Stream_Name) %>%
  arrange(`year(Date)`) %>%
  summarise(
    tau = MannKendall(mean_Q)$tau,
    p_value = MannKendall(mean_Q)$sl,
    .groups = "drop"
  )

Q_trend<-Q_trend %>%
  mutate(sig = ifelse(p_value <= 0.1, "yes", "no"))

stoich_spatial<-stoich_spatial %>%
  left_join(Q_trend)

p4<-stoich_spatial %>%
  filter(sig=="yes") %>%
  ggplot(aes(change_class, tau))+
  geom_boxplot(aes(group = change_class), outliers = F)+
  facet_wrap(~solute_pair, nrow=3)+
  geom_jitter()+theme_classic()+
  labs(x="")+theme(text = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust=1))+
  geom_hline(yintercept = 0)

p4

p5<-ggplot(stoich_spatial, aes(change_class, fill=sig))+
  geom_bar(position = "fill")+
  facet_wrap(~solute_pair, nrow=3)+
  theme_classic()+
  labs(x="", y="Proportion of Sites with Q Trend", fill="Sig")+theme(text = element_text(size = 20), 
                                                         legend.position = "right",
                   axis.text.x = element_text(angle = 45, hjust=1))+
  geom_hline(yintercept = 0)+
  scale_fill_manual(values = c("yes"="black","no"="grey"))+
  labs(fill="Significant\nTrend")

p5

pdf("WatershedChars_Stationarity.pdf", width = 12, height = 10)

ggarrange(p2, p1, p5, nrow = 1, widths = c(1,1,1.4))

p2

dev.off()
