####individual nutrient DTW####
require(dtw)
setwd("/Users/keirajohnson/Library/CloudStorage/Box-Box/Keira_Johnson/SiSyn/NutrientSynchrony")

chem<-read.csv("WRTDS_Outputs_Clean_02022026.csv")

stoich_class<-read.csv("Annual_Seasonality_Differences_Summary.csv")

chem_dist<-chem %>%
  filter(Stream_Name %in% stoich_class$Stream_Name) %>%
  group_by(Stream_Name, chemical, Year) %>%
  mutate(norm_conc=scale(FNConc_mgL), norm_Q=scale(Discharge_cms))

chem_dist_mean <- chem_dist %>%
  filter(chemical=="P") %>%
  group_by(chemical, Stream_Name, Month) %>%
  summarise(mean_conc=mean(norm_conc))

# chem_dist_mean <- chem_dist %>%
#   #filter(chemical=="DSi") %>%
#   group_by(chemical, Stream_Name, Month) %>%
#   summarise(mean_conc=mean(norm_conc), mean_Q=mean(norm_Q))

streams<-unique(chem_dist_mean$Stream_Name)

dist_mat_list<-list()

for(i in 1:length(streams)){
  
  print(i)
  
  chem_dist_test <- chem_dist %>%
    filter(chemical=="P") %>%
    filter(Stream_Name==streams[i]) %>%
    select(Stream_Name, Month, Year, norm_conc) %>%
    pivot_wider(names_from = Month, values_from = norm_conc, values_fn = mean)

  chem_dist_mean_test<-chem_dist_mean %>%
    filter(Stream_Name==streams[i]) %>%
    select(Stream_Name, Month, mean_conc) %>%
    pivot_wider(names_from = Month, values_from = mean_conc, values_fn = mean)
  # 
  # chem_dist_test <- chem_dist %>%
  #   filter(Stream_Name == streams[i]) %>%
  #   group_by(Stream_Name, Year, Month) %>%
  #   summarise(norm_Q = mean(norm_Q, na.rm = TRUE), .groups = "drop") %>%
  #   pivot_wider(names_from = Month, values_from = norm_Q)
  # 
  # chem_dist_mean_test<-chem_dist_mean %>%
  #   filter(Stream_Name==streams[i]) %>%
  #   select(Stream_Name, Month, mean_Q) %>%
  #   pivot_wider(names_from = Month, values_from = mean_Q, values_fn = mean)
  
  years<-unique(chem_dist_test$Year)
  
  norm_distance<-list()
  
  for (k in 1:length(years)) {
    
    aln <- dtw(
      as.matrix(chem_dist_test[k,4:15]), as.matrix(chem_dist_mean_test[1,3:14]), 
      window.type="sakoechiba", window.size=0L
    )
    
    norm_distance[[k]]<-aln$distance
    
  }
  
  dist_mat<-tibble(distance = unlist(norm_distance))
  
  dist_mat$Stream_Name<-chem_dist_test$Stream_Name
  dist_mat$Year<-chem_dist_test$Year
  
  dist_mat_list[[i]]<-dist_mat
  
  
}

dist_mat_df<-do.call(bind_rows, dist_mat_list)

dist_mat_df_sum<-dist_mat_df %>%
  group_by(Stream_Name) %>%
  summarise(mean_dist=mean(distance),
            sd_dist=sd(distance))

dist_mat_df_sum_Si<-dist_mat_df_sum
dist_mat_df_sum_N<-dist_mat_df_sum
dist_mat_df_sum_P<-dist_mat_df_sum
dist_mat_df_sum_Q<-dist_mat_df_sum

dist_mat_df_sum_Si$chemical<-"Si"
dist_mat_df_sum_N$chemical<-"N"
dist_mat_df_sum_P$chemical<-"P"
dist_mat_df_sum_Q$chemical<-"Q"
#dist_mat_df_sum_Q<-dist_mat_df_sum_Q[,c(1,2,3,5)]

dist_mat_all<-bind_rows(dist_mat_df_sum_Si, dist_mat_df_sum_N, dist_mat_df_sum_P, dist_mat_df_sum_Q)

write.csv(dist_mat_all, "Interannual_Distance_NonFN.csv")

ggplot(dist_mat_all, aes(chemical, mean_dist))+geom_boxplot(outliers = F)

stoich_class_dist<-stoich_class %>%
  left_join(dist_mat_all) %>%
  mutate(solute_pair=paste(Solute1, "-", Solute2)) %>%
  group_by(Stream_Name, solute_pair) %>%
  filter(chemical %in% c(Solute1, Solute2, "Q"))

pdf("DTW_Dist_Interannual_Nutrient.pdf", width = 15, height = 6)

g1<-stoich_class_dist %>%
  mutate(solute_pair=factor(solute_pair, levels = c("Si - N", "Si - P", "N - P"))) %>%
  mutate(chemical=factor(chemical, levels = c("Si","N", "P", "Q"))) %>%
  mutate(change_class=factor(change_class, levels=c("stable", "timing stable", "magnitude stable", "unstable"))) %>%
  ggplot(aes(x=change_class, y=mean_dist, fill=chemical))+
  geom_boxplot(outliers = F, alpha=0.8)+
  facet_wrap(~solute_pair, nrow=3)+
  labs(y="Mean Distance from Average Seasonal Pattern (unitless)", x="Stationarity Class", fill="")+
  theme_classic()+
  theme(text = element_text(size = 20))+
  scale_fill_manual(values = c("P"="black", "N"="grey40", "Si"="grey80", "Q"="lightblue"))

g1

dev.off()

#difference between chemicals across classes

m1 <- lmer(mean_dist ~ chemical * change_class + (1 | Stream_Name), data = filter(stoich_class_dist, chemical != "Q", solute_pair=="Si - N"))

anova(m1, type = 3)

emmeans(m1, pairwise ~ chemical | change_class, adjust = "tukey")

#differences across classes by chemical

# ANOVA with repeated observations accounted for by stream
m2 <- lmer(mean_dist ~ chemical * change_class + (1 | Stream_Name), data = filter(stoich_class_dist, chemical != "Q", solute_pair=="N - P"))

anova(m2, type = 3)

em2<-emmeans(m2, pairwise ~ change_class | chemical, adjust = "tukey")

cld(em2, adjust = "tukey", Letters = letters)

#test Q across stationarity classes

m3 <- aov(mean_dist ~ change_class, data = filter(stoich_class_dist, chemical == "Q", solute_pair == "N - P"))

summary(m3)

em3 <- emmeans(m3, ~ change_class)

cld(em3, adjust = "tukey", Letters = letters)



p1<-stoich_class_dist %>%
  mutate(solute_pair=factor(solute_pair, levels = c("Si - N", "Si - P", "N - P"))) %>%
  mutate(change_class=factor(change_class, levels=c("unstable", "magnitude stable", "timing stable", "stable"))) %>%
  filter(solute_pair=="Si - N") %>%
  dplyr::select(-sd_dist) %>%
  pivot_wider(names_from = chemical, values_from = mean_dist) %>%
  ggplot(aes(y=N, x=Si, col=change_class))+
  geom_point(size=2)+
  #facet_wrap(~solute_pair)+
  labs(x="Si Interannual Distance", y="N Interannual Distance", col="")+
  theme_classic()+
  theme(text = element_text(size = 20), legend.position = "null")+
  scale_color_manual(values = carto_pal(4, "Antique"))
#scale_fill_manual(values = c("P"="black", "N"="grey40", "Si"="grey80"))

p2<-stoich_class_dist %>%
  mutate(solute_pair=factor(solute_pair, levels = c("Si - N", "Si - P", "N - P"))) %>%
  mutate(change_class=factor(change_class, levels=c("unstable", "magnitude stable", "timing stable", "stable"))) %>%
  filter(solute_pair=="Si - P") %>%
  dplyr::select(-sd_dist) %>%
  pivot_wider(names_from = chemical, values_from = mean_dist) %>%
  ggplot(aes(y=P, x=Si, col=change_class))+
  geom_point(size=2)+
  #facet_wrap(~solute_pair)+
  labs(x="Si Interannual Distance", y="P Interannual Distance", col="")+
  theme_classic()+
  theme(text = element_text(size = 20), legend.position = "null")+
  scale_color_manual(values = carto_pal(4, "Antique"))
#scale_fill_manual(values = c("P"="black", "N"="grey40", "Si"="grey80"))

p3<-stoich_class_dist %>%
  mutate(solute_pair=factor(solute_pair, levels = c("Si - N", "Si - P", "N - P"))) %>%
  mutate(change_class=factor(change_class, levels=c("unstable", "magnitude stable", "timing stable", "stable"))) %>%
  filter(solute_pair=="N - P") %>%
  dplyr::select(-sd_dist) %>%
  pivot_wider(names_from = chemical, values_from = mean_dist) %>%
  ggplot(aes(y=P, x=N, col=change_class))+
  geom_point(size=2)+
  #facet_wrap(~solute_pair)+
  labs(x="N Interannual Distance", y="P Interannual Distance", col="")+
  theme_classic()+
  theme(text = element_text(size = 20))+
  scale_color_manual(values = carto_pal(4, "Antique"))+
  guides(color = guide_legend(reverse = TRUE))


pdf("Individual_DTW_Biplot.pdf", width = 15, height = 4.5)

l1<-ggarrange(p1, p2, p3, nrow = 1, widths = c(1,1,1.5))

dev.off()

l1

chem_wide_Q <- chem %>%
  filter(Stream_Name %in% stoich_class$Stream_Name) %>%
  group_by(Stream_Name, Month, Year) %>%
  mutate(
    Discharge_cms = mean(Discharge_cms, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  dplyr::select(
    Stream_Name, Month, Year,
    Discharge_cms, chemical, GenConc_mgL
  ) %>%
  pivot_wider(
    names_from = chemical,
    values_from = GenConc_mgL,
    values_fn = mean
  )
  
chem_norm_mean <- chem_wide_Q %>%
  pivot_longer(cols = c(4:7), values_to = "value", names_to = "chemical") %>%
  filter(!is.na(value)) %>%
  group_by(Stream_Name, chemical) %>%
  mutate(norm_conc=scale(value)) %>%
  ungroup() %>%
  group_by(Stream_Name, chemical, Year) %>%
  summarise(mean_conc = mean(norm_conc)) %>%
  summarise(sd_conc=sd(mean_conc)) %>%
  mutate(chemical=ifelse(chemical=="DSi", "Si", chemical)) %>%
  mutate(chemical=ifelse(chemical=="Discharge_cms", "Q", chemical))

stoich_class_conc_sd<-stoich_class %>%
  left_join(chem_norm_mean) %>%
  mutate(solute_pair=paste(Solute1, "-", Solute2)) %>%
  group_by(Stream_Name, solute_pair) %>%
  filter(chemical %in% c(Solute1, Solute2, "Q"))

g2<-stoich_class_conc_sd %>%
  mutate(solute_pair=factor(solute_pair, levels = c("Si - N", "Si - P", "N - P"))) %>%
  mutate(chemical=factor(chemical, levels = c("Si", "N", "P", "Q"))) %>%
  mutate(change_class=factor(change_class, levels=c("stable", "timing stable", "magnitude stable", "unstable"))) %>%
  #left_join(sig_sd) %>%
  ggplot(aes(x=change_class, y=sd_conc, fill=chemical))+
  geom_boxplot(outliers = F, alpha=0.8)+
  facet_wrap(~solute_pair, nrow=3)+
  labs(y="Standard Deviation(Mean Annual Z-Score Normalized Concentrations)", x="Stationarity Class", fill="")+
  theme_classic()+
  theme(text = element_text(size = 20))+
  scale_fill_manual(values = c("P"="black", "N"="grey40", "Si"="grey80", "Q"="lightblue"))

g2

#difference between chemicals across classes

m1 <- lmer(sd_conc ~ chemical * change_class + (1 | Stream_Name), data = filter(stoich_class_conc_sd, chemical != "Q", solute_pair=="Si - N"))

anova(m1, type = 3)

emmeans(m1, pairwise ~ chemical | change_class, adjust = "tukey")

#differences across classes by chemical

# ANOVA with repeated observations accounted for by stream
m2 <- lmer(sd_conc ~ chemical * change_class + (1 | Stream_Name), data = filter(stoich_class_conc_sd, chemical != "Q", solute_pair=="N - P"))

anova(m2, type = 3)

em2<-emmeans(m2, pairwise ~ change_class | chemical, adjust = "tukey")

cld(em2, adjust = "tukey", Letters = letters)

#test Q across stationarity classes

m3 <- aov(sd_conc ~ change_class, data = filter(stoich_class_conc_sd, chemical == "Q", solute_pair == "Si - N"))

summary(m3)

em3 <- emmeans(m3, ~ change_class)

cld(em3, adjust = "tukey", Letters = letters)

stoich_class_conc_sd %>%
  filter(chemical != "Q", solute_pair == "Si - N") %>%
  count(Stream_Name) %>%
  count(n, name = "number_of_streams")


k1<-stoich_class_conc_sd %>%
  mutate(solute_pair=factor(solute_pair, levels = c("Si - N", "Si - P", "N - P"))) %>%
  mutate(change_class=factor(change_class, levels=c("unstable", "magnitude stable", "timing stable", "stable"))) %>%
  filter(solute_pair=="Si - N") %>%
  #select(-sd_dist) %>%
  pivot_wider(names_from = chemical, values_from = sd_conc) %>%
  ggplot(aes(y=N, x=Si, col=change_class))+
  geom_point(size=2)+
  #facet_wrap(~solute_pair)+
  labs(x="Si SD Normalized Concentration", y="N SD Normalized Concentration", col="")+
  theme_classic()+
  theme(text = element_text(size = 20), legend.position = "null")+
  scale_color_manual(values = carto_pal(4, "Antique"))

k1

k2<-stoich_class_conc_sd %>%
  mutate(solute_pair=factor(solute_pair, levels = c("Si - N", "Si - P", "N - P"))) %>%
  mutate(change_class=factor(change_class, levels=c("unstable", "magnitude stable", "timing stable", "stable"))) %>%
  filter(solute_pair=="Si - P") %>%
  #select(-sd_dist) %>%
  pivot_wider(names_from = chemical, values_from = sd_conc) %>%
  ggplot(aes(y=P, x=Si, col=change_class))+
  geom_point(size=2)+
  #facet_wrap(~solute_pair)+
  labs(x="Si SD Normalized Concentration", y="P SD Normalized Concentration", col="")+
  theme_classic()+
  theme(text = element_text(size = 20), legend.position = "null")+
  scale_color_manual(values = carto_pal(4, "Antique"))

k2

k3<-stoich_class_conc_sd %>%
  mutate(solute_pair=factor(solute_pair, levels = c("Si - N", "Si - P", "N - P"))) %>%
  mutate(change_class=factor(change_class, levels=c("unstable", "magnitude stable", "timing stable", "stable"))) %>%
  filter(solute_pair=="N - P") %>%
  #select(-sd_dist) %>%
  pivot_wider(names_from = chemical, values_from = sd_conc) %>%
  ggplot(aes(y=P, x=N, col=change_class))+
  geom_point(size=2)+
  #facet_wrap(~solute_pair)+
  labs(x="N SD Normalized Concentration", y="P SD Normalized Concentration", col="")+
  theme_classic()+
  theme(text = element_text(size = 20), legend.position = "right")+
  scale_color_manual(values = carto_pal(4, "Antique"))+
  guides(color = guide_legend(reverse = TRUE))

k3

pdf("Individual_SD_Biplot.pdf", width = 18, height = 6)

l2<-ggarrange(k1, k2, k3, nrow = 1, widths = c(1,1,1.5))

dev.off()

pdf("DTW_Dist_SD_Interannual_Nutrient_nonFN_V2.pdf", width = 10, height = 10)

ggarrange(g1, g2, nrow = 1, widths = c(1, 1.1))

dev.off()

pdf("Individual_Nutrient_DTW_Conc_Biplot_Q_V2.pdf", width = 10, height = 18)

ggarrange(l1, l2, nrow = 2)

dev.off()

sig_sd<-stoich_class_conc_sd %>%
  group_by(change_class, solute_pair) %>%
  summarise(p.value = t.test(sd_conc ~ chemical)$p.value) %>%
  mutate(sig=ifelse(p.value < 0.05, "true", "false"))

sig_dtw<-stoich_class_dist %>%
  group_by(change_class, solute_pair) %>%
  summarise(p.value = t.test(mean_dist ~ chemical)$p.value) %>%
  mutate(sig=ifelse(p.value < 0.05, "true", "false"))
