setwd("/Users/keirajohnson/Library/CloudStorage/Box-Box/Keira_Johnson/SiSyn/NutrientSynchrony")

stoich_class<-read.csv("Annual_Seasonality_Differences_Summary.csv")

chem<-read.csv("WRTDS_Outputs_Clean_02022026.csv")

chem <- chem %>%
  filter(Stream_Name %in% stoich_class$Stream_Name)

chem_wide <- chem %>%
  dplyr::select(Stream_Name, Month, Year, chemical, GenConc_mgL) %>%
  pivot_wider(names_from = chemical, values_from = GenConc_mgL, values_fn = mean)

chem_wide<-chem_wide %>%
  mutate(Si_N = DSi/N, Si_P=DSi/P, N_P=N/P)

chem_wide_stoich <- chem %>%
  dplyr::select(Stream_Name, Month, Year, chemical, FNConc_uM) %>%
  pivot_wider(names_from = chemical, values_from = FNConc_uM, values_fn = mean)

chem_wide_stoich<-chem_wide_stoich %>%
  mutate(Si_N = DSi/N, Si_P=DSi/P, N_P=N/P)

stoich_class_Si_N <- stoich_class %>%
  filter(Solute1=="Si" & Solute2=="N")

p1<-chem_wide_stoich %>%
  left_join(stoich_class_Si_N[,c(2,14)]) %>%
  filter(!is.na(change_class)) %>%
  group_by(Stream_Name, change_class, Month) %>%
  summarise(mean_Si_N=mean(Si_N, na.rm = T)) %>%
  ungroup() %>%
  mutate(change_class=factor(change_class, levels = c("stable", "timing stable", "magnitude stable", "unstable"))) %>%
  ggplot(aes(change_class, log10(mean_Si_N)))+
  geom_boxplot(outliers = F, aes(fill=change_class))+
  scale_fill_manual(values = c("#7C7C7C", "#AF6458", "#D9AF6B", "#855C75"))+
  geom_hline(yintercept = log10(1), linetype=2)+
  theme_classic()+
  labs(x="", y="log(Mean Si:N)")+
  theme(text = element_text(size = 15), axis.text.x = element_blank(),
        legend.position = "null")
p1

Si_N_dat<-chem_wide_stoich %>%
  left_join(stoich_class_Si_N[,c(2,14)]) %>%
  filter(!is.na(change_class)) %>%
  group_by(Stream_Name, change_class) %>%
  summarise(mean_Si_N=log10(mean(Si_N, na.rm = T))) %>%
  ungroup() %>%
  dplyr::filter(mean_Si_N < mean(mean_Si_N)+4*sd(mean_Si_N)) 

mod <- aov(mean_Si_N ~ change_class, data = Si_N_dat)

cld_df_a <- emmeans(mod, ~ change_class) %>%
  multcomp::cld(Letters = letters, alpha=0.05)


stoich_class_Si_P <- stoich_class %>%
  filter(Solute1=="Si" & Solute2=="P")

p2<-chem_wide_stoich %>%
  left_join(stoich_class_Si_P[,c(2,14)]) %>%
  filter(!is.na(change_class)) %>%
  group_by(Stream_Name, change_class, Month) %>%
  summarise(mean_Si_P=mean(Si_P, na.rm=T)) %>%
  ungroup() %>%
  #filter(!Stream_Name %in% UK_sites$Stream_Name) %>%
  mutate(change_class=factor(change_class, levels = c("stable", "timing stable", "magnitude stable", "unstable"))) %>%
  ggplot(aes(change_class, log10(mean_Si_P)))+
  geom_boxplot(outliers = F, aes(fill = change_class))+
  #geom_jitter(alpha=0.5, aes(col=as.factor(Month)))+
  geom_hline(yintercept = log10(16), linetype=2)+
  theme_classic()+
  scale_fill_manual(values = c("#7C7C7C", "#AF6458", "#D9AF6B", "#855C75"))+
  labs(x="", y="log(Mean Si:P)", col="Month")+
  theme(text = element_text(size = 15), axis.text.x = element_blank(),
        legend.position = "null")

p2

Si_P_dat<-chem_wide_stoich %>%
  left_join(stoich_class_Si_P[,c(2,14)]) %>%
  filter(!is.na(change_class)) %>%
  group_by(Stream_Name, change_class) %>%
  summarise(mean_Si_P=log10(mean(Si_P, na.rm = T)))

mod <- aov(mean_Si_P ~ change_class, data = Si_P_dat)

cld_df_b <- emmeans(mod, ~ change_class) %>%
  multcomp::cld(Letters = letters)

stoich_class_N_P <- stoich_class %>%
  filter(Solute1=="N" & Solute2=="P")

p3<-chem_wide_stoich %>%
  left_join(stoich_class_N_P[,c(2,14)]) %>%
  filter(!is.na(change_class)) %>%
  group_by(Stream_Name, change_class, Month) %>%
  summarise(mean_N_P=mean(N_P, na.rm=T)) %>%
  ungroup() %>%
  #filter(!Stream_Name %in% UK_sites$Stream_Name) %>%
  mutate(change_class=factor(change_class, levels = c("stable", "timing stable", "magnitude stable", "unstable"))) %>%
  ggplot(aes(change_class, log10(mean_N_P)))+
  geom_boxplot(outliers = F, aes(fill = change_class))+
  geom_hline(yintercept = log10(16), linetype=2)+
  theme_classic()+
  scale_fill_manual(values = c("#7C7C7C", "#AF6458", "#D9AF6B", "#855C75"))+
  labs(x="", y="log(Mean N:P)")+
  theme(text = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "null")

p3

N_P_dat<-chem_wide_stoich %>%
  left_join(stoich_class_N_P[,c(2,14)]) %>%
  filter(!is.na(change_class)) %>%
  group_by(Stream_Name, change_class) %>%
  summarise(mean_N_P=log10(mean(N_P, na.rm = T)))

mod <- aov(mean_N_P ~ change_class, data = N_P_dat)

cld_df_c <- emmeans(mod, ~ change_class) %>%
  multcomp::cld(Letters = letters)

pdf("Monthly_Stoich_Class_Log_Molar.pdf", width = 3.5, height = 12)

k1<-ggarrange(p1, p2, p3, nrow = 3, heights = c(1,1,1.4), align = "v")

k1

dev.off()

k1

p1<-chem_wide_stoich %>%
  left_join(stoich_class_Si_N[,c(2,14)]) %>%
  filter(!is.na(Si_N), Si_N > 0) %>%
  group_by(Stream_Name, change_class, Month) %>%
  arrange(Year, .by_group = TRUE) %>%
  summarise(
    slope = if(n() >= 3) sens.slope(Si_N)$estimates else NA_real_,
    p = if(n() >= 3) mk.test(Si_N)$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  group_by(change_class) %>%
  mutate(
    mean_slope = mean(slope, na.rm = TRUE),
    sd_slope = sd(slope, na.rm = TRUE)
  ) %>%
  filter(abs(slope - mean_slope) <= 3 * sd_slope) %>%
  ungroup() %>%
  mutate(change_class=factor(change_class, levels = c("stable", "timing stable", "magnitude stable", "unstable"))) %>%
  ggplot(aes(change_class, slope)) +
  #geom_jitter(alpha = 0.4) +
  geom_boxplot(outliers = F, aes(fill=change_class)) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = c("#7C7C7C", "#AF6458", "#D9AF6B", "#855C75"))+
  theme_classic()+
  labs(x="", y="Sens Slope Si:N")+
  theme(text = element_text(size = 15), axis.text.x = element_blank(),
        legend.position = "null")

p1

Si_N_slope<-chem_wide_stoich %>%
  left_join(stoich_class_Si_N[,c(2,14)]) %>%
  filter(!is.na(Si_N), Si_N > 0) %>%
  group_by(Stream_Name, change_class, Month) %>%
  arrange(Year, .by_group = TRUE) %>%
  summarise(
    slope = if(n() >= 3) sens.slope(Si_N)$estimates else NA_real_,
    p = if(n() >= 3) mk.test(Si_N)$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  group_by(change_class) %>%
  mutate(
    mean_slope = mean(slope, na.rm = TRUE),
    sd_slope = sd(slope, na.rm = TRUE)
  ) %>%
  filter(abs(slope - mean_slope) <= 3 * sd_slope) %>%
  ungroup()
  #filter(p < 0.05)

mod <- aov(slope ~ change_class, data = Si_N_slope)

cld_df_d <- emmeans(mod, ~ change_class) %>%
  multcomp::cld(Letters = letters)

slope_sig_d<-Si_N_slope %>%
  dplyr::group_by(change_class) %>%
  summarise(
    p.value = t.test(slope, mu = 0)$p.value,
    .groups = "drop"
  )

p2<-chem_wide_stoich %>%
  left_join(stoich_class_Si_P[,c(2,14)]) %>%
  filter(!is.na(Si_P), Si_P > 0) %>%
  group_by(Stream_Name, change_class, Month) %>%
  arrange(Year, .by_group = TRUE) %>%
  summarise(
    slope = if(n() >= 3) sens.slope(Si_P)$estimates else NA_real_,
    p = if(n() >= 3) mk.test(Si_P)$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  group_by(change_class) %>%
  mutate(
    mean_slope = mean(slope, na.rm = TRUE),
    sd_slope = sd(slope, na.rm = TRUE)
  ) %>%
  filter(abs(slope - mean_slope) <= 3 * sd_slope) %>%
  #filter(p < 0.05) %>%
  ungroup() %>%
  mutate(change_class=factor(change_class, levels = c("stable", "timing stable", "magnitude stable", "unstable"))) %>%
  ggplot(aes(change_class, (slope))) +
  #geom_jitter(alpha = 0.8, aes(col=as.factor(Month))) +
  geom_boxplot(outliers = F, aes(fill=change_class)) +
  #geom_violin()+
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = c("#7C7C7C", "#AF6458", "#D9AF6B", "#855C75"))+
  #facet_wrap(~Month, scales="free_y") +
  #geom_jitter()+
  labs(x="", y="Sens Slope Si:P")+
  theme_classic()+
  #facet_wrap(~Month)+
  theme(text = element_text(size = 15), axis.text.x = element_blank(),
        legend.position = "null")
  #ylim(-1, 1)

p2

Si_P_slope<-chem_wide_stoich %>%
  left_join(stoich_class_Si_P[,c(2,14)]) %>%
  filter(!is.na(Si_P), Si_P > 0) %>%
  group_by(Stream_Name, change_class, Month) %>%
  arrange(Year, .by_group = TRUE) %>%
  summarise(
    slope = if(n() >= 3) sens.slope(Si_P)$estimates else NA_real_,
    p = if(n() >= 3) mk.test(Si_P)$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  group_by(change_class) %>%
  mutate(
    mean_slope = mean(slope, na.rm = TRUE),
    sd_slope = sd(slope, na.rm = TRUE)
  ) %>%
  filter(abs(slope - mean_slope) <= 3 * sd_slope) %>%
  ungroup()
  #filter(p < 0.05)

mod <- aov(slope ~ change_class, data = Si_P_slope)

cld_df_e <- emmeans(mod, ~ change_class) %>%
  multcomp::cld(Letters = letters)

slope_sig_e<-Si_P_slope %>%
  dplyr::group_by(change_class) %>%
  summarise(
    p.value = t.test(slope, mu = 0)$p.value,
    .groups = "drop"
  )

stoich_class_N_P <- stoich_class %>%
  filter(Solute1=="N" & Solute2=="P")

p3<-chem_wide_stoich %>%
  left_join(stoich_class_N_P[,c(2,14)]) %>%
  filter(!is.na(N_P), N_P > 0) %>%
  group_by(Stream_Name, change_class, Month) %>%
  arrange(Year, .by_group = TRUE) %>%
  summarise(
    slope = if(n() >= 3) sens.slope(N_P)$estimates else NA_real_,
    p = if(n() >= 3) mk.test(N_P)$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  #filter(p < 0.05) %>%
  #filter(!Stream_Name %in% UK_sites$Stream_Name) %>%
  group_by(change_class) %>%
  mutate(
    mean_slope = mean(slope, na.rm = TRUE),
    sd_slope = sd(slope, na.rm = TRUE)
  ) %>%
  filter(abs(slope - mean_slope) <= 3 * sd_slope) %>%
  ungroup() %>%
  mutate(change_class=factor(change_class, levels = c("stable", "timing stable", "magnitude stable", "unstable"))) %>%
  ggplot(aes(change_class, (slope))) +
  #geom_jitter(alpha = 0.8, aes(col=as.factor(Month))) +
  geom_boxplot(outliers = F, aes(fill=change_class)) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = c("#7C7C7C", "#AF6458", "#D9AF6B", "#855C75"))+
  theme_classic()+
  labs(x="", y="Sen's Slope N:P")+
  theme(text = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "null")
  #ylim(-0.2, 1)

p3

N_P_slope<-chem_wide_stoich %>%
  left_join(stoich_class_N_P[,c(2,14)]) %>%
  filter(!is.na(N_P), N_P > 0) %>%
  group_by(Stream_Name, change_class, Month) %>%
  arrange(Year, .by_group = TRUE) %>%
  summarise(
    slope = if(n() >= 3) sens.slope(N_P)$estimates else NA_real_,
    p = if(n() >= 3) mk.test(N_P)$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  group_by(change_class) %>%
  mutate(
    mean_slope = mean(slope, na.rm = TRUE),
    sd_slope = sd(slope, na.rm = TRUE)
  ) %>%
  filter(abs(slope - mean_slope) <= 3 * sd_slope) %>%
  ungroup()
#filter(p < 0.05)

q3<-N_P_slope %>%
  mutate(sig = ifelse(p <=0.05, "yes", "no")) %>%
  mutate(change_class=factor(change_class, levels = c("stable", "timing stable", "magnitude stable", "unstable"))) %>%
  ggplot(aes(change_class, fill=sig))+geom_bar(position = "fill", col="black")+theme_classic()+
  theme(text = element_text(size = 20), legend.position = "null")+labs(x="", y="Proportion of Sites")+
  scale_fill_manual(values = c("yes"="black", no="white"))+ggtitle("N-P")

q2<-Si_P_slope %>%
  mutate(sig = ifelse(p <=0.05, "yes", "no")) %>%
  mutate(change_class=factor(change_class, levels = c("stable", "timing stable", "magnitude stable", "unstable"))) %>%
  ggplot(aes(change_class, fill=sig))+geom_bar(position = "fill", col="black")+theme_classic()+
  theme(text = element_text(size = 20), legend.position = "null")+labs(x="", y="Proportion of Sites")+
  scale_fill_manual(values = c("yes"="black", no="white"))+ggtitle("Si-P")

q1<-Si_N_slope %>%
  mutate(sig = ifelse(p <=0.05, "yes", "no")) %>%
  mutate(change_class=factor(change_class, levels = c("stable", "timing stable", "magnitude stable", "unstable"))) %>%
  ggplot(aes(change_class, fill=sig))+geom_bar(position = "fill", col="black")+theme_classic()+
  theme(text = element_text(size = 20), legend.position = "null")+labs(x="", y="Proportion of Sites")+
  scale_fill_manual(values = c("yes"="black", no="white"))+ggtitle("Si-N")

ggarrange(q1, q2, q3, nrow = 1)

mod <- aov(slope ~ change_class, data = N_P_slope)

cld_df_f <- emmeans(mod, ~ change_class) %>%
  multcomp::cld(Letters = letters)

slope_sig_f<-N_P_slope %>%
  dplyr::group_by(change_class) %>%
  summarise(
    p.value = t.test(slope, mu = 0)$p.value,
    .groups = "drop"
  )

k2<-ggarrange(p1, p2, p3, nrow = 1, align = "h")

k2

ggarrange(k1, k2, nrow = 1)

pdf("Monthly_Stoich_Slope_Class.pdf", width = 5, height = 12)

k2<-ggarrange(p1, p2, p3, nrow = 3, heights = c(1,1,1.4), align = "v")

dev.off()

pdf("Monthly_Stoich_Mean_Slope_Class.pdf", width = 7, height = 11)

ggarrange(k1, k2, align = "h", widths = c(1, 1.04))


streams<-c("Gide alv Gideabacka","WILLAMETTE RIVER", #stationary
           "ws3","NEUSE RIVER", #magnitude non-stationary
           "MERCED R", "Motala Strom Norrkoping", #timing non-stationary
           "I007.0W","OUSE AT SKELTON" #jointly non-stationary
)

stream_names_clean<-c("Gide River (Sweden)", "Willamette River (OR, USA)",
                      "Hubbard Brook (NH, USA)", "Neuse River (NC, USA)",
                      "Merced River (CA, USA)", "Motala River (Sweden)",
                      "Illinois River (IL, USA)", "Ouse River (UK)")

stream_names_clean<-c("Gide alv Gideabacka" = "Gide River (Sweden)","WILLAMETTE RIVER"="Willamette River (OR, USA)",
                      "ws3"="Hubbard Brook (NH, USA)","NEUSE RIVER"="Neuse River (NC, USA)",
                      "MERCED R"="Merced River (CA, USA)","Motala Strom Norrkoping"="Motala River (Sweden)",
                      "I007.0W"="Illinois River (IL, USA)", "OUSE AT SKELTON"="Ouse River (UK)")

# chem_wide_stoich %>%
#   left_join(stoich_class_Si_N[,c(2,13)]) %>%
#   filter(!is.na(Si_N), Si_N > 0) %>%
#   filter(Stream_Name %in% streams) %>%
#   mutate(Stream_Name=factor(Stream_Name, levels = streams)) %>%
#   ggplot(aes(Year, log(Si_N)))+geom_line(aes(colour = Month, group = Month), size=1)+
#   facet_wrap(~Stream_Name, axes="all_y", scales="free", nrow=4, ncol=2, labeller=labeller(Stream_Name=stream_names_clean))+
#   theme_classic()+
#   theme(text = element_text(size = 15))+
#   scale_color_gradientn(colors=c("lightblue", "forestgreen","goldenrod","lightblue"))+
#   labs(x="Time", y="log10(Si:N)")

chem_Q <- chem %>%
  dplyr::select(Stream_Name, Month, Year, Discharge_cms) %>%
  filter(Stream_Name %in% streams) %>%
  filter(!duplicated(Discharge_cms))

chem_wide_stoich %>%
  left_join(stoich_class_Si_N[,c(2,13)]) %>%
  filter(!is.na(Si_N), Si_N > 0) %>%
  filter(Stream_Name %in% streams) %>%
  mutate(Stream_Name=factor(Stream_Name, levels = streams)) %>%
  ggplot(aes(Month, log(Si_N)))+geom_line(aes(colour = Year, group = Year), size=1, alpha=0.7)+
  facet_wrap(~Stream_Name, axes="all_y", scales="free", nrow=4, ncol=2, labeller=labeller(Stream_Name=stream_names_clean))+
  theme_classic()+
  theme(text = element_text(size = 15))+
  scale_color_viridis_c(option = "F")
  labs(x="Time", y="log(Si:N)")
  
cols<-c("#7c7c7c", "#7c7c7c", "#af6458", "#af6458","#d9af6b", "#d9af6b", "#855c75", "#855c75")
plot_list_1<-list()
plot_list_2<-list()


for (i in seq_along(streams)) {
    
  one_stream <- chem_wide_stoich %>%
    left_join(stoich_class_Si_N[, c(2, 13)]) %>%
    left_join(chem_Q) %>%
    filter(!is.na(Si_N), Si_N > 0) %>%
    filter(Stream_Name %in% streams[[i]])
    
    # plot_list_1[[i]] <- ggplot(one_stream, aes(Month)) +
    #   geom_line(aes(y = log(Si_N), colour = Year, group = Year),linewidth = 0.5,alpha = 0.7) +
    #   scale_color_viridis_c(option = "F") +
    #   ggnewscale::new_scale_color() +
    #   geom_line(aes(y = scale_factor * log(Discharge_cms) + scale_offset,colour = Year,group = Year),
    #     linewidth = 0.5,alpha = 0.6) +
    #   scale_color_viridis_c(option = "G") +
    #   scale_y_continuous(
    #     name = "log(Si:N)",
    #     sec.axis = sec_axis(
    #       ~ (. - scale_offset) / scale_factor,
    #       name = "Normalized discharge")) +
    #   scale_x_continuous(breaks = 1:12, labels = 1:12) +
    #   theme_classic() +theme(text = element_text(size = 15), legend.position = "none",
    #                          plot.title = element_text(color = cols[[i]])) +
    #   ggtitle(stream_names_clean[i]) +
    #   labs(x = "Time")
    # 
    plot_list_1[[i]] <- ggplot(one_stream, aes(Month)) +
      geom_line(aes(y = log10(Discharge_cms), colour = Year, group = Year),linewidth = 0.5,alpha = 0.7) +
      scale_color_gradient(high = cols[[i]], low = "grey")+
      scale_x_continuous(breaks = 1:12, labels = 1:12) +
      theme_classic() +theme(text = element_text(size = 15), legend.position = "none",
                             plot.title = element_text(color = cols[[i]])) +
      ggtitle(stream_names_clean[i]) +
      labs(x = "Month", y="log(Q)")

    plot_list_2[[i]] <- ggplot(one_stream, aes(Month)) +
      geom_line(aes(y = log10(Si_N), colour = Year, group = Year),linewidth = 0.5,alpha = 0.7) +
      scale_color_gradient(high = cols[[i]], low = "grey")+
      scale_x_continuous(breaks = 1:12, labels = 1:12) +
      theme_classic() +theme(text = element_text(size = 15), legend.position = "none",
                             plot.title = element_text(color = cols[[i]])) +
      ggtitle(stream_names_clean[i]) +
      labs(x = "Month", y="log(Si:N)")

    
  }
  
p2<-ggarrange(plotlist=plot_list_1, nrow = 4, ncol=2)

p2
  
p3<-ggarrange(plotlist=plot_list_2, nrow = 4, ncol=2, align = "hv")

p3

for (i in seq_along(streams)) {
  
  plot_list_2[[i]] <- chem_wide_stoich %>%
    left_join(stoich_class_Si_N[, c(2, 13)]) %>%
    left_join(chem_Q) %>%
    filter(!is.na(Si_N), Si_N > 0) %>%
    filter(Stream_Name %in% streams[[i]]) %>%
    ggplot(aes(x=log(Discharge_cms), y=log10(Si_N)))+
    geom_smooth(aes(col=Year, group = Year), method = "lm", se=F, size=0.5)+
    geom_point(aes(col=Year), alpha=0.1)+
    scale_color_gradient(high = cols[[i]], low = "grey")+
    theme_classic()+
    theme(legend.position = "null", text = element_text(size = 15))+
    labs(y="", x="log(Q)")
   
  
}

p4<-ggarrange(plotlist=plot_list_2, nrow = 4, ncol=2, align = "hv")

p4

h1<-ggarrange(p3, p4, widths = c(1, 0.6), nrow = 1)

pdf("Stoich_Change_Examples_Rivers.pdf", width = 11.5, height = 10)

p3

dev.off()  


