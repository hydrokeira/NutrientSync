require(dplyr)
require(ggplot2)

setwd("/Users/keirajohnson/Library/CloudStorage/Box-Box/Keira_Johnson/SiSyn/NutrientSynchrony")

chem<-read.csv("WRTDS_Outputs_Clean_02022026.csv")

stoich_class<-read.csv("Annual_Seasonality_Differences_Summary.csv")

Si_N_example<-stoich_class %>%
  mutate(solute_pair=paste(Solute1, "-", Solute2)) %>%
  filter(solute_pair=="Si - N")

streams<-c("Gide alv Gideabacka","WILLAMETTE RIVER", #stationary
           "ws3","NEUSE RIVER", #magnitude non-stationary
           "MERCED R", "Motala Strom Norrkoping", #timing non-stationary
           "I007.0W","OUSE AT SKELTON" #jointly non-stationary
           )

chem_norm<-chem %>%
  dplyr::filter(chemical %in% c("DSi", "N")) %>%
  dplyr::filter(Stream_Name %in% streams) %>%
  dplyr::group_by(chemical, Stream_Name) %>%
  mutate(norm_conc=scale(FNConc_mgL))

plot_list<-list()

#cols<-c("#66C5CC","#F6CF71","#F89C74","#DCB0F2","#87C55F","#9EB9F3","#FE88B1","#C9DB74")

cols<-c("#7c7c7c", "#7c7c7c", "#af6458", "#af6458","#d9af6b", "#d9af6b", "#855c75", "#855c75")

stream_names_clean<-c("Gide River (Sweden)", "Willamette River (OR, USA)",
                      "Hubbard Brook (NH, USA)", "Neuse River (NC, USA)",
                      "Merced River (CA, USA)", "Motala River (Sweden)",
                      "Illinois River (IL, USA)", "Ouse River (UK)")

for (i in seq_along(streams)) {
  
  one_stream<-chem_norm %>%
    filter(Stream_Name==streams[i])
  
  scale_factor=1
  
  plot_list[[i]]<-ggplot(one_stream,
         aes(Month, norm_conc,
             color = Year,
             group = interaction(Year, chemical),
             linetype = chemical)) +
    geom_line() +
    scale_color_gradient(low = "grey", high = "black") +
    scale_y_continuous(
      name = "",
      sec.axis = sec_axis(~ . / scale_factor)) +
    theme_classic()+theme(text = element_text(size = 15), legend.position = "null",
                          plot.title = element_text(color = cols[i]))+
    scale_x_continuous(breaks = seq(1,12,1), labels = seq(1,12,1))+
    ggtitle(paste(stream_names_clean[i]))
  
  
}

p2<-ggarrange(plotlist=plot_list, nrow = 4, ncol=2)

p2

notable_rivers<-stoich_class %>%
  filter(Stream_Name %in% streams) %>%
  mutate(solute_pair=paste(Solute1, "-", Solute2)) %>%
  filter(solute_pair=="Si - N") %>%
  mutate(Stream_Name=factor(Stream_Name, levels = streams))

stoich_class<-stoich_class %>%
  mutate(solute_pair=paste(Solute1, "-", Solute2)) %>%
  filter(solute_pair=="Si - N")

p1<-ggplot()+
  geom_point(notable_rivers, mapping = aes(sd_sync, s_tot, col=Stream_Name), size=4)+
  theme_classic()+
  geom_hline(yintercept = 0.8)+geom_vline(xintercept = 0.1)+
  facet_wrap(~solute_pair)+
  labs(x="SD of Shape Similarity", y="Offset Stability", col="River")+
  theme(text = element_text(size = 20), legend.position = "right")+
  #scale_color_manual(values = carto_pal(n=9, "Pastel"))+
  scale_color_manual(values = c("#7c7c7c", "#7c7c7c", "#af6458", "#af6458",
                                "#d9af6b", "#d9af6b", "#855c75", "#855c75"))

p1

pdf("River_Stationarity_Examples_V3.pdf", width = 8, height = 12)

ggarrange(p1, p2, heights = c(0.3, 0.7), nrow = 2)

dev.off()



