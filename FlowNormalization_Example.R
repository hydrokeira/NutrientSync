setwd("/Users/keirajohnson/Library/CloudStorage/Box-Box/Keira_Johnson/SiSyn/NutrientSynchrony")

chem<-read.csv("WRTDS_Outputs_Clean_02022026.csv")

W<-chem %>%
  filter(Stream_Name=="WILLAMETTE RIVER")

p2<-W %>%
  ggplot(aes(Month, GenConc_mgL))+geom_line(aes(col=Year, group = Year))+
  facet_wrap(~chemical, scales="free")+theme_classic()+
  theme(text = element_text(size=20))+
  ggtitle("Willamette River - Not Flow Normalized")+
  scale_x_continuous(labels = seq(1,12,1), breaks = seq(1,12,1))+
  labs(y="Non-Flow Normalized\nConcentration (mg/L)", x="Month")

p1<-W %>%
  ggplot(aes(Month, FNConc_mgL))+geom_line(aes(col=Year, group = Year))+
  facet_wrap(~chemical, scales="free")+theme_classic()+
  theme(text = element_text(size=20))+
  ggtitle("Willamette River - Flow Normalized")+
  labs(y="Flow Normalized\nConcentration (mg/L)", x="")

k1<-ggarrange(p1, p2, nrow = 2)


O<-chem %>%
  filter(Stream_Name=="OUSE AT SKELTON")

p2<-O %>%
  ggplot(aes(Month, GenConc_mgL))+geom_line(aes(col=Year, group = Year))+
  facet_wrap(~chemical, scales="free")+theme_classic()+
  theme(text = element_text(size=20))+
  ggtitle("Ouse River - Not Flow Normalized")+
  scale_x_continuous(labels = seq(1,12,1), breaks = seq(1,12,1))+
  labs(y="Non-Flow Normalized\nConcentration (mg/L)", x="Month")

p1<-O %>%
  ggplot(aes(Month, FNConc_mgL))+geom_line(aes(col=Year, group = Year))+
  facet_wrap(~chemical, scales="free")+theme_classic()+
  theme(text = element_text(size=20))+
  ggtitle("Ouse River - Flow Normalized")+
  labs(y="Flow Normalized\nConcentration (mg/L)", x="")

k2<-ggarrange(p1, p2, nrow = 2)

pdf("FlowNormalization_Example_Rivers.pdf", width = 14, height = 16)

ggarrange(k1, k2, nrow = 2)

dev.off()
