setwd("/Users/keirajohnson/Library/CloudStorage/Box-Box/Keira_Johnson/SiSyn/NutrientSynchrony")

chem<-read.csv("Annual_Seasonality_Differences.csv")

stoich_class<-read.csv("Annual_Seasonality_Differences_Summary.csv")

stoich_class <- stoich_class %>%
  filter(!Latitude==-1.9359)

chem_comb <- chem %>%
  mutate(Solute1=ifelse(Solute1=="DSi","Si",Solute1)) %>%
  left_join(stoich_class[2:17]) %>%
  mutate(solute_pair=paste(Solute1, "-", Solute2))

chem_comb %>%
  mutate(change_class=case_when(
    change_class=="stable"~"stationary",
    change_class=="timing stable"~"magnitude non-stationary",
    change_class=="magnitude stable"~"timing non-stationary",
    .default = "jointly non-stationary"
  )) %>%
  mutate(change_class=factor(change_class, levels = c("stationary", "magnitude non-stationary", 
                                                      "timing non-stationary", "jointly non-stationary"))) %>%
  mutate(solute_pair=factor(solute_pair, levels = c("Si - N", "Si - P", "N - P"))) %>%
  group_by(Stream_Name, solute_pair) %>%
  mutate(first_year = min(Year, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    Stream_Name = reorder_within(
      Stream_Name,
      -first_year,
      solute_pair
    )
  ) %>%
  ggplot(aes(Year, Stream_Name))+geom_tile(aes(fill=change_class))+facet_grid(rows=vars(solute_pair), cols=vars(change_class), scales="free")+
  theme_classic()+theme(axis.text.y = element_blank(), text=element_text(size = 20))+
  scale_fill_manual(values = c("#7c7c7c", "#af6458",
                                "#d9af6b", "#855c75"))+
  labs(x="Time", y="", fill="")
