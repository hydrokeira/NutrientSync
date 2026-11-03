setwd("/Users/keirajohnson/Library/CloudStorage/Box-Box/Keira_Johnson/SiSyn/NutrientSynchrony")

stoich<-read.csv("Annual_Seasonality_Differences.csv")

pdf("POR_Nutrients.pdf", width = 12, height = 7)

stoich %>%
  mutate(solute_pair = paste(Solute1, "-", Solute2)) %>%
  group_by(Stream_Name, solute_pair) %>%
  mutate(n_years = n_distinct(Year), min_year=min(Year)) %>%
  ungroup() %>%
  mutate(Stream_Name = reorder_within(Stream_Name, -min_year, solute_pair)) %>%
  ggplot(aes(Year, Stream_Name, fill=n_years))+geom_tile()+facet_wrap(~solute_pair, scales="free_y")+
  theme_classic()+theme(text = element_text(size = 15), axis.text.y = element_blank())

dev.off()  

