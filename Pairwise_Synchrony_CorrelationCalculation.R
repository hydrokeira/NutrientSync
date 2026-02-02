require(dplyr)
require(ggplot2)
require(tidyr)
require(ggpubr)

#### define functions ####

compute_corr_synchrony <- function(P, Q) {
  # Step 1: compute mean monthly values (length 12 each)
  rm <- tapply(P, rep(1:12, length.out = length(P)), mean, na.rm = TRUE)
  km <- tapply(Q, rep(1:12, length.out = length(Q)), mean, na.rm = TRUE)
  
  # Step 2: ensure both are numeric vectors of same length
  rm <- as.numeric(rm)
  km <- as.numeric(km)
  
  # Step 3: circularly shift km and compute correlation at each shift
  shifts <- 0:11
  cors <- sapply(shifts, function(w) {
    q_shift <- km[((seq_along(km) - w - 1) %% length(km)) + 1]
    cor(rm, q_shift, use = "pairwise.complete.obs")
  })
  
  base_corr<-abs(cors[1])
  
  # Step 4: find best alignment and return results
  best_shift <- shifts[which.max(cors)]
  best_cor <- max(cors, na.rm=TRUE)
  
  # amplitude ratio
  sd1<-abs(sd(rm, na.rm=TRUE))
  sd2<-abs(sd(km, na.rm=TRUE))
  
  sd_max<-max(sd1, sd2)
  sd_min<-min(sd1, sd2)
  
  amp_ratio <- sd_min/sd_max
  
  synchrony <- best_cor * amp_ratio
  
  list(synchrony = synchrony,
       base_cor = base_corr,
       shape_cor = best_cor,
       amp_ratio = amp_ratio,
       best_shift = best_shift)
}

shift_vector <- function(x, shift) {
  n <- length(x)
  shift <- ((shift %% n) + n) %% n
  x[(seq_len(n) - shift - 1) %% n + 1]
}

#### read in, clean data ####

setwd("/Users/keirajohnson/Library/CloudStorage/Box-Box/Keira_Johnson/SiSyn/NutrientSynchrony")

chem<-read.csv("WRTDS_Outputs_Clean_02022026.csv")

chem_norm<-chem %>%
  group_by(Stream_Name, chemical) %>%
  mutate(norm_conc=scale(FNConc_mgL))

chem_avg<-chem_norm %>%
  group_by(Stream_Name, chemical, Month) %>%
  summarise(mean_conc=mean(norm_conc, na.rm=T))

chem_wide<-chem_avg %>%
  pivot_wider(names_from = Month, values_from = mean_conc)

cluster_si<-chem_wide %>%
  filter(chemical=="DSi")

cluster_n<-chem_wide %>%
  filter(chemical=="N")

cluster_p<-chem_wide %>%
  filter(chemical=="P")

#### for Si - N pairs ####

results_Si_N <- data.frame(Solute1 = character(0), Solute2 = character(0), corr = numeric(0))

Si_N_data <- chem_avg %>%
  filter(chemical %in% c("DSi", "N")) %>%
  group_by(Stream_Name) %>%
  mutate(num_chem=n_distinct(chemical)) %>%
  filter(num_chem==2)

Si_N_sites<-unique(Si_N_sites$Stream_Name)

#pdf("ASI_DSi_N_plots_CORRsynchrony.pdf", width = 10, height = 6)

# Compute ASI for each pair of lines
for (i in 1:length(Si_N_sites)) {
  print(i)
  line1 <- cluster_si %>%
    filter(Stream_Name==Si_N_sites[i]) %>%
    ungroup() %>%
    select(,3:14) %>%
    as.matrix()
  
  line2 <- cluster_n %>%
    filter(Stream_Name==Si_N_sites[i]) %>%
    ungroup() %>%
    select(,3:14) %>%
    as.matrix()
  
  # Compute ASI for the pair
  #asi_value <- compute_asi(line1, line2)
  
  corr_value<-compute_corr_synchrony(line1, line2)
  
  # apply the shift
  #shifted_line2 <- shift_vector(line2, 12-asi_value$best_shift)
  
  # shifted_line2 <- shift_vector(line2, corr_value$best_shift)
  # 
  # p1<-ggplot()+geom_line(aes(x=seq(1,12,1), y=line1))+geom_line(aes(x=seq(1,12,1), y=line2), col="red")+
  #   geom_line(aes(x = seq(1, 12, 1), y = shifted_line2), color = "red", linetype = "dashed")+
  #   labs(title=sites[i], subtitle=paste("Synchrony=", round(corr_value$synchrony, 2),
  #                                       "\nPhase Offset=",corr_value$best_shift),
  #        x="Month", y="Normalized Concentration")+theme_classic()+
  #   theme(text = element_text(size = 20))+
  #   scale_x_continuous(breaks = seq(1,12,1))
  # 
  # print(p1)
  
  # Store the result
  results_Si_N <- rbind(results_Si_N, 
                        data.frame(Stream_Name = Si_N_sites[i],
                                   Solute1 = "DSi", Solute2 = "N", corr = corr_value))
}

#dev.off()

results_Si_N$shift<-ifelse(results_Si_N$corr.best_shift <=6, results_Si_N$corr.best_shift, abs(12-results_Si_N$corr.best_shift))

# Create a regular 0–1 grid
grid <- expand.grid(
  corr.shape_cor = seq(0, 1, length.out = 200),
  corr.amp_ratio = seq(0, 1, length.out = 200)
)

line_df <- data.frame(
  corr.shape_cor = seq(0.001, 1, length.out = 200),
  corr.amp_ratio = 0.5 / seq(0.001, 1, length.out = 200)
)

line_df <- subset(line_df, corr.amp_ratio <= 1)

# Define what the background gradient should represent:
# For example, a simple function of both axes:
grid$corr.synchrony <- (grid$corr.shape_cor*grid$corr.amp_ratio)
# (you can replace this with another function, e.g., corr.shape_cor * corr.amp_ratio)

z1<-ggplot() +
  geom_raster(data = grid, mapping=aes(corr.shape_cor, corr.amp_ratio, fill = corr.synchrony)) +
  geom_point(results_Si_N, mapping=aes(corr.shape_cor, corr.amp_ratio), 
             color = "black", size=0.5) +
  geom_line(data = line_df, aes(corr.shape_cor, corr.amp_ratio),
            color = "black") +
  scale_fill_gradientn(
    colours = c("grey", "grey", "blue"),  # from blue to grey
    values = c(0, 0.47, 1),        # everything above 0.5 stays grey
    limits = c(0, 1),
    name = "Synchrony"
  ) +
  coord_fixed()+
  theme_classic()+
  labs(x="Correlation Coefficient", y="Amplitude Ratio")+
  ggtitle("DSi-N")+
  theme(text = element_text(size = 20), legend.position = "null")

z1

p1<-ggplot(results_Si_N, aes(shift, corr.synchrony))+
  geom_hline(yintercept = mean(results_Si_N$corr.synchrony), col="red")+
  #geom_hline(yintercept = 0.85, col="grey")+
  #geom_vline(xintercept = 0.5, col="grey")+geom_vline(xintercept = 5.5, col="grey")+
  geom_boxplot(aes(group=shift), fill=NA)+geom_jitter()+
  theme_classic()+theme(text = element_text(size = 20))+
  labs(x="Shift (Months)", y="Synchrony")+ggtitle("Si-N")

p1

p4<-ggplot(results_Si_N, aes(shift))+
  geom_bar(fill="blue")+
  theme_classic()+theme(text = element_text(size = 20))+
  labs(x="Shift (Months)", y="Number of Sites")
p4

results_Si_N<-results_Si_N %>%
  mutate(shift_class = case_when(
    shift==0 ~ "in phase",
    shift > 5 ~ "out of phase",
    .default = "offset"
  ),
  sync_class = case_when(
    corr.synchrony >=0.8~"synchronous",
    .default = "asynchronous"
  ),
  sync_shift=paste0(sync_class, ", ", shift_class)
  )

prop_table_Si_N <- results_Si_N %>%
  count(sync_class, shift_class) %>%
  mutate(prop = n / sum(n))

p4<-ggplot(prop_table_Si_N, aes(y = sync_class, x = shift_class, fill = prop)) +
  geom_tile(color = "white", alpha=0.8) +
  geom_text(aes(label = scales::percent(prop, accuracy = .1)), color = "black", size = 5) +
  scale_fill_viridis_c(option = "plasma", name = "Proportion", limits=c(0,0.5)) +
  theme_minimal(base_size = 18) +
  theme(legend.position = "null")+
  labs(x = "", y = "")

p4


#### now for Si - P pairs ####

results_Si_P <- data.frame(Solute1 = character(0), Solute2 = character(0), corr = numeric(0))

Si_P_data <- chem_avg %>%
  filter(chemical %in% c("DSi", "P")) %>%
  group_by(Stream_Name) %>%
  mutate(num_chem=n_distinct(chemical)) %>%
  filter(num_chem==2)

Si_P_sites<-unique(Si_P_data$Stream_Name)

# Compute ASI for each pair of lines
for (i in 1:length(Si_P_sites)) {
  print(i)
  
  line1 <- cluster_si %>%
    filter(Stream_Name==Si_P_sites[i]) %>%
    ungroup() %>%
    select(,3:14) %>%
    as.matrix()
  
  line2 <- cluster_p %>%
    filter(Stream_Name==Si_P_sites[i]) %>%
    ungroup() %>%
    select(,3:14) %>%
    as.matrix()
  
  corr_value<-compute_corr_synchrony(line1, line2)
  
  # apply the shift
  #shifted_line2 <- shift_vector(line2, 12-asi_value$best_shift)
  
  # shifted_line2 <- shift_vector(line2, corr_value$best_shift)
  # 
  # p1<-ggplot()+geom_line(aes(x=seq(1,12,1), y=line1))+geom_line(aes(x=seq(1,12,1), y=line2), col="red")+
  #   geom_line(aes(x = seq(1, 12, 1), y = shifted_line2), color = "red", linetype = "dashed")+
  #   labs(title=sites[i], subtitle=paste("Synchrony=", round(corr_value$synchrony, 2),
  #                                       "\nPhase Offset=",corr_value$best_shift),
  #        x="Month", y="Normalized Concentration")+theme_classic()+
  #   theme(text = element_text(size = 20))+
  #   scale_x_continuous(breaks = seq(1,12,1))
  # 
  # print(p1)
  
  # Store the result
  results_Si_P <- rbind(results_Si_P, 
                        data.frame(Stream_Name = Si_P_sites[i],
                                   Solute1 = "DSi", Solute2 = "P", corr = corr_value))
}

results_Si_P$shift<-ifelse(results_Si_P$corr.best_shift <=6, results_Si_P$corr.best_shift, abs(12-results_Si_P$corr.best_shift))

# Create a regular 0–1 grid
grid <- expand.grid(
  corr.shape_cor = seq(0, 1, length.out = 200),
  corr.amp_ratio = seq(0, 1, length.out = 200)
)

line_df <- data.frame(
  corr.shape_cor = seq(0.001, 1, length.out = 200),
  corr.amp_ratio = 0.5 / seq(0.001, 1, length.out = 200)
)

line_df <- subset(line_df, corr.amp_ratio <= 1)

# Define what the background gradient should represent:
# For example, a simple function of both axes:
grid$corr.synchrony <- (grid$corr.shape_cor*grid$corr.amp_ratio)
# (you can replace this with another function, e.g., corr.shape_cor * corr.amp_ratio)

z2<-ggplot() +
  geom_raster(data = grid, mapping=aes(corr.shape_cor, corr.amp_ratio, fill = corr.synchrony)) +
  geom_point(results_Si_P, mapping=aes(corr.shape_cor, corr.amp_ratio), 
             color = "black", size=0.5) +
  geom_line(data = line_df, aes(corr.shape_cor, corr.amp_ratio),
            color = "black") +
  scale_fill_gradientn(
    colours = c("grey", "grey", "blue"),  # from blue to grey
    values = c(0, 0.47, 1),        # everything above 0.5 stays grey
    limits = c(0, 1),
    name = "Synchrony"
  ) +
  coord_fixed()+
  theme_classic()+
  labs(x="Correlation Coefficient", y="")+
  ggtitle("DSi-P")+
  theme(text = element_text(size = 20), legend.position = "null")

z2

p2<-ggplot(results_Si_P, aes(shift, corr.synchrony))+
  geom_hline(yintercept = mean(results_Si_P$corr.synchrony), col="red")+
  #geom_hline(yintercept = 0.85, col="grey")+
  #geom_vline(xintercept = .5, col="grey")+geom_vline(xintercept = 5.5, col="grey")+
  geom_boxplot(aes(group=shift), fill=NA)+geom_jitter()+
  theme_classic()+theme(text = element_text(size = 20))+
  labs(x="Shift (Months)", y="")+ggtitle("Si-P")

p2

p5<-ggplot(results_Si_P, aes(shift))+
  geom_bar(fill="blue")+
  theme_classic()+theme(text = element_text(size = 20))+
  labs(x="Shift (Months)", y="Number of Sites")
p5

results_Si_P<-results_Si_P %>%
  mutate(shift_class = case_when(
    shift==0 ~ "in phase",
    shift > 5 ~ "out of phase",
    .default = "offset"
  ),
  sync_class = case_when(
    corr.synchrony >=0.85~"synchronous",
    .default = "asynchronous"
  ),
  sync_shift=paste0(sync_class, ", ", shift_class)
  )

prop_table_Si_P <- results_Si_P %>%
  count(sync_class, shift_class) %>%
  mutate(prop = n / sum(n))

p5<-ggplot(prop_table_Si_P, aes(y = sync_class, x = shift_class, fill = prop)) +
  geom_tile(color = "white", alpha=0.8) +
  geom_text(aes(label = scales::percent(prop, accuracy = .1)), color = "black", size = 5) +
  scale_fill_viridis_c(option = "plasma", name = "Proportion", limits=c(0,0.5)) +
  theme_minimal(base_size = 18) +
  theme(legend.position = "null", axis.text.y = element_blank())+
  labs(x = "", y = "")

p5

#### now for N - P pairs ####

results_N_P <- data.frame(Solute1 = character(0), Solute2 = character(0), corr = numeric(0))

N_P_data <- chem_avg %>%
  filter(chemical %in% c("N", "P")) %>%
  group_by(Stream_Name) %>%
  mutate(num_chem=n_distinct(chemical)) %>%
  filter(num_chem==2)

N_P_sites<-unique(N_P_data$Stream_Name)

# Compute ASI for each pair of lines
for (i in 1:length(N_P_sites)) {
  print(i)
  
  line1 <- cluster_n %>%
    filter(Stream_Name==N_P_sites[i]) %>%
    ungroup() %>%
    select(,3:14) %>%
    as.matrix()
  
  line2 <- cluster_p %>%
    filter(Stream_Name==N_P_sites[i]) %>%
    ungroup() %>%
    select(,3:14) %>%
    as.matrix()
  
  corr_value<-compute_corr_synchrony(line1, line2)
  
  # apply the shift
  #shifted_line2 <- shift_vector(line2, 12-asi_value$best_shift)
  
  # shifted_line2 <- shift_vector(line2, corr_value$best_shift)
  # 
  # p1<-ggplot()+geom_line(aes(x=seq(1,12,1), y=line1))+geom_line(aes(x=seq(1,12,1), y=line2), col="red")+
  #   geom_line(aes(x = seq(1, 12, 1), y = shifted_line2), color = "red", linetype = "dashed")+
  #   labs(title=sites[i], subtitle=paste("Synchrony=", round(corr_value$synchrony, 2),
  #                                       "\nPhase Offset=",corr_value$best_shift),
  #        x="Month", y="Normalized Concentration")+theme_classic()+
  #   theme(text = element_text(size = 20))+
  #   scale_x_continuous(breaks = seq(1,12,1))
  # 
  # print(p1)
  
  # Store the result
  results_N_P <- rbind(results_N_P, 
                        data.frame(Stream_Name = N_P_sites[i],
                                   Solute1 = "N", Solute2 = "P", corr = corr_value))
}

results_N_P$shift<-ifelse(results_N_P$corr.best_shift <=6, results_N_P$corr.best_shift, abs(12-results_N_P$corr.best_shift))

# Create a regular 0–1 grid
grid <- expand.grid(
  corr.shape_cor = seq(0, 1, length.out = 200),
  corr.amp_ratio = seq(0, 1, length.out = 200)
)

line_df <- data.frame(
  corr.shape_cor = seq(0.001, 1, length.out = 200),
  corr.amp_ratio = 0.5 / seq(0.001, 1, length.out = 200)
)

line_df <- subset(line_df, corr.amp_ratio <= 1)

# Define what the background gradient should represent:
# For example, a simple function of both axes:
grid$corr.synchrony <- (grid$corr.shape_cor*grid$corr.amp_ratio)
# (you can replace this with another function, e.g., corr.shape_cor * corr.amp_ratio)

z3<-ggplot() +
  geom_raster(data = grid, mapping=aes(corr.shape_cor, corr.amp_ratio, fill = corr.synchrony)) +
  geom_point(results_N_P, mapping=aes(corr.shape_cor, corr.amp_ratio), 
             color = "black", size=0.5) +
  geom_line(data = line_df, aes(corr.shape_cor, corr.amp_ratio),
            color = "black") +
  scale_fill_gradientn(
    colours = c("grey", "grey", "blue"),  # from blue to grey
    values = c(0, 0.48, 1),        # everything above 0.5 stays grey
    limits = c(0, 1),
    name = "Synchrony"
  ) +
  coord_fixed()+
  theme_classic()+
  labs(x="Correlation Coefficient", y="")+
  ggtitle("N-P")+
  theme(text = element_text(size = 20))

z3

pdf("Synchrony_Gradient_NoPoints.pdf", width = 8, height = 6)

ggplot() +
  geom_raster(data = grid, mapping=aes(corr.shape_cor, corr.amp_ratio, fill = corr.synchrony)) +
  geom_line(data = line_df, aes(corr.shape_cor, corr.amp_ratio),
            color = "black") +
  scale_fill_gradientn(
    colours = c("grey", "grey", "blue"),  # from blue to grey
    values = c(0, 0.48, 1),        # everything above 0.5 stays grey
    limits = c(0, 1),
    name = "Synchrony"
  ) +
  coord_fixed()+
  theme_classic()+
  labs(x="Correlation Coefficient", y="Amplitude Ratio")+
  theme(text = element_text(size = 20))

dev.off()

pdf("Corr_Amp_Synchrony.pdf", width = 15, height = 6)

ggarrange(z1, z2, z3, nrow = 1, widths = c(0.4, 0.4, 0.54))

dev.off()

p3<-ggplot(results_N_P, aes(shift, corr.synchrony))+
  geom_hline(yintercept = mean(results_N_P$corr.synchrony), col="red")+
  #geom_hline(yintercept = 0.85, col="grey")+
  #geom_vline(xintercept = .5, col="grey")+geom_vline(xintercept = 5.5, col="grey")+
  geom_boxplot(aes(group=shift), fill=NA)+geom_jitter()+
  theme_classic()+theme(text = element_text(size = 20))+
  labs(x="Shift (Months)", y="")+ggtitle("N-P")
p3

results_N_P %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous")) %>%
  ggplot(aes(x=1))+geom_bar(aes(fill=sync_char))+theme_classic()+labs(x="", y="Number of Sites")+
  theme(text = element_text(size = 20))+ggtitle("N-P")

p6<-ggplot(results_N_P, aes(shift))+
  geom_bar(fill="blue")+
  theme_classic()+theme(text = element_text(size = 20))+
  labs(x="Shift (Months)", y="Number of Sites")
p6


ggarrange(p1, p2, p3, p4, p5, p6, nrow=2, ncol = 3, align = "hv")

results_N_P<-results_N_P %>%
  mutate(shift_class = case_when(
    shift==0 ~ "in phase",
    shift > 5 ~ "out of phase",
    .default = "offset"
  ),
  sync_class = case_when(
    1-ASI.JS_min >=0.85~"synchronous",
    .default = "asynchronous"
  ),
  sync_shift=paste0(sync_class, ", ", shift_class)
  )

prop_table_N_P <- results_N_P %>%
  count(sync_class, shift_class) %>%
  mutate(prop = n / sum(n))

p6<-ggplot(prop_table_N_P, aes(y = sync_class, x = shift_class, fill = prop)) +
  geom_tile(color = "white", alpha=0.8) +
  geom_text(aes(label = scales::percent(prop, accuracy = .1)), color = "black", size = 5) +
  scale_fill_viridis_c(option = "plasma", name = "Proportion", limits=c(0,0.5)) +
  theme_minimal(base_size = 18) +
  theme(legend.position = "null", axis.text.y = element_blank())+
  labs(x = "", y = "")

p6

results_all<-bind_rows(results_Si_N, results_Si_P, results_N_P)

plot_a<-results_all %>%
  mutate(solute_pair=paste0(Solute1, "-", Solute2)) %>%
  ggplot(aes(solute_pair, corr.synchrony))+
  geom_hline(yintercept = 0.5, col="red")+
  geom_boxplot(fill=NA)+geom_jitter(size=0.5)+
  theme_classic()+labs(x="", y="Synchrony")+
  theme(text = element_text(size = 20))

plot_a

plot_b<-results_all %>%
  mutate(solute_pair=paste0(Solute1, "-", Solute2)) %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous")) %>%
  ggplot(aes(solute_pair))+geom_bar(aes(fill=sync_char), position = "fill")+
  theme_classic()+labs(x="", y="Proportion of Sites", fill="")+
  scale_fill_manual(values = c("grey", "blue"))+
  theme(text = element_text(size = 20), legend.position = "top")

plot_b

plot_c<-results_all %>%
  mutate(solute_pair=paste0(Solute1, "-", Solute2)) %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous")) %>%
  filter(sync_char=="synchronous") %>%
  ggplot(aes(shift))+geom_bar(fill="blue")+
  theme_classic()+labs(x="Shift (Months)", y="Number of Sites", fill="")+
  ggtitle("Distribution of Phase Offset for Synchronous Sites")+
  scale_x_continuous(labels = seq(0,6,1), breaks = seq(0,6,1))+
  theme(text = element_text(size = 20))+
  facet_wrap(~solute_pair)

plot_c

set1<-ggarrange(plot_a, plot_b)

pdf("Synchrony_ThreePanel_Plot.pdf", width = 10, height = 10)

ggarrange(set1, plot_c, nrow = 2, heights = c(0.5, 0.7))

dev.off()

write.csv(results_all, "Synchrony_Shift_AllSolutePairs.csv")
