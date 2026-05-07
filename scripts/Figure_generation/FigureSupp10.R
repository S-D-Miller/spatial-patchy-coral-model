rm(list=ls())

##############
##############
#Loads in necessary packages
##############
##############
library(tidyverse)
source("scripts\\model_functions.R")

##############
##############
#Loads in data
##############
##############
#DEFAULT VALUES
hysteresis_outputs <- read.csv("data\\nocluster_outputs_default_params.csv")
hysteresis_outputs$unique_id_rep <- paste(hysteresis_outputs$prop_cluster, hysteresis_outputs$prop_structure, hysteresis_outputs$rep, sep = "_")
hysteresis_outputs_cluster <- read.csv("data\\cluster_outputs_default_params.csv")
hysteresis_outputs_cluster$unique_id_rep <- paste(hysteresis_outputs_cluster$prop_cluster, hysteresis_outputs_cluster$prop_structure, hysteresis_outputs_cluster$rep, sep = "_")
hysteresis_outputs <- rbind(hysteresis_outputs, hysteresis_outputs_cluster)

#DOUBLE a 
hysteresis_outputs_a15 <- read.csv("data\\nocluster_outputs_a15.csv")
hysteresis_outputs_a15$unique_id_rep <- paste(hysteresis_outputs_a15$prop_cluster, hysteresis_outputs_a15$prop_structure, hysteresis_outputs_a15$rep, sep = "_")
hysteresis_outputs_cluster_a15 <- read.csv("data\\cluster_outputs_a15.csv")
hysteresis_outputs_cluster_a15$unique_id_rep <- paste(hysteresis_outputs_cluster_a15$prop_cluster, hysteresis_outputs_cluster_a15$prop_structure, hysteresis_outputs_cluster_a15$rep, sep = "_")
hysteresis_outputs_a15 <- rbind(hysteresis_outputs_a15, hysteresis_outputs_cluster_a15)

#Uniform landscape
uniform_hysteresis <- read.csv("data\\uniform_outputs_default_params.csv")

#Plot locations for all data points
plot_locs <- data.frame("quadrant" = c("0", "1", "2", "3", "4"),
                        "plot_loc" = c(1, 2, 3, 4, 5))

##############
##############
#CALCULATE BISTABLE REGIONS
##############
##############
#######
#######
#Default values
#First  calculate equilibrium values so we  can scale the bistability magnitude
equilibrium_characteristics_spatial <- hysteresis_outputs %>%
  group_by(unique_id_rep) %>%
  reframe(p_struc = mean(p_struc),
          clump = mean(clump),
          high_C = mean(high_C_init_high_C),
          low_C = mean(low_C_init_low_C),
          high_M = mean(high_M_init_low_C),
          low_M = mean(low_M_init_high_C))
#Do the same for the uniform values
uniform_characteristics <- uniform_hysteresis %>%
  group_by(unique_id_rep) %>%
  reframe(p_struc = mean(p_struc),
          clump = mean(clump),
          high_C = mean(high_C_init_high_C),
          low_C = mean(low_C_init_low_C),
          high_M = mean(high_M_init_low_C),
          low_M = mean(low_M_init_high_C))

equilibrium_characteristics <- rbind(equilibrium_characteristics_spatial, uniform_characteristics)

#Now we calculate the bistable region of the uniform landscape
bistable_uniform <- calculate_bistable_region(uniform_hysteresis, unique(uniform_hysteresis$unique_id_rep))[[1]] %>%
  mutate(quadrant = 0)
#...and now for the the first landscape
bistable_dat <- calculate_bistable_region(hysteresis_outputs, unique(hysteresis_outputs$unique_id_rep)[1])[[1]]
#Then loop over all other landscapes
for(i in 2:length(unique(hysteresis_outputs$unique_id_rep))){
  
  bistable_dat <- rbind(bistable_dat, 
                        calculate_bistable_region(hysteresis_outputs, unique(hysteresis_outputs$unique_id_rep)[i])[[1]])
  
}

#Creates a "quadrant" column that converts the amount of habitat into a category based on 0-25%, 25-50%, etc.
bistable_dat <- bistable_dat %>%
  mutate(quadrant = ifelse(p_struc <= 0.25, "4",
                           ifelse(p_struc > 0.25 & p_struc <= 0.50, "3",
                                  ifelse(p_struc > 0.5 & p_struc <= 0.75, "2", "1")))) %>%
  rbind(bistable_uniform) %>%
  inner_join(equilibrium_characteristics) %>%
  mutate(prop_possible = median_diff / (high_C - low_C)) %>%
  left_join(plot_locs, by = "quadrant")

#######
#######
#Double dispersal --> a = 1.5
#First  calculate equilibrium values so we  can scale the bistability magnitude
equilibrium_characteristics_a15 <- hysteresis_outputs_a15 %>%
  group_by(unique_id_rep) %>%
  reframe(p_struc = mean(p_struc),
          clump = mean(clump),
          high_C = mean(high_C_init_high_C),
          low_C = mean(low_C_init_low_C),
          high_M = mean(high_M_init_low_C),
          low_M = mean(low_M_init_high_C))

#Do the same for the uniform values
uniform_characteristics <- uniform_hysteresis %>%
  group_by(unique_id_rep) %>%
  reframe(p_struc = mean(p_struc),
          clump = mean(clump),
          high_C = mean(high_C_init_high_C),
          low_C = mean(low_C_init_low_C),
          high_M = mean(high_M_init_low_C),
          low_M = mean(low_M_init_high_C))

#Now we calculate the bistable region of the uniform landscape
bistable_uniform <- calculate_bistable_region(uniform_hysteresis, unique(uniform_hysteresis$unique_id_rep))[[1]] %>%
  mutate(quadrant = 0)
#And the same for the spatial data
bistable_dat_a15 <- calculate_bistable_region(hysteresis_outputs_a15, unique(hysteresis_outputs_a15$unique_id_rep)[1])[[1]]
for(i in 2:length(unique(hysteresis_outputs_a15$unique_id_rep))){
  
  bistable_dat_a15 <- rbind(bistable_dat_a15, 
                            calculate_bistable_region(hysteresis_outputs_a15, unique(hysteresis_outputs_a15$unique_id_rep)[i])[[1]])
  
}

#Creates a "quadrant" column that converts the amount of habitat into a category based on 0-25%, 25-50%, etc.
bistable_dat_a15 <- bistable_dat_a15 %>%
  mutate(quadrant = ifelse(p_struc <= 0.25, "4",
                           ifelse(p_struc > 0.25 & p_struc <= 0.50, "3",
                                  ifelse(p_struc > 0.5 & p_struc <= 0.75, "2", "1")))) %>%
  rbind(bistable_uniform) %>%#combine uniform and spatial
  left_join(plot_locs, by = "quadrant") 

##############
##############
#MAKE PLOTS
##############
##############
#Figure S7a --> Bistable range a = 1.5
dev.new(width = 8, height = 7, noRStudioGD = TRUE)
par(oma = c(0,2,0,0), xpd = NA)
boxplot(bistable_range ~ quadrant, data = bistable_dat_a15,
        xlab = "Amount of reef habitat", ylab = "Range of bistable region (dV)",
        cex.lab = 2, cex.axis = 1.35, names = c("100%", "75-100%", "50-75%", "25-50%", "0-25%"))
text(0.5, 2.5, "a)", cex = 2.5)
points(jitter(bistable_range) ~ jitter(plot_loc), data = bistable_dat_a15[bistable_dat_a15$quadrant != "0",], cex = 1.25, lwd = 2)
points(bistable_range ~ plot_loc, data = bistable_dat_a15[bistable_dat_a15$quadrant == "0",], cex = 1.25, lwd = 2)

#Figure S7b --> Bistable magnitude a = 1.5
dev.new(width = 8, height = 7, noRStudioGD = TRUE)
par(oma = c(0,2,0,0), xpd = NA)
boxplot(median_diff ~ quadrant, data = bistable_dat_a15,
        xlab = "Amount of reef habitat", ylab = "Magnitude of bistable region",
        cex.lab = 2, cex.axis = 1.35, names = c("100%", "75-100%", "50-75%", "25-50%", "0-25%"),
        ylim = c(0,0.5))
text(0.5, 0.5, "b)", cex = 2.5)
points(jitter(median_diff) ~ jitter(plot_loc), data = bistable_dat_a15[bistable_dat_a15$quadrant != "0",], cex = 1.25, lwd = 2)
points(median_diff ~ plot_loc, data = bistable_dat_a15[bistable_dat_a15$quadrant == "0",], cex = 1.25, lwd = 2)

#Figure S7c --> CHANGE in range relative to default a
plot_spat <- bistable_dat[bistable_dat$quadrant != "0",]
plot_nonspat <- bistable_dat[bistable_dat$quadrant == "0",]

plot_spat_15 <- bistable_dat_a15[bistable_dat_a15$quadrant != "0",]
plot_nonspat_15 <- bistable_dat_a15[bistable_dat_a15$quadrant == "0",]

dev.new(width = 8, height = 7, noRStudioGD = TRUE)
par(oma = c(0,2,0,0), xpd = NA)
boxplot(bistable_dat_a15$bistable_range - bistable_dat$bistable_range ~ bistable_dat_a15$quadrant,
        xlab = "Amount of reef habitat", ylab = "Difference in bistable range",
        cex.lab = 2, cex.axis = 1.35, names = c("100%", "75-100%", "50-75%", "25-50%", "0-25%"))
text(0.5, 0.5, "c)", cex = 2.5)
points(jitter(plot_spat_15$bistable_range - plot_spat$bistable_range) ~ jitter(plot_spat_15$plot_loc), cex = 1.25, lwd = 2)
points(plot_nonspat$bistable_range - plot_nonspat$bistable_range ~ plot_nonspat$plot_loc, cex = 1.25, lwd = 2)

#Figure S7d --> CHANGE in magnitude relative to default a
dev.new(width = 8, height = 7, noRStudioGD = TRUE)
par(oma = c(0,2,0,0), xpd = NA)
boxplot(bistable_dat_a15$median_diff - bistable_dat$median_diff ~ bistable_dat_a15$quadrant,
        xlab = "Amount of reef habitat", ylab = "Difference in magnitude",
        cex.lab = 2, cex.axis = 1.35, names = c("100%", "75-100%", "50-75%", "25-50%", "0-25%"))
text(0.5, 0.33, "d)", cex = 2.5)
points(jitter(plot_spat_15$median_diff - plot_spat$median_diff) ~ jitter(plot_spat_15$plot_loc), cex = 1.25, lwd = 2)
points(plot_nonspat$median_diff - plot_nonspat$median_diff ~ plot_nonspat$plot_loc, cex = 1.25, lwd = 2)