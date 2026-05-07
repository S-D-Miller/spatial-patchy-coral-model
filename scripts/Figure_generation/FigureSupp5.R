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
hysteresis_outputs <- read.csv("data\\nocluster_outputs_default_params.csv")
hysteresis_outputs$unique_id_rep <- paste(hysteresis_outputs$prop_cluster, hysteresis_outputs$prop_structure, hysteresis_outputs$rep, sep = "_")
hysteresis_outputs_cluster <- read.csv("data\\cluster_outputs_default_params.csv")
hysteresis_outputs_cluster$unique_id_rep <- paste(hysteresis_outputs_cluster$prop_cluster, hysteresis_outputs_cluster$prop_structure, hysteresis_outputs_cluster$rep, sep = "_")
hysteresis_outputs <- rbind(hysteresis_outputs, hysteresis_outputs_cluster)

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

##############
##############
#MAKE PLOTS
##############
##############
#Fig S5a --> Bistable range
dev.new(width = 8, height = 7, noRStudioGD = TRUE)
par(oma = c(0,2,0,0), xpd = NA)
boxplot(bistable_range ~ quadrant, data = bistable_dat,
        xlab = "Amount of reef habitat", ylab = "Range of bistable region (dV)",
        cex.lab = 2, cex.axis = 1.35, names = c("100%", "75-100%", "50-75%", "25-50%", "0-25%"))
text(0.5, 2.97, "a)", cex = 2.5)
points(jitter(bistable_range) ~ jitter(plot_loc), data = bistable_dat[bistable_dat$quadrant != "0",], cex = 1.25, lwd = 2)
points(bistable_range ~ plot_loc, data = bistable_dat[bistable_dat$quadrant == "0",], cex = 1.25, lwd = 2)

#Fig S5b --> Bistable magnitude
dev.new(width = 8, height = 7, noRStudioGD = TRUE)
par(oma = c(0,2,0,0), xpd = NA)
boxplot(median_diff ~ quadrant, data = bistable_dat,
        xlab = "Amount of reef habitat", ylab = "Magnitude of bistable region",
        cex.lab = 2, cex.axis = 1.35, names = c("100%", "75-100%", "50-75%", "25-50%", "0-25%"),
        ylim = c(0,0.52))
text(0.5, 0.51, "b)", cex = 2.5)
points(jitter(median_diff) ~ jitter(plot_loc), data = bistable_dat[bistable_dat$quadrant != "0",], cex = 1.25, lwd = 2)
points(median_diff ~ plot_loc, data = bistable_dat[bistable_dat$quadrant == "0",], cex = 1.25, lwd = 2)