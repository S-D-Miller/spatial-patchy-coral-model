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
#Sensitivity metadata
sensitivity_reps <- read.csv("data\\sensitivity_landscape_metadata.csv")
reps <- sensitivity_reps$unique_id_rep

plot_locs <- data.frame("perc_diff_param" = as.integer(c(-20.0, -15.0, -10.0, -5.0, 5.0, 10.0, 15.0, 20.0, 1000)),
                        "plot_loc" = 1:9)

#DEFAULT VALUES
hysteresis_outputs <- read.csv("data\\nocluster_outputs_default_params.csv")
hysteresis_outputs$unique_id_rep <- paste(hysteresis_outputs$prop_cluster, hysteresis_outputs$prop_structure, hysteresis_outputs$rep, sep = "_")
hysteresis_outputs_cluster <- read.csv("data\\cluster_outputs_default_params.csv")
hysteresis_outputs_cluster$unique_id_rep <- paste(hysteresis_outputs_cluster$prop_cluster, hysteresis_outputs_cluster$prop_structure, hysteresis_outputs_cluster$rep, sep = "_")

default_outputs <- rbind(hysteresis_outputs, hysteresis_outputs_cluster) %>%
  dplyr::filter(unique_id_rep %in% sensitivity_reps$unique_id_rep) %>%
  mutate(alpha = 0.75,
         beta = 3,
         coral_r = 0.0165,
         macro_r = 0.013)

##Sensitivity --> coral growth
rc_outputs <- read.csv("data\\sensitivity_outputs_coral_r.csv")
colnames(rc_outputs)[1] <- "unique_id_rep"
rc_outputs <- merge(rc_outputs, sensitivity_reps, by = "unique_id_rep")

##Sensitivity --> macro growth
rm_outputs <- read.csv("data\\sensitivity_outputs_macro_r.csv")
colnames(rm_outputs)[1] <- "unique_id_rep"
rm_outputs <- merge(rm_outputs, sensitivity_reps, by = "unique_id_rep")

##Sensitivity --> alpha
alpha_outputs <- read.csv("data\\sensitivity_outputs_alpha.csv")
colnames(alpha_outputs)[1] <- "unique_id_rep"
alpha_outputs <- merge(alpha_outputs, sensitivity_reps, by = "unique_id_rep")

##Sensitivity --> beta
beta_outputs <- read.csv("data\\sensitivity_outputs_beta.csv")
colnames(beta_outputs)[1] <- "unique_id_rep"
beta_outputs <- merge(beta_outputs, sensitivity_reps, by = "unique_id_rep")

##Sensitivity --> rms
rms_outputs <- read.csv("data\\sensitivity_outputs_rms.csv")
colnames(rms_outputs)[1] <- "unique_id_rep"
rms_outputs <- merge(rms_outputs, sensitivity_reps, by = "unique_id_rep")

##############
##############
#DATA PROCESSING
##############
##############
#######
#######
#Start by calculating equilibria for each landscape for each value of changed parameter
#######
#######
###
#Coral growth
###
equilibrium_rc <- rc_outputs %>%
  group_by(unique_id_rep, coral_r) %>%
  reframe(high_C = mean(high_C_init_high_C),
          low_C = mean(low_C_init_low_C),
          high_M = mean(high_M_init_low_C),
          low_M = mean(low_M_init_high_C)) %>%
  mutate(perc_change_r = (coral_r / 0.0165) * 100) %>%
  group_by(unique_id_rep) %>%
  mutate(default_high_C = high_C[coral_r == 0.0165],
         default_low_C = low_C[coral_r == 0.0165],
         default_high_M = high_M[coral_r == 0.0165],
         default_low_M = low_M[coral_r == 0.0165]) %>%
  ungroup() %>%
  mutate(perc_change_high_C = (high_C / default_high_C) * 100,
         perc_change_low_C = (low_C / default_low_C) * 100,
         perc_change_high_M = (high_M / default_high_M) * 100,
         perc_change_low_M = (low_M / default_low_M) * 100)
###
#Macroalgae growth
###
equilibrium_rm <- rm_outputs %>%
  group_by(unique_id_rep, macro_r) %>%
  reframe(high_C = mean(high_C_init_high_C),
          low_C = mean(low_C_init_low_C),
          high_M = mean(high_M_init_low_C),
          low_M = mean(low_M_init_high_C)) %>%
  mutate(perc_change_r = (macro_r / 0.0130) * 100) %>%
  group_by(unique_id_rep) %>%
  mutate(default_high_C = high_C[macro_r == 0.0130],
         default_low_C = low_C[macro_r == 0.0130],
         default_high_M = high_M[macro_r == 0.0130],
         default_low_M = low_M[macro_r == 0.0130]) %>%
  ungroup() %>%
  mutate(perc_change_high_C = (high_C / default_high_C) * 100,
         perc_change_low_C = (low_C / default_low_C) * 100,
         perc_change_high_M = (high_M / default_high_M) * 100,
         perc_change_low_M = (low_M / default_low_M) * 100)

###
#Alpha
###
equilibrium_alpha <- alpha_outputs %>%
  group_by(unique_id_rep, alpha) %>%
  reframe(high_C = mean(high_C_init_high_C),
          low_C = mean(low_C_init_low_C),
          high_M = mean(high_M_init_low_C),
          low_M = mean(low_M_init_high_C)) %>%
  mutate(perc_change_alpha = (alpha / 0.75) * 100) %>%
  group_by(unique_id_rep) %>%
  mutate(default_high_C = high_C[alpha == 0.75],
         default_low_C = low_C[alpha == 0.75],
         default_high_M = high_M[alpha == 0.75],
         default_low_M = low_M[alpha == 0.75]) %>%
  ungroup() %>%
  mutate(perc_change_high_C = (high_C / default_high_C) * 100,
         perc_change_low_C = (low_C / default_low_C) * 100,
         perc_change_high_M = (high_M / default_high_M) * 100,
         perc_change_low_M = (low_M / default_low_M) * 100)

###
#Beta
###
equilibrium_beta <- beta_outputs %>%
  group_by(unique_id_rep, beta) %>%
  reframe(high_C = mean(high_C_init_high_C),
          low_C = mean(low_C_init_low_C),
          high_M = mean(high_M_init_low_C),
          low_M = mean(low_M_init_high_C)) %>%
  mutate(perc_change_beta = (beta / 3) * 100) %>%
  group_by(unique_id_rep) %>%
  mutate(default_high_C = high_C[beta == 3],
         default_low_C = low_C[beta == 3],
         default_high_M = high_M[beta == 3],
         default_low_M = low_M[beta == 3]) %>%
  ungroup() %>%
  mutate(perc_change_high_C = (high_C / default_high_C) * 100,
         perc_change_low_C = (low_C / default_low_C) * 100,
         perc_change_high_M = (high_M / default_high_M) * 100,
         perc_change_low_M = (low_M / default_low_M) * 100)

###
#Rms
###
equilibrium_rms <- rms_outputs %>%
  group_by(unique_id_rep, rms) %>%
  reframe(high_C = mean(high_C_init_high_C),
          low_C = mean(low_C_init_low_C),
          high_M = mean(high_M_init_low_C),
          low_M = mean(low_M_init_high_C)) %>%
  mutate(perc_change_r = (rms / 0.5) * 100) %>%
  group_by(unique_id_rep) %>%
  mutate(default_high_C = high_C[rms == 0.5],
         default_low_C = low_C[rms == 0.5],
         default_high_M = high_M[rms == 0.5],
         default_low_M = low_M[rms == 0.5]) %>%
  ungroup() %>%
  mutate(perc_change_high_C = (high_C / default_high_C) * 100,
         perc_change_low_C = (low_C / default_low_C) * 100,
         perc_change_high_M = (high_M / default_high_M) * 100,
         perc_change_low_M = (low_M / default_low_M) * 100)


#######
#######
#Calculate the range and magnitude of bistability for each landscape under each parameter value
#######
#######
###
#Coral growth
###
rc_list <- unique(rc_outputs$coral_r) #Each unique value of rc
#Calculate the bistable region for the first value on the first landscape
bistable_rc <- calculate_bistable_region(rc_outputs[rc_outputs$coral_r == rc_list[1],], reps[1])[[1]] %>%
  mutate(coral_r = rc_list[1])

#Then loop for all other values on that landscape
for(i in 2:length(rc_list)){
  
  bistable_rc <- rbind(bistable_rc, 
                        calculate_bistable_region(rc_outputs[rc_outputs$coral_r == rc_list[i],], reps[1])[[1]] %>% mutate(coral_r = rc_list[i]))
  
}

#Now we loop through all value for all other landscapes
for(j in 2:length(reps)){
  for(i in 1:length(rc_list)){
    
    bistable_rc <- rbind(bistable_rc, 
                          calculate_bistable_region(rc_outputs[rc_outputs$coral_r == rc_list[i],], reps[j])[[1]] %>% mutate(coral_r = rc_list[i]))
  }
} #end with 88 rows...one for each parameter value (n = 11) and landscape (n = 8) combo

#Convert into percent differences from the default
bistable_rc <- bistable_rc %>%
  group_by(unique_id_rep) %>%
  mutate(perc_diff_param = round(((coral_r / 0.0165) * 100)-100, 1),
         perc_diff_median = ((median_diff / median_diff[coral_r == 0.0165])*100)-100,
         perc_diff_range = ((bistable_range / bistable_range[coral_r == 0.0165])*100)-100)

#Filter unused values
coral_r_boxplot <- bistable_rc %>%
  filter(perc_diff_param != 0,
         perc_diff_param != 25) %>%
  left_join(plot_locs)
coral_r_boxplot$plot_loc[is.na(coral_r_boxplot$plot_loc)] <- 9

#Now we repeat for other spatial parameters
###
#Macroalgae growth
###
rm_list <- unique(rm_outputs$macro_r)

#Calculate the bistable region for the first value on the first landscape
bistable_rm <- calculate_bistable_region(rm_outputs[rm_outputs$macro_r == rm_list[1],], reps[1])[[1]] %>%
  mutate(macro_r = rm_list[1])

#Then loop for all other values on that landscape
for(i in 2:length(rm_list)){
  
  bistable_rm <- rbind(bistable_rm, 
                        calculate_bistable_region(rm_outputs[rm_outputs$macro_r == rm_list[i],], reps[1])[[1]] %>% mutate(macro_r = rm_list[i]))
  
}
#Now we loop through all value for all other landscapes
for(j in 2:length(reps)){
  for(i in 1:length(rm_list)){
    
    bistable_rm <- rbind(bistable_rm, 
                          calculate_bistable_region(rm_outputs[rm_outputs$macro_r == rm_list[i],], reps[j])[[1]] %>% mutate(macro_r = rm_list[i]))
  }
}

#Convert into percent differences from the default
bistable_rm <- bistable_rm %>%
  group_by(unique_id_rep) %>%
  mutate(perc_diff_param = round(((macro_r / 0.0130) * 100)-100, 1),
         perc_diff_median = ((median_diff / median_diff[macro_r == 0.0130])*100)-100,
         perc_diff_range = ((bistable_range / bistable_range[macro_r == 0.0130])*100)-100)

#Filter unused values
macro_r_boxplot <- bistable_rm %>%
  filter(perc_diff_param != 0,
         perc_diff_param != 25,
         macro_r != 0.04) %>%
  left_join(plot_locs)
macro_r_boxplot$plot_loc[is.na(macro_r_boxplot$plot_loc)] <- 9

###
#Alpha
###
alpha_list <- unique(alpha_outputs$alpha)

#Calculate the bistable region for the first value on the first landscape
bistable_alpha <- calculate_bistable_region(alpha_outputs[alpha_outputs$alpha == alpha_list[1],], reps[1])[[1]] %>%
  mutate(alpha = alpha_list[1])

#Then loop for all other values on that landscape
for(i in 2:length(alpha_list)){
  
  bistable_alpha <- rbind(bistable_alpha, 
                        calculate_bistable_region(alpha_outputs[alpha_outputs$alpha == alpha_list[i],], reps[1])[[1]] %>% mutate(alpha = alpha_list[i]))
  
}

#Now we loop through all value for all other landscapes
for(j in 2:length(reps)){
  for(i in 1:length(alpha_list)){
    
    bistable_alpha <- rbind(bistable_alpha, 
                          calculate_bistable_region(alpha_outputs[alpha_outputs$alpha == alpha_list[i],], reps[j])[[1]] %>% mutate(alpha = alpha_list[i]))
  }
}

#Convert into percent differences from the default
bistable_alpha <- bistable_alpha %>%
  group_by(unique_id_rep) %>%
  mutate(perc_diff_param = round(((alpha / 0.75) * 100)-100, 1),
         perc_diff_median = ((median_diff / median_diff[alpha == 0.75])*100)-100,
         perc_diff_range = ((bistable_range / bistable_range[alpha == 0.75])*100)-100)

#Filter unused values
alpha_boxplot <- bistable_alpha %>%
  filter(perc_diff_param != 0,
         perc_diff_param != 25,
         alpha != 1.5) %>%
  left_join(plot_locs)
alpha_boxplot$plot_loc[is.na(alpha_boxplot$plot_loc)] <- 9

###
#Beta
###
beta_list <- unique(beta_outputs$beta)

#Calculate the bistable region for the first value on the first landscape
bistable_beta <- calculate_bistable_region(beta_outputs[beta_outputs$beta == beta_list[1],], reps[1])[[1]] %>%
  mutate(beta = beta_list[1])

#Then loop for all other values on that landscape
for(i in 2:length(beta_list)){
  
  bistable_beta <- rbind(bistable_beta, 
                        calculate_bistable_region(beta_outputs[beta_outputs$beta == beta_list[i],], reps[1])[[1]] %>% mutate(beta = beta_list[i]))
  
}

#Now we loop through all value for all other landscapes
for(j in 2:length(reps)){
  for(i in 1:length(beta_list)){
    
    bistable_beta <- rbind(bistable_beta, 
                          calculate_bistable_region(beta_outputs[beta_outputs$beta == beta_list[i],], reps[j])[[1]] %>% mutate(beta = beta_list[i]))
  }
}


#Convert into percent differences from the default
bistable_beta <- bistable_beta %>%
  group_by(unique_id_rep) %>%
  mutate(perc_diff_param = round(((beta / 3) * 100)-100, 1),
         perc_diff_median = ((median_diff / median_diff[beta == 3])*100)-100,
         perc_diff_range = ((bistable_range / bistable_range[beta == 3])*100)-100)

#Filter unused values
beta_boxplot <- bistable_beta %>%
  filter(perc_diff_param != 0,
         perc_diff_param != 25,
         beta != 5,
         beta != 2) %>%
  left_join(plot_locs)
beta_boxplot$plot_loc[is.na(beta_boxplot$plot_loc)] <- 9
beta_boxplot$perc_diff_param[beta_boxplot$beta == 1] <- 200

###
#Rms
###
rms_list <- unique(rms_outputs$rms) #Each unique value of rms
#Calculate the bistable region for the first value on the first landscape
bistable_rms <- calculate_bistable_region(rms_outputs[rms_outputs$rms == rms_list[1],], reps[1])[[1]] %>%
  mutate(rms = rms_list[1])

#Then loop for all other values on that landscape
for(i in 2:length(rms_list)){
  
  bistable_rms <- rbind(bistable_rms, 
                         calculate_bistable_region(rms_outputs[rms_outputs$rms == rms_list[i],], reps[1])[[1]] %>% mutate(rms = rms_list[i]))
  
}

#Now we loop through all value for all other landscapes
for(j in 2:length(reps)){
  for(i in 1:length(rms_list)){
    
    bistable_rms <- rbind(bistable_rms, 
                           calculate_bistable_region(rms_outputs[rms_outputs$rms == rms_list[i],], reps[j])[[1]] %>% mutate(rms = rms_list[i]))
  }
} #end with 88 rows...one for each parameter value (n = 11) and landscape (n = 8) combo

#Convert into percent differences from the default
bistable_rms <- bistable_rms %>%
  group_by(unique_id_rep) %>%
  mutate(perc_diff_param = round(((rms / 0.5) * 100)-100, 1),
         perc_diff_median = ((median_diff / median_diff[rms == 0.5])*100)-100,
         perc_diff_range = ((bistable_range / bistable_range[rms == 0.5])*100)-100)

#Filter unused values
rms_boxplot <- bistable_rms %>%
  filter(perc_diff_param != 0,
         perc_diff_param != 25,
         perc_diff_param != 50,
         perc_diff_param != -50) %>%
  left_join(plot_locs)
rms_boxplot$plot_loc[is.na(rms_boxplot$plot_loc)] <- 9
rms_boxplot$perc_diff_param[rms_boxplot$rms == 0] <- 200 #Just to make it so it plots all the way to the right


##############
##############
#MAKE PLOTS
##############
##############
#######
#######
#Coral growth
#######
#######
#Change in bistable magnitude
dev.new(width = 10, height = 4.5, noRStudioGD = TRUE)
boxplot(perc_diff_median ~ perc_diff_param, data = coral_r_boxplot,
        names = c("-20%", "-15%", "-10%", "-5%", "5%", "10%", "15%", "20%", "Extreme"),
        xlab = "Percent change in sc", ylab = "Percent change in bistable magnitude", outline = F,
        ylim = c(-35,22))
points(perc_diff_median ~ jitter(plot_loc), data = coral_r_boxplot, cex = 1.25, lwd = 2)

#Change in bistable range
dev.new(width = 10, height = 4.5, noRStudioGD = TRUE)
boxplot(perc_diff_range ~ perc_diff_param, data = coral_r_boxplot,
        names = c("-20%", "-15%", "-10%", "-5%", "5%", "10%", "15%", "20%", "Extreme"),
        xlab = "Percent change in sc", ylab = "Percent change in bistable range", outline = F)
points(perc_diff_range ~ jitter(plot_loc), data = coral_r_boxplot, cex = 1.25, lwd = 2)

#######
#######
#Macroalgae growth
#######
#######
#Change in bistable magnitude
dev.new(width = 10, height = 4.5, noRStudioGD = TRUE)
boxplot(perc_diff_median ~ perc_diff_param, data = macro_r_boxplot,
        names = c("-20%", "-15%", "-10%", "-5%", "5%", "10%", "15%", "20%", "Extreme"),
        xlab = "Percent change in parameter", ylab = "Percent change in bistable magnitude", outline = F,
        ylim = c(-28, 20))
points(perc_diff_median ~ jitter(plot_loc), data = macro_r_boxplot, cex = 1.25, lwd = 2)

#Change in bistable range
dev.new(width = 10, height = 4.5, noRStudioGD = TRUE)
boxplot(perc_diff_range ~ perc_diff_param, data = macro_r_boxplot,
        names = c("-20%", "-15%", "-10%", "-5%", "5%", "10%", "15%", "20%", "Extreme"),
        xlab = "Percent change in parameter", ylab = "Percent change in bistable range", outline = F)
points(perc_diff_range ~ jitter(plot_loc), data = macro_r_boxplot, cex = 1.25, lwd = 2)

#######
#######
#Alpha
#######
#######
#Change in bistable magnitude
dev.new(width = 10, height = 4.5, noRStudioGD = TRUE)
boxplot(perc_diff_median ~ perc_diff_param, data = alpha_boxplot,
        names = c("-20%", "-15%", "-10%", "-5%", "5%", "10%", "15%", "20%", "Extreme"),
        xlab = "Percent change in parameter", ylab = "Percent change in bistable magnitude", outline = F)
points(perc_diff_median ~ jitter(plot_loc), data = alpha_boxplot, cex = 1.25, lwd = 2)

#Change in bistable range
dev.new(width = 10, height = 4.5, noRStudioGD = TRUE)
boxplot(perc_diff_range ~ perc_diff_param, data = alpha_boxplot,
        names = c("-20%", "-15%", "-10%", "-5%", "5%", "10%", "15%", "20%", "Extreme"),
        xlab = "Percent change in parameter", ylab = "Percent change in bistable range", outline = F)
points(perc_diff_range ~ jitter(plot_loc), data = alpha_boxplot, cex = 1.25, lwd = 2)

#######
#######
#Beta
#######
#######
#Change in bistable magnitude
dev.new(width = 10, height = 4.5, noRStudioGD = TRUE)
boxplot(perc_diff_median ~ perc_diff_param, data = beta_boxplot,
        names = c("-20%", "-15%", "-10%", "-5%", "5%", "10%", "15%", "20%", "Extreme"),
        xlab = "Percent change in parameter", ylab = "Percent change in bistable magnitude", outline = F)
points(perc_diff_median ~ jitter(plot_loc), data = beta_boxplot, cex = 1.25, lwd = 2)

#Change in bistable range
dev.new(width = 10, height = 4.5, noRStudioGD = TRUE)
boxplot(perc_diff_range ~ perc_diff_param, data = beta_boxplot,
        names = c("-20%", "-15%", "-10%", "-5%", "5%", "10%", "15%", "20%", "Extreme"),
        xlab = "Percent change in parameter", ylab = "Percent change in bistable range", outline = F)
points(perc_diff_range ~ jitter(plot_loc), data = beta_boxplot, cex = 1.25, lwd = 2)

#######
#######
#Rms
#######
#######
#Change in bistable magnitude
dev.new(width = 10, height = 4.5, noRStudioGD = TRUE)
boxplot(perc_diff_median ~ perc_diff_param, data = rms_boxplot,
        names = c("-20%", "-15%", "-10%", "-5%", "5%", "10%", "15%", "20%", "Extreme"),
        xlab = "Percent change in parameter", ylab = "Percent change in bistable magnitude", outline = F)
points(perc_diff_median ~ jitter(plot_loc), data = rms_boxplot, cex = 1.25, lwd = 2)

#Change in bistable range
dev.new(width = 10, height = 4.5, noRStudioGD = TRUE)
boxplot(perc_diff_range ~ perc_diff_param, data = rms_boxplot,
        names = c("-20%", "-15%", "-10%", "-5%", "5%", "10%", "15%", "20%", "Extreme"),
        xlab = "Percent change in parameter", ylab = "Percent change in bistable range", outline = F, ylim = c(-75, 20))
points(perc_diff_range ~ jitter(plot_loc), data = rms_boxplot, cex = 1.25, lwd = 2)