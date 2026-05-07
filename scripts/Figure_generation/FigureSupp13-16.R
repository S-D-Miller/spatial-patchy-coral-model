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

plot_locs <- data.frame("perc_diff_param" = as.integer(c(-20.0, -15.0, -10.0, -5.0, 5.0, 10.0, 15.0, 20.0)),
                        "plot_loc" = 1:8)

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

##Sensitivity --> coral recruitment
rc_outputs <- read.csv("data\\sensitivity_analysis_Rc.csv")
colnames(rc_outputs)[1] <- "unique_id_rep"
rc_outputs <- merge(rc_outputs, sensitivity_reps, by = "unique_id_rep")

##Sensitivity --> macro recruitment
rm_outputs <- read.csv("data\\sensitivity_analysis_Rm.csv")
colnames(rm_outputs)[1] <- "unique_id_rep"
rm_outputs <- merge(rm_outputs, sensitivity_reps, by = "unique_id_rep")

##Sensitivity --> macro production
pm_outputs <- read.csv("data\\sensitivity_analysis_pm.csv")
colnames(pm_outputs)[1] <- "unique_id_rep"
pm_outputs <- merge(pm_outputs, sensitivity_reps, by = "unique_id_rep")

##Sensitivity --> overgrowth
gamma_outputs <- read.csv("data\\sensitivity_analysis_y.csv")
colnames(gamma_outputs)[1] <- "unique_id_rep"
gamma_outputs <- merge(gamma_outputs, sensitivity_reps, by = "unique_id_rep")

##Sensitivity --> omega
omega_outputs <- read.csv("data\\sensitivity_analysis_w.csv")
colnames(omega_outputs)[1] <- "unique_id_rep"
omega_outputs <- merge(omega_outputs, sensitivity_reps, by = "unique_id_rep")

##Sensitivity --> coral mortality
dc_outputs <- read.csv("data\\sensitivity_analysis_Dc.csv")
colnames(dc_outputs)[1] <- "unique_id_rep"
dc_outputs <- merge(dc_outputs, sensitivity_reps, by = "unique_id_rep")

##Sensitivity --> macroalgal mortality
di_outputs <- read.csv("data\\sensitivity_analysis_Di.csv")
colnames(di_outputs)[1] <- "unique_id_rep"
di_outputs <- merge(di_outputs, sensitivity_reps, by = "unique_id_rep")

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
#Coral recruitment
###
equilibrium_rc <- rc_outputs %>%
  group_by(unique_id_rep, Rc) %>%
  reframe(high_C = mean(high_C_init_high_C),
          low_C = mean(low_C_init_low_C),
          high_M = mean(high_M_init_low_C),
          low_M = mean(low_M_init_high_C)) %>%
  mutate(perc_change_rc = (Rc / 0.001) * 100) %>%
  group_by(unique_id_rep) %>%
  mutate(default_high_C = high_C[Rc == 0.001],
         default_low_C = low_C[Rc == 0.001],
         default_high_M = high_M[Rc == 0.001],
         default_low_M = low_M[Rc == 0.001]) %>%
  ungroup() %>%
  mutate(perc_change_high_C = (high_C / default_high_C) * 100,
         perc_change_low_C = (low_C / default_low_C) * 100,
         perc_change_high_M = (high_M / default_high_M) * 100,
         perc_change_low_M = (low_M / default_low_M) * 100)
###
#Macroalgae recruitment
###
equilibrium_rm <- rm_outputs %>%
  group_by(unique_id_rep, Rm) %>%
  reframe(high_C = mean(high_C_init_high_C),
          low_C = mean(low_C_init_low_C),
          high_M = mean(high_M_init_low_C),
          low_M = mean(low_M_init_high_C)) %>%
  mutate(perc_change_rm = (Rm / 0.0001) * 100) %>%
  group_by(unique_id_rep) %>%
  mutate(default_high_C = high_C[Rm == 0.0001],
         default_low_C = low_C[Rm == 0.0001],
         default_high_M = high_M[Rm == 0.0001],
         default_low_M = low_M[Rm == 0.0001]) %>%
  ungroup() %>%
  mutate(perc_change_high_C = (high_C / default_high_C) * 100,
         perc_change_low_C = (low_C / default_low_C) * 100,
         perc_change_high_M = (high_M / default_high_M) * 100,
         perc_change_low_M = (low_M / default_low_M) * 100)

###
#Macroalgal production
###
equilibrium_pm <- pm_outputs %>%
  group_by(unique_id_rep, pm) %>%
  reframe(high_C = mean(high_C_init_high_C),
          low_C = mean(low_C_init_low_C),
          high_M = mean(high_M_init_low_C),
          low_M = mean(low_M_init_high_C)) %>%
  mutate(perc_change_pm = (pm / 0.5) * 100) %>%
  group_by(unique_id_rep) %>%
  mutate(default_high_C = high_C[pm == 0.5],
         default_low_C = low_C[pm == 0.5],
         default_high_M = high_M[pm == 0.5],
         default_low_M = low_M[pm == 0.5]) %>%
  ungroup() %>%
  mutate(perc_change_high_C = (high_C / default_high_C) * 100,
         perc_change_low_C = (low_C / default_low_C) * 100,
         perc_change_high_M = (high_M / default_high_M) * 100,
         perc_change_low_M = (low_M / default_low_M) * 100)

###
#Macroalgal overgrowth
###
equilibrium_gamma <- gamma_outputs %>%
  group_by(unique_id_rep, y) %>%
  reframe(high_C = mean(high_C_init_high_C),
          low_C = mean(low_C_init_low_C),
          high_M = mean(high_M_init_low_C),
          low_M = mean(low_M_init_high_C)) %>%
  mutate(perc_change_y = (y / 0.4) * 100) %>%
  group_by(unique_id_rep) %>%
  mutate(default_high_C = high_C[y == 0.4],
         default_low_C = low_C[y == 0.4],
         default_high_M = high_M[y == 0.4],
         default_low_M = low_M[y == 0.4]) %>%
  ungroup() %>%
  mutate(perc_change_high_C = (high_C / default_high_C) * 100,
         perc_change_low_C = (low_C / default_low_C) * 100,
         perc_change_high_M = (high_M / default_high_M) * 100,
         perc_change_low_M = (low_M / default_low_M) * 100)

###
#Omega
###
equilibrium_omega <- omega_outputs %>%
  group_by(unique_id_rep, w) %>%
  reframe(high_C = mean(high_C_init_high_C),
          low_C = mean(low_C_init_low_C),
          high_M = mean(high_M_init_low_C),
          low_M = mean(low_M_init_high_C)) %>%
  mutate(perc_change_w = (w / 2) * 100) %>%
  group_by(unique_id_rep) %>%
  mutate(default_high_C = high_C[w == 2],
         default_low_C = low_C[w == 2],
         default_high_M = high_M[w == 2],
         default_low_M = low_M[w == 2]) %>%
  ungroup() %>%
  mutate(perc_change_high_C = (high_C / default_high_C) * 100,
         perc_change_low_C = (low_C / default_low_C) * 100,
         perc_change_high_M = (high_M / default_high_M) * 100,
         perc_change_low_M = (low_M / default_low_M) * 100)

###
#Coral mortality
###
equilibrium_dc <- dc_outputs %>%
  group_by(unique_id_rep, Dc) %>%
  reframe(high_C = mean(high_C_init_high_C),
          low_C = mean(low_C_init_low_C),
          high_M = mean(high_M_init_low_C),
          low_M = mean(low_M_init_high_C)) %>%
  mutate(perc_change_dc = (Dc / 0.05) * 100) %>%
  group_by(unique_id_rep) %>%
  mutate(default_high_C = high_C[Dc == 0.05],
         default_low_C = low_C[Dc == 0.05],
         default_high_M = high_M[Dc == 0.05],
         default_low_M = low_M[Dc == 0.05]) %>%
  ungroup() %>%
  mutate(perc_change_high_C = (high_C / default_high_C) * 100,
         perc_change_low_C = (low_C / default_low_C) * 100,
         perc_change_high_M = (high_M / default_high_M) * 100,
         perc_change_low_M = (low_M / default_low_M) * 100)

###
#Macroalgal mortality
###
equilibrium_di <- di_outputs %>%
  group_by(unique_id_rep, Di) %>%
  reframe(high_C = mean(high_C_init_high_C),
          low_C = mean(low_C_init_low_C),
          high_M = mean(high_M_init_low_C),
          low_M = mean(low_M_init_high_C)) %>%
  mutate(perc_change_di = (Di / 0.4) * 100) %>%
  group_by(unique_id_rep) %>%
  mutate(default_high_C = high_C[Di == 0.4],
         default_low_C = low_C[Di == 0.4],
         default_high_M = high_M[Di == 0.4],
         default_low_M = low_M[Di == 0.4]) %>%
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
rc_list <- unique(rc_outputs$Rc) #Each unique value of rc
#Calculate the bistable region for the first value on the first landscape
bistable_rc <- calculate_bistable_region(rc_outputs[rc_outputs$Rc == rc_list[1],], reps[1])[[1]] %>%
  mutate(Rc = rc_list[1])

#Then loop for all other values on that landscape
for(i in 2:length(rc_list)){
  
  bistable_rc <- rbind(bistable_rc, 
                       calculate_bistable_region(rc_outputs[rc_outputs$Rc == rc_list[i],], reps[1])[[1]] %>% mutate(Rc = rc_list[i]))
  
}

#Now we loop through all value for all other landscapes
for(j in 2:length(reps)){
  for(i in 1:length(rc_list)){
    
    bistable_rc <- rbind(bistable_rc, 
                         calculate_bistable_region(rc_outputs[rc_outputs$Rc == rc_list[i],], reps[j])[[1]] %>% mutate(Rc = rc_list[i]))
  }
} #end with 88 rows...one for each parameter value (n = 11) and landscape (n = 8) combo

#Convert into percent differences from the default
bistable_rc <- bistable_rc %>%
  group_by(unique_id_rep) %>%
  mutate(perc_diff_param = round(((Rc / 0.001) * 100)-100, 1),
         perc_diff_median = ((median_diff / median_diff[Rc == 0.001])*100)-100,
         perc_diff_range = ((bistable_range / bistable_range[Rc == 0.001])*100)-100)

#Filter unused values
rc_boxplot <- bistable_rc %>%
  filter(perc_diff_param != 0,
         perc_diff_param != 25) %>%
  mutate(perc_diff_param = as.integer(perc_diff_param)) %>%
  left_join(plot_locs)

#Now we repeat for other spatial parameters
###
#Macroalgae growth
###
rm_list <- unique(rm_outputs$Rm)

#Calculate the bistable region for the first value on the first landscape
bistable_rm <- calculate_bistable_region(rm_outputs[rm_outputs$Rm == rm_list[1],], reps[1])[[1]] %>%
  mutate(macro_r = rm_list[1])

#Then loop for all other values on that landscape
for(i in 2:length(rm_list)){
  
  bistable_rm <- rbind(bistable_rm, 
                       calculate_bistable_region(rm_outputs[rm_outputs$Rm == rm_list[i],], reps[1])[[1]] %>% mutate(macro_r = rm_list[i]))
  
}
#Now we loop through all value for all other landscapes
for(j in 2:length(reps)){
  for(i in 1:length(rm_list)){
    
    bistable_rm <- rbind(bistable_rm, 
                         calculate_bistable_region(rm_outputs[rm_outputs$Rm == rm_list[i],], reps[j])[[1]] %>% mutate(macro_r = rm_list[i]))
  }
}

#Convert into percent differences from the default
bistable_rm <- bistable_rm %>%
  group_by(unique_id_rep) %>%
  mutate(perc_diff_param = round(((macro_r / 0.0001) * 100)-100, 1),
         perc_diff_median = ((median_diff / median_diff[macro_r == 0.0001])*100)-100,
         perc_diff_range = ((bistable_range / bistable_range[macro_r == 0.0001])*100)-100)

#Filter unused values
rm_boxplot <- bistable_rm %>%
  filter(perc_diff_param != 0,
         perc_diff_param != 25,
         macro_r != 0.04) %>%
  left_join(plot_locs)

###
#Macroalgal production
###
pm_list <- unique(pm_outputs$pm)

#Calculate the bistable region for the first value on the first landscape
bistable_pm <- calculate_bistable_region(pm_outputs[pm_outputs$pm == pm_list[1],], reps[1])[[1]] %>%
  mutate(pm = pm_list[1])

#Then loop for all other values on that landscape
for(i in 2:length(pm_list)){
  
  bistable_pm <- rbind(bistable_pm, 
                          calculate_bistable_region(pm_outputs[pm_outputs$pm == pm_list[i],], reps[1])[[1]] %>% mutate(pm = pm_list[i]))
  
}

#Now we loop through all value for all other landscapes
for(j in 2:length(reps)){
  for(i in 1:length(pm_list)){
    
    bistable_pm <- rbind(bistable_pm, 
                            calculate_bistable_region(pm_outputs[pm_outputs$pm == pm_list[i],], reps[j])[[1]] %>% mutate(pm = pm_list[i]))
  }
}

#Convert into percent differences from the default
bistable_pm <- bistable_pm %>%
  group_by(unique_id_rep) %>%
  mutate(perc_diff_param = round(((pm / 0.5) * 100)-100, 1),
         perc_diff_median = ((median_diff / median_diff[pm == 0.5])*100)-100,
         perc_diff_range = ((bistable_range / bistable_range[pm == 0.5])*100)-100)

#Filter unused values
pm_boxplot <- bistable_pm %>%
  filter(perc_diff_param != 0,
         perc_diff_param != 25,
         pm != 1.5) %>%
  left_join(plot_locs)

###
#Gamma
###
y_list <- unique(gamma_outputs$y)

#Calculate the bistable region for the first value on the first landscape
bistable_y <- calculate_bistable_region(gamma_outputs[gamma_outputs$y == y_list[1],], reps[1])[[1]] %>%
  mutate(y = y_list[1])

#Then loop for all other values on that landscape
for(i in 2:length(y_list)){
  
  bistable_y <- rbind(bistable_y, 
                         calculate_bistable_region(gamma_outputs[gamma_outputs$y == y_list[i],], reps[1])[[1]] %>% mutate(y = y_list[i]))
  
}

#Now we loop through all value for all other landscapes
for(j in 2:length(reps)){
  for(i in 1:length(y_list)){
    
    bistable_y <- rbind(bistable_y, 
                           calculate_bistable_region(gamma_outputs[gamma_outputs$y == y_list[i],], reps[j])[[1]] %>% mutate(y = y_list[i]))
  }
}


#Convert into percent differences from the default
bistable_y <- bistable_y %>%
  group_by(unique_id_rep) %>%
  mutate(perc_diff_param = round(((y / 0.4) * 100)-100, 1),
         perc_diff_median = ((median_diff / median_diff[y == 0.4])*100)-100,
         perc_diff_range = ((bistable_range / bistable_range[y == 0.4])*100)-100)

#Filter unused values
y_boxplot <- bistable_y %>%
  filter(perc_diff_param != 0,
         perc_diff_param != 25,
         y != 5,
         y != 2) %>%
  left_join(plot_locs)
y_boxplot$perc_diff_param[y_boxplot$y == 1] <- 200

###
#Omega
###
w_list <- unique(omega_outputs$w) #Each unique value of rms
#Calculate the bistable region for the first value on the first landscape
bistable_w <- calculate_bistable_region(omega_outputs[omega_outputs$w == w_list[1],], reps[1])[[1]] %>%
  mutate(w = w_list[1])

#Then loop for all other values on that landscape
for(i in 2:length(w_list)){
  
  bistable_w <- rbind(bistable_w, 
                        calculate_bistable_region(omega_outputs[omega_outputs$w == w_list[i],], reps[1])[[1]] %>% mutate(w = w_list[i]))
  
}

#Now we loop through all value for all other landscapes
for(j in 2:length(reps)){
  for(i in 1:length(w_list)){
    
    bistable_w <- rbind(bistable_w, 
                          calculate_bistable_region(omega_outputs[omega_outputs$w == w_list[i],], reps[j])[[1]] %>% mutate(w = w_list[i]))
  }
} #end with 88 rows...one for each parameter value (n = 11) and landscape (n = 8) combo

#Convert into percent differences from the default
bistable_w <- bistable_w %>%
  group_by(unique_id_rep) %>%
  mutate(perc_diff_param = round(((w / 2) * 100)-100, 1),
         perc_diff_median = ((median_diff / median_diff[w == 2])*100)-100,
         perc_diff_range = ((bistable_range / bistable_range[w == 2])*100)-100)

#Filter unused values
w_boxplot <- bistable_w %>%
  filter(perc_diff_param != 0,
         perc_diff_param != 25,
         perc_diff_param != 50,
         perc_diff_param != -50) %>%
  left_join(plot_locs)
w_boxplot$perc_diff_param[w_boxplot$w == 0] <- 200 #Just to make it so it plots all the way to the right

###
#Dc
###
dc_list <- unique(dc_outputs$Dc) #Each unique value of rms
#Calculate the bistable region for the first value on the first landscape
bistable_dc <- calculate_bistable_region(dc_outputs[dc_outputs$Dc == dc_list[1],], reps[1])[[1]] %>%
  mutate(dc = dc_list[1])

#Then loop for all other values on that landscape
for(i in 2:length(dc_list)){
  
  bistable_dc <- rbind(bistable_dc, 
                      calculate_bistable_region(dc_outputs[dc_outputs$Dc == dc_list[i],], reps[1])[[1]] %>% mutate(dc = dc_list[i]))
  
}

#Now we loop through all value for all other landscapes
for(j in 2:length(reps)){
  for(i in 1:length(dc_list)){
    
    bistable_dc <- rbind(bistable_dc, 
                        calculate_bistable_region(dc_outputs[dc_outputs$Dc == dc_list[i],], reps[j])[[1]] %>% mutate(dc = dc_list[i]))
  }
} #end with 88 rows...one for each parameter value (n = 11) and landscape (n = 8) combo

#Convert into percent differences from the default
bistable_dc <- bistable_dc %>%
  group_by(unique_id_rep) %>%
  mutate(perc_diff_param = round(((dc / 0.05) * 100)-100, 1),
         perc_diff_median = ((median_diff / median_diff[dc == 0.05])*100)-100,
         perc_diff_range = ((bistable_range / bistable_range[dc == 0.05])*100)-100)

#Filter unused values
dc_boxplot <- bistable_dc %>%
  filter(perc_diff_param != 0,
         perc_diff_param != 25,
         perc_diff_param != 50,
         perc_diff_param != -50) %>%
  left_join(plot_locs)
#dc_boxplot$perc_diff_param[w_boxplot$dc == 0] <- 200 #Just to make it so it plots all the way to the right

###
#Di
###
di_list <- unique(di_outputs$Di) #Each unique value of rms
#Calculate the bistable region for the first value on the first landscape
bistable_di <- calculate_bistable_region(di_outputs[di_outputs$Di == di_list[1],], reps[1])[[1]] %>%
  mutate(di = di_list[1])

#Then loop for all other values on that landscape
for(i in 2:length(di_list)){
  
  bistable_di <- rbind(bistable_di, 
                       calculate_bistable_region(di_outputs[di_outputs$Di == di_list[i],], reps[1])[[1]] %>% mutate(di = di_list[i]))
  
}

#Now we loop through all value for all other landscapes
for(j in 2:length(reps)){
  for(i in 1:length(di_list)){
    
    bistable_di <- rbind(bistable_di, 
                         calculate_bistable_region(di_outputs[di_outputs$Di == di_list[i],], reps[j])[[1]] %>% mutate(di = di_list[i]))
  }
} #end with 88 rows...one for each parameter value (n = 11) and landscape (n = 8) combo

#Convert into percent differences from the default
bistable_di <- bistable_di %>%
  group_by(unique_id_rep) %>%
  mutate(perc_diff_param = round(((di / 0.4) * 100)-100, 1),
         perc_diff_median = ((median_diff / median_diff[di == 0.4])*100)-100,
         perc_diff_range = ((bistable_range / bistable_range[di == 0.4])*100)-100)

#Filter unused values
di_boxplot <- bistable_di %>%
  filter(perc_diff_param != 0,
         perc_diff_param != 25,
         perc_diff_param != 50,
         perc_diff_param != -50) %>%
  left_join(plot_locs)
#dc_boxplot$perc_diff_param[w_boxplot$dc == 0] <- 200 #Just to make it so it plots all the way to the right


##############
##############
#MAKE PLOTS
##############
##############
#######
#######
#Coral recruitment
#######
#######
#Change in bistable magnitude
dev.new(width = 10, height = 4.5, noRStudioGD = TRUE)
boxplot(perc_diff_median ~ perc_diff_param, data = rc_boxplot,
        names = c("-20%", "-15%", "-10%", "-5%", "5%", "10%", "15%", "20%"),
        xlab = "Percent change in parameter", ylab = "Percent change in bistable magnitude", outline = F)
points(perc_diff_median ~ jitter(plot_loc), data = rc_boxplot, cex = 1.25, lwd = 2)

#Change in bistable range
dev.new(width = 10, height = 4.5, noRStudioGD = TRUE)
boxplot(perc_diff_range ~ perc_diff_param, data = rc_boxplot,
        names = c("-20%", "-15%", "-10%", "-5%", "5%", "10%", "15%", "20%"),
        xlab = "Percent change in parameter", ylab = "Percent change in bistable range", outline = F, ylim = c(-25, 20))
points(perc_diff_range ~ jitter(plot_loc), data = rc_boxplot, cex = 1.25, lwd = 2)

#######
#######
#Macroalgae recruitment
#######
#######
#Change in bistable magnitude
dev.new(width = 10, height = 4.5, noRStudioGD = TRUE)
boxplot(perc_diff_median ~ perc_diff_param, data = rm_boxplot,
        names = c("-20%", "-15%", "-10%", "-5%", "5%", "10%", "15%", "20%"),
        xlab = "Percent change in parameter", ylab = "Percent change in bistable magnitude", outline = F)
points(perc_diff_median ~ jitter(plot_loc), data = rm_boxplot, cex = 1.25, lwd = 2)

#Change in bistable range
dev.new(width = 10, height = 4.5, noRStudioGD = TRUE)
boxplot(perc_diff_range ~ perc_diff_param, data = rm_boxplot,
        names = c("-20%", "-15%", "-10%", "-5%", "5%", "10%", "15%", "20%"),
        xlab = "Percent change in parameter", ylab = "Percent change in bistable range", outline = F)
points(perc_diff_range ~ jitter(plot_loc), data = rm_boxplot, cex = 1.25, lwd = 2)

#######
#######
#Macroalgal productivity
#######
#######
#Change in bistable magnitude
dev.new(width = 10, height = 4.5, noRStudioGD = TRUE)
boxplot(perc_diff_median ~ perc_diff_param, data = pm_boxplot,
        names = c("-20%", "-15%", "-10%", "-5%", "5%", "10%", "15%", "20%"),
        xlab = "Percent change in parameter", ylab = "Percent change in bistable magnitude", outline = F)
points(perc_diff_median ~ jitter(plot_loc), data = pm_boxplot, cex = 1.25, lwd = 2)

#Change in bistable range
dev.new(width = 10, height = 4.5, noRStudioGD = TRUE)
boxplot(perc_diff_range ~ perc_diff_param, data = pm_boxplot,
        names = c("-20%", "-15%", "-10%", "-5%", "5%", "10%", "15%", "20%"),
        xlab = "Percent change in parameter", ylab = "Percent change in bistable range", outline = F)
points(perc_diff_range ~ jitter(plot_loc), data = pm_boxplot, cex = 1.25, lwd = 2)

#######
#######
#Gamma
#######
#######
#Change in bistable magnitude
dev.new(width = 10, height = 4.5, noRStudioGD = TRUE)
boxplot(perc_diff_median ~ perc_diff_param, data = y_boxplot,
        names = c("-20%", "-15%", "-10%", "-5%", "5%", "10%", "15%", "20%"),
        xlab = "Percent change in parameter", ylab = "Percent change in bistable magnitude", outline = F)
points(perc_diff_median ~ jitter(plot_loc), data = y_boxplot, cex = 1.25, lwd = 2)

#Change in bistable range
dev.new(width = 10, height = 4.5, noRStudioGD = TRUE)
boxplot(perc_diff_range ~ perc_diff_param, data = y_boxplot,
        names = c("-20%", "-15%", "-10%", "-5%", "5%", "10%", "15%", "20%"),
        xlab = "Percent change in parameter", ylab = "Percent change in bistable range", outline = F)
points(perc_diff_range ~ jitter(plot_loc), data = y_boxplot, cex = 1.25, lwd = 2)

#######
#######
#Omega
#######
#######
#Change in bistable magnitude
dev.new(width = 10, height = 4.5, noRStudioGD = TRUE)
boxplot(perc_diff_median ~ perc_diff_param, data = w_boxplot,
        names = c("-20%", "-15%", "-10%", "-5%", "5%", "10%", "15%", "20%"),
        xlab = "Percent change in parameter", ylab = "Percent change in bistable magnitude", outline = F)
points(perc_diff_median ~ jitter(plot_loc), data = w_boxplot, cex = 1.25, lwd = 2)

#Change in bistable range
dev.new(width = 10, height = 4.5, noRStudioGD = TRUE)
boxplot(perc_diff_range ~ perc_diff_param, data = w_boxplot,
        names = c("-20%", "-15%", "-10%", "-5%", "5%", "10%", "15%", "20%"),
        xlab = "Percent change in parameter", ylab = "Percent change in bistable range", outline = F)
points(perc_diff_range ~ jitter(plot_loc), data = w_boxplot, cex = 1.25, lwd = 2)

#######
#######
#Coral mortality
#######
#######
#Change in bistable magnitude
dev.new(width = 10, height = 4.5, noRStudioGD = TRUE)
boxplot(perc_diff_median ~ perc_diff_param, data = dc_boxplot,
        names = c("-20%", "-15%", "-10%", "-5%", "5%", "10%", "15%", "20%"),
        xlab = "Percent change in parameter", ylab = "Percent change in bistable magnitude", outline = F)
points(perc_diff_median ~ jitter(plot_loc), data = dc_boxplot, cex = 1.25, lwd = 2)

#Change in bistable range
dev.new(width = 10, height = 4.5, noRStudioGD = TRUE)
boxplot(perc_diff_range ~ perc_diff_param, data = dc_boxplot,
        names = c("-20%", "-15%", "-10%", "-5%", "5%", "10%", "15%", "20%"),
        xlab = "Percent change in parameter", ylab = "Percent change in bistable range", outline = F)
points(perc_diff_range ~ jitter(plot_loc), data = dc_boxplot, cex = 1.25, lwd = 2)

#######
#######
#Macroalgal mortality
#######
#######
#Change in bistable magnitude
#We need to make a new variable for the boxplots where we set NA values to -999 (so they plot, albeit outside of the visible range)
di_plot <- di_boxplot %>%
  mutate(perc_diff_median = ifelse(is.na(perc_diff_median), -999, perc_diff_median)) %>%
  mutate(perc_diff_range = ifelse(is.na(perc_diff_range), -999, perc_diff_range))
di_plot$perc_diff_median[di_plot$di < 0.37] <- -999
di_plot$perc_diff_range[di_plot$di < 0.37] <- -999
  
dev.new(width = 10, height = 4.5, noRStudioGD = TRUE)
boxplot(perc_diff_median ~ perc_diff_param, data = di_plot,
       names = c("-20%", "-15%", "-10%", "-5%", "5%", "10%", "15%", "20%"),
        xlab = "Percent change in parameter", ylab = "Percent change in bistable magnitude", outline = F,
       ylim = c(-65, 65))
#And we plot the original points values
points(perc_diff_median ~ jitter(plot_loc), data = di_boxplot, cex = 1.25, lwd = 2)

#Change in bistable range
dev.new(width = 10, height = 4.5, noRStudioGD = TRUE)
boxplot(perc_diff_range ~ perc_diff_param, data = di_plot,
        names = c("-20%", "-15%", "-10%", "-5%", "5%", "10%", "15%", "20%"),
        xlab = "Percent change in parameter", ylab = "Percent change in bistable range", outline = F,
        ylim = c(-100,50))
points(perc_diff_range ~ jitter(plot_loc), data = di_boxplot, cex = 1.25, lwd = 2)