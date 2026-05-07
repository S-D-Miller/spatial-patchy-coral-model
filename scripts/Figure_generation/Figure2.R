rm(list=ls())

##############
##############
#Loads in necessary packages
##############
##############
library(tidyverse)

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

#Uniform landscape for reference lines
uniform_hysteresis <- read.csv("data\\uniform_outputs_default_params.csv")

##############
##############
#Processes data
##############
##############
output_sum_rep <- hysteresis_outputs %>%
  group_by(unique_id_rep, high_initial_var, dV) %>%
  summarize(prop_structure = mean(p_struc),
            clumpiness = mean(clump),
            prop_C = mean(prop_C),
            end_C = mean(end_C),
            prop_M = mean(prop_M),
            end_M = mean(end_Mi + end_Mv),
            prop_ratio_C = mean(prop_ratio_C),
            high_C_init_high_C = mean(high_C_init_high_C),
            high_C_init_low_C = mean(high_C_init_low_C),
            low_C_init_high_C = mean(low_C_init_high_C),
            low_C_init_low_C = mean(low_C_init_low_C),
            high_M_init_high_C = mean(high_M_init_high_C),
            high_M_init_low_C = mean(high_M_init_low_C),
            low_M_init_high_C = mean(low_M_init_high_C),
            low_M_init_low_C = mean(low_M_init_low_C))

######
#Filters and summarizes landscapes with 0 - 25% habitat
######
herbivory_ranges_q1 <- output_sum_rep %>%
  filter(prop_structure <= 0.25) %>%
  group_by(high_initial_var, dV) %>%
  summarize(mean_end_C = mean(end_C),
            sd_end_C = sd(end_C),
            mean_end_M = mean(end_M),
            sd_end_M = sd(end_M),
            end_C_95h = quantile(end_C, probs = 0.975),
            end_C_95l = quantile(end_C, probs = 0.025))

######
#Filters and summarizes landscapes with 25-50% habitat
######
herbivory_ranges_q2 <- output_sum_rep %>%
  filter(prop_structure > 0.25 & prop_structure <= 0.5) %>%
  group_by(high_initial_var, dV) %>%
  summarize(mean_end_C = mean(end_C),
            median_end_C = median(end_C),
            sd_end_C = sd(end_C),
            mean_end_M = mean(end_M),
            median_end_M = median(end_M),
            sd_end_M = sd(end_M),
            end_C_95h = quantile(end_C, probs = 0.975),
            end_C_95l = quantile(end_C, probs = 0.025))

######
#Filters and summarizes landscapes with 50-75% habitat
######
herbivory_ranges_q3 <- output_sum_rep %>%
  filter(prop_structure > 0.5 & prop_structure <= 0.75) %>%
  group_by(high_initial_var, dV) %>%
  summarize(mean_end_C = mean(end_C),
            sd_end_C = sd(end_C),
            mean_end_M = mean(end_M),
            sd_end_M = sd(end_M),
            end_C_95h = quantile(end_C, probs = 0.975),
            end_C_95l = quantile(end_C, probs = 0.025))

######
#Filters and summarizes landscapes with 75-100% habitat
######
herbivory_ranges_q4 <- output_sum_rep %>%
  filter(prop_structure > 0.75 & prop_structure <= 1) %>%
  group_by(high_initial_var, dV) %>%
  summarize(mean_end_C = mean(end_C),
            sd_end_C = sd(end_C),
            mean_end_M = mean(end_M),
            sd_end_M = sd(end_M),
            end_C_95h = quantile(end_C, probs = 0.975),
            end_C_95l = quantile(end_C, probs = 0.025))

####
#Sets colors
####
coral_col <- "firebrick2"
macro_col <- "deepskyblue4"

###############
###############
#Figure 3a --> 75-100% habitat
###############
###############
dev.new()
par(oma = c(1,2,0,0))
#Starts by creating an empty plot
plot(uniform_hysteresis$dV[uniform_hysteresis$high_initial_var == "C"],
     uniform_hysteresis$end_C[uniform_hysteresis$high_initial_var == "C"],
     type = 'n', lwd = 2.5, ylim = c(0,0.6), xlab = "Herbivory (dV)",
     ylab = "Mean ending coral cover", main = "", cex.axis = 1.75, cex.lab = 2)
#Then adds ribbons
attach(herbivory_ranges_q4)
polygon(c(dV[high_initial_var == "C"],
          rev(dV[high_initial_var == "C"])),
          c(end_C_95h[high_initial_var == "C"], 
            rev(end_C_95l[high_initial_var == "C"])),
          col = adjustcolor(coral_col, alpha.f=0.3), border = NA)

polygon(c(dV[high_initial_var == "Mi"],
          rev(dV[high_initial_var == "Mi"])),
        c(end_C_95h[high_initial_var == "Mi"], 
          rev(end_C_95l[high_initial_var == "Mi"])),
        col = adjustcolor(macro_col, alpha.f=0.3), border = NA)
#Then plot the gray reference lines for uniform landscape
points(uniform_hysteresis$dV[uniform_hysteresis$high_initial_var == "C"],
       uniform_hysteresis$end_C[uniform_hysteresis$high_initial_var == "C"],
       type = 'l', col = 'gray65', lwd = 2.5)
points(uniform_hysteresis$dV[uniform_hysteresis$high_initial_var == "Mi"],
       uniform_hysteresis$end_C[uniform_hysteresis$high_initial_var == "Mi"],
       type = 'l', col = 'gray65', lty = 2, lwd = 2.5)
#Then plots lines for fragmented habitat
points(dV[high_initial_var == "C"],
       mean_end_C[high_initial_var == "C"],
       type = 'l', col = coral_col, lwd = 3)
points(dV[high_initial_var == "Mi"],
       mean_end_C[high_initial_var == "Mi"],
       type = 'l', col = macro_col, lwd = 3, lty = 2)
text(0.15, 0.59, labels = "a)", cex = 2.5)
legend('bottomright', col = c(coral_col, macro_col, 'gray65', 'gray65'), legend = c("High C", "High M",
                                                                                    "High C (uniform)", "High M (uniform)"),
       lty = c(1,2, 1, 2), lwd = 2.5, cex = 1)

###############
###############
#Figure 3b --> 50-75% habitat
###############
###############
dev.new()
par(oma = c(1,2,0,0))
plot(uniform_hysteresis$dV[uniform_hysteresis$high_initial_var == "C"],
     uniform_hysteresis$end_C[uniform_hysteresis$high_initial_var == "C"],
     type = 'n', lwd = 2.5, ylim = c(0,0.6), xlab = "Herbivory (dV)",
     ylab = "Mean ending coral cover", main = "", cex.axis = 1.75, cex.lab = 2)

#Then adds ribbons
attach(herbivory_ranges_q3)
polygon(c(dV[high_initial_var == "C"],
          rev(dV[high_initial_var == "C"])),
        c(end_C_95h[high_initial_var == "C"], 
          rev(end_C_95l[high_initial_var == "C"])),
        col = adjustcolor(coral_col, alpha.f=0.3), border = NA)

polygon(c(dV[high_initial_var == "Mi"],
          rev(dV[high_initial_var == "Mi"])),
        c(end_C_95h[high_initial_var == "Mi"], 
          rev(end_C_95l[high_initial_var == "Mi"])),
        col = adjustcolor(macro_col, alpha.f=0.3), border = NA)
#Then plot the gray reference lines for uniform landscape
points(uniform_hysteresis$dV[uniform_hysteresis$high_initial_var == "C"],
       uniform_hysteresis$end_C[uniform_hysteresis$high_initial_var == "C"],
       type = 'l', col = 'gray65', lwd = 2.5)
points(uniform_hysteresis$dV[uniform_hysteresis$high_initial_var == "Mi"],
       uniform_hysteresis$end_C[uniform_hysteresis$high_initial_var == "Mi"],
       type = 'l', col = 'gray65', lty = 2, lwd = 2.5)
#Then plots lines for fragmented habitat
points(dV[high_initial_var == "C"],
       mean_end_C[high_initial_var == "C"],
       type = 'l', col = coral_col, lwd = 3)
points(dV[high_initial_var == "Mi"],
       mean_end_C[high_initial_var == "Mi"],
       type = 'l', col = macro_col, lwd = 3, lty = 2)
text(0.15, 0.59, labels = "b)", cex = 2.5)

###############
###############
#Figure 3c --> 25-50% habitat
###############
###############
dev.new()
par(oma = c(1,2,0,0))
plot(uniform_hysteresis$dV[uniform_hysteresis$high_initial_var == "C"],
     uniform_hysteresis$end_C[uniform_hysteresis$high_initial_var == "C"],
     type = 'n', lwd = 2.5, ylim = c(0,0.6), xlab = "Herbivory (dV)",
     ylab = "Mean ending coral cover", main = "", cex.axis = 1.75, cex.lab = 2)

#Then adds ribbons
attach(herbivory_ranges_q2)
polygon(c(dV[high_initial_var == "C"],
          rev(dV[high_initial_var == "C"])),
        c(end_C_95h[high_initial_var == "C"], 
          rev(end_C_95l[high_initial_var == "C"])),
        col = adjustcolor(coral_col, alpha.f=0.3), border = NA)

polygon(c(dV[high_initial_var == "Mi"],
          rev(dV[high_initial_var == "Mi"])),
        c(end_C_95h[high_initial_var == "Mi"], 
          rev(end_C_95l[high_initial_var == "Mi"])),
        col = adjustcolor(macro_col, alpha.f=0.3), border = NA)
#Then plot the gray reference lines for uniform landscape
points(uniform_hysteresis$dV[uniform_hysteresis$high_initial_var == "C"],
       uniform_hysteresis$end_C[uniform_hysteresis$high_initial_var == "C"],
       type = 'l', col = 'gray65', lwd = 2.5)
points(uniform_hysteresis$dV[uniform_hysteresis$high_initial_var == "Mi"],
       uniform_hysteresis$end_C[uniform_hysteresis$high_initial_var == "Mi"],
       type = 'l', col = 'gray65', lty = 2, lwd = 2.5)
#Then plots lines for fragmented habitat
points(dV[high_initial_var == "C"],
       mean_end_C[high_initial_var == "C"],
       type = 'l', col = coral_col, lwd = 3)
points(dV[high_initial_var == "Mi"],
       mean_end_C[high_initial_var == "Mi"],
       type = 'l', col = macro_col, lwd = 3, lty = 2)
text(0.15, 0.59, labels = "c)", cex = 2.5)

###############
###############
#Figure 3d --> 0-25% habitat
###############
###############
dev.new()
par(oma = c(1,2,0,0))
plot(uniform_hysteresis$dV[uniform_hysteresis$high_initial_var == "C"],
     uniform_hysteresis$end_C[uniform_hysteresis$high_initial_var == "C"],
     type = 'n', lwd = 2.5, ylim = c(0,0.6), xlab = "Herbivory (dV)",
     ylab = "Mean ending coral cover", main = "", cex.axis = 1.75, cex.lab = 2)

#Then adds ribbons
attach(herbivory_ranges_q1)
polygon(c(dV[high_initial_var == "C"],
          rev(dV[high_initial_var == "C"])),
        c(end_C_95h[high_initial_var == "C"], 
          rev(end_C_95l[high_initial_var == "C"])),
        col = adjustcolor(coral_col, alpha.f=0.3), border = NA)

polygon(c(dV[high_initial_var == "Mi"],
          rev(dV[high_initial_var == "Mi"])),
        c(end_C_95h[high_initial_var == "Mi"], 
          rev(end_C_95l[high_initial_var == "Mi"])),
        col = adjustcolor(macro_col, alpha.f=0.3), border = NA)
#Then plot the gray reference lines for uniform landscape
points(uniform_hysteresis$dV[uniform_hysteresis$high_initial_var == "C"],
       uniform_hysteresis$end_C[uniform_hysteresis$high_initial_var == "C"],
       type = 'l', col = 'gray65', lwd = 2.5)
points(uniform_hysteresis$dV[uniform_hysteresis$high_initial_var == "Mi"],
       uniform_hysteresis$end_C[uniform_hysteresis$high_initial_var == "Mi"],
       type = 'l', col = 'gray65', lty = 2, lwd = 2.5)
#Then plots lines for fragmented habitat
points(dV[high_initial_var == "C"],
       mean_end_C[high_initial_var == "C"],
       type = 'l', col = coral_col, lwd = 3)
points(dV[high_initial_var == "Mi"],
       mean_end_C[high_initial_var == "Mi"],
       type = 'l', col = macro_col, lwd = 3, lty = 2)
text(0.15, 0.59, labels = "d)", cex = 2.5)