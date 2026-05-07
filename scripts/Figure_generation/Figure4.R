rm(list=ls())
##############
##############
#Loads in necessary packages
##############
##############
library(tidyverse)
library(betareg)
library(pscl)
library(gplots)
library(boot)

##############
##############
#Loads in data
##############
##############

#######
#MODEL DATA
#######
hysteresis_outputs <- read.csv("data\\nocluster_outputs_default_params.csv")
hysteresis_outputs$unique_id_rep <- paste(hysteresis_outputs$prop_cluster, hysteresis_outputs$prop_structure, hysteresis_outputs$rep, sep = "_")
hysteresis_outputs_cluster <- read.csv("data\\cluster_outputs_default_params.csv")
hysteresis_outputs_cluster$unique_id_rep <- paste(hysteresis_outputs_cluster$prop_cluster, hysteresis_outputs_cluster$prop_structure, hysteresis_outputs_cluster$rep, sep = "_")
hysteresis_outputs <- rbind(hysteresis_outputs, hysteresis_outputs_cluster)

#######
#FIELD DATA
#######
dat <- read.csv("data\\benthic_transect.csv")
dat_min <- read.csv("data\\benthic_segment.csv")
dat_min$c_dom <- ifelse(dat_min$CoralStd > dat_min$AlgaeStd, 1, 0)
dat_min$m_dom <- ifelse(dat_min$AlgaeStd > dat_min$CoralStd, 1, 0)

n_segments <- dat_min %>%
  group_by(UniqueCode) %>%
  reframe(n_segments = n())

##############
##############
#Processes model data for heatmap
##############
##############
#Groups by individual rep
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
            low_M_init_low_C = mean(low_M_init_low_C)) %>%
  mutate(n_cells = prop_structure * 225,
         prop_eq_M = end_M / high_M_init_low_C,
         prop_eq_C = end_C / high_C_init_high_C,
         perc_hab = prop_structure * 100) %>%
  filter(dV < 6 & dV >= 0.25) #Gets rid of herbivory values beyond range of bistability

################
################
#6a --> PROPORTION PLOT
################
################
###
#Create heat map from model outcomes
###
gray_vals <- colorRampPalette(c('gray90', 'gray25'))(14)
col_breaks <- c(0, seq(0,0.7, by = 0.05))
legend_vals <- seq(0, 70, length.out = 15)

#Create the 2d histogram values with counts
proportion_heat <- hist2d(output_sum_rep$perc_hab, output_sum_rep$prop_C, show=FALSE, nbins = 5)

# Scale counts by habitat bin to correct for unequal outcomes in each bin
scaled_counts_proportion <- sweep(proportion_heat$counts, 1, rowSums(proportion_heat$counts), "/")

####
#Prepare field data
####
#Process data for plotting
var_plot <- dat_min %>%
  filter(!is.na(c_dom),
         !is.na(m_dom)) %>%
  group_by(UniqueCode) %>%
  summarize(prop_c = sum(c_dom) / (sum(c_dom) + sum(m_dom)),
            tot_successes = sum(c_dom),
            tot_fails = sum(m_dom),
            tot_trials = sum(c_dom) + sum(m_dom))%>%
  left_join(dat) %>%
  mutate(perc_hab = 100 - SoftSub) %>%
  mutate(quadrant_struc = as.factor(ifelse(perc_hab <= 25, "1",
                                           ifelse(perc_hab > 25 & perc_hab <= 50, "2",
                                                  ifelse(perc_hab > 50 & perc_hab <= 75, "3", "4")))))
#Run logistic regression
var_binom <- glm(prop_c ~ perc_hab, family = "binomial", data = var_plot, weights = tot_trials)
x_var <- seq(0,100,by=0.01)
#Make predictions for the line and calculate confidence intervals
pred.logit.y <- predict.glm(var_binom, newdata=list(perc_hab=x_var, type = "response"), se.fit = T)
y_var <- predict.glm(var_binom, list(perc_hab = x_var),type="response")

boot_predict_prop <- function(data, indices) {
  boot_data <- data[indices, ]
  boot_model <- glm(prop_c ~ perc_hab, family = "binomial", data = boot_data, weights = tot_trials)
  predict(boot_model, newdata=list(perc_hab=x_var, type = "response"), type = "response")
}

# Run bootstrap
boot_results_prop <- boot(var_plot, boot_predict_prop, R = 1000)

# Calculate percentile intervals
pred_intervals_prop <- t(apply(boot_results_prop$t, 2, quantile, probs = c(0.025, 0.975)))

###########
#Makes the combined plot for proportions
###########
dev.new()
par(oma=c(0,1,0,0))
#Start with heat map
image(proportion_heat$x, proportion_heat$y, scaled_counts_proportion, 
      xlab="", ylab="", 
      main="", col = c('white', gray_vals),
      breaks = col_breaks, cex.axis = 1.85,
      ylim = c(-0.05, 1.05), xlim = c(0, 100))
#Now add field data
points(var_plot$perc_hab, jitter(var_plot$prop_c), pch = 21, cex = 1.75, bg = 'deepskyblue2')
lines(x_var, y_var, lwd = 7.5, col = 'black') #adds black outline to regression line
lines(x_var, y_var, lwd = 5, col = 'deepskyblue2')
lines(x_var, pred_intervals_prop[,1], lty = 2, col = 'deepskyblue2', lwd = 3)
lines(x_var, pred_intervals_prop[,2], lty = 2, col = 'deepskyblue2', lwd = 3)
mtext("Percent of reef habitat", side = 1, cex = 2.5, line = 3)
mtext("Proportion dominated by coral", side = 2, cex = 2.5, line = 3)

################
################
#6b --> ALGAE PLOT
################
################
###
#Create heat map from model outcomes
###
gray_vals_algae <- colorRampPalette(c('gray90', 'gray25'))(14)
col_breaks_algae <- c(0, seq(0,0.7, by = 0.05))
legend_vals_algae <- seq(0, 70, length.out = 15)

#Create the 2d histogram values with counts
algae_heat <- hist2d(output_sum_rep$perc_hab, output_sum_rep$end_M, show=FALSE, nbins = 5)

# Scale counts by habitat bin to correct for unequal outcomes in each bin
scaled_counts_algae <- sweep(algae_heat$counts, 1, rowSums(algae_heat$counts), "/")

####
#Prepare field data
####
#Process data for plotting
plot_dat <- dat %>%
  mutate(perc_hab = 100 - SoftSub) %>%
  left_join(n_segments) %>% 
  mutate(quadrant_struc = as.factor(ifelse(perc_hab <= 25, "1",
                                           ifelse(perc_hab > 25 & perc_hab <= 50, "2",
                                                  ifelse(perc_hab > 50 & perc_hab <= 75, "3", "4")))))
#Run beta regression
algae_beta <- betareg(AlgaeStd ~ perc_hab, data = plot_dat)
y_algae_beta <- betareg::predict(algae_beta, data.frame("perc_hab" = x_var))

#Bootstrap confidence intervals
boot_predict_algae <- function(data, indices) {
  boot_data <- data[indices, ]
  boot_model <- betareg(AlgaeStd ~ perc_hab, data = boot_data)
  predict(boot_model, newdata = data.frame("perc_hab" = x_var), type = "response")
}

# Run bootstrap
boot_results_algae <- boot(plot_dat, boot_predict_algae, R = 1000)

# Calculate percentile intervals
pred_intervals_algae <- t(apply(boot_results_algae$t, 2, quantile, probs = c(0.025, 0.975)))

###############
#Makes combined plot for algae
###############
dev.new()
par(oma=c(0,1,0,0))
#Start with heat map
image(algae_heat$x, algae_heat$y, scaled_counts_algae, 
      xlab="", ylab="", 
      main="", col = c('white', gray_vals_algae),
      breaks = col_breaks_algae, cex.axis = 1.85,
      ylim = c(-0.05, 0.7), xlim = c(0, 100))
#Now add field data
points(plot_dat$perc_hab, plot_dat$AlgaeStd, pch = 21, cex = 1.75, bg = 'deepskyblue2')
lines(x_var, y_algae_beta, lwd = 7.5, col = 'black')
lines(x_var, y_algae_beta, lwd = 5, col = 'deepskyblue2')
lines(x_var, pred_intervals_algae[,1], lty = 2, col = 'deepskyblue2', lwd = 3)
lines(x_var, pred_intervals_algae[,2], lty = 2, col = 'deepskyblue2', lwd = 3)
mtext("Percent of reef habitat", side = 1, cex = 2.5, line = 3)
mtext("Proportional cover of macroalgae", side = 2, cex = 2.5, line = 3)

#Makes the heatmap legend which we will crop later
dev.new()
image(1, legend_vals_algae, matrix(legend_vals_algae, nrow = 1), 
      col = c('white',gray_vals_algae), xaxt = "n", xlab = "", ylab = "",
      cex.axis = 2)

################
################
#6c --> CORAL PLOT
################
################
###
#Create heat map from model outcomes
###
gray_vals_coral <- colorRampPalette(c('gray90', 'gray25'))(14)
col_breaks_coral <- c(0, seq(0,0.7, by = 0.05))
legend_vals_coral <- seq(0, 70, length.out = 15)

#Create the 2d histogram values with counts
coral_heat <- hist2d(output_sum_rep$perc_hab, output_sum_rep$end_C, show=FALSE, nbins = 5)

# Scale counts by habitat bin to correct for unequal outcomes in each bin
scaled_counts_coral <- sweep(coral_heat$counts, 1, rowSums(coral_heat$counts), "/")

####
#Prepare field data
####
#Process data for plotting
plot_dat <- dat %>%
  mutate(perc_hab = 100 - SoftSub) %>%
  left_join(n_segments) %>% 
  mutate(quadrant_struc = as.factor(ifelse(perc_hab <= 25, "1",
                                           ifelse(perc_hab > 25 & perc_hab <= 50, "2",
                                                  ifelse(perc_hab > 50 & perc_hab <= 75, "3", "4")))))
#Run beta regression
coral_beta <- betareg(CoralStd ~ perc_hab, data = plot_dat)
y_coral_beta <- betareg::predict(coral_beta, data.frame("perc_hab" = x_var))

#Bootstrap confidence intervals
boot_predict_coral <- function(data, indices) {
  boot_data <- data[indices, ]
  boot_model <- betareg(CoralStd ~ perc_hab, data = boot_data)
  predict(boot_model, newdata = data.frame("perc_hab" = x_var), type = "response")
}

# Run bootstrap
boot_results_coral <- boot(plot_dat, boot_predict_coral, R = 1000)

# Calculate percentile intervals
pred_intervals_coral <- t(apply(boot_results_coral$t, 2, quantile, probs = c(0.025, 0.975)))

###########
#Makes combined figure for coral
###########
dev.new()
par(oma=c(0,1,0,0))
#Start with heat map
image(coral_heat$x, coral_heat$y, scaled_counts_coral, 
      xlab="", ylab="", 
      main="", col = c('white', gray_vals_coral),
      breaks = col_breaks_coral, cex.axis = 1.85,
      ylim = c(-0.05, 0.8), xlim = c(0, 100))
#Now add field data
points(plot_dat$perc_hab, plot_dat$CoralStd, pch = 21, cex = 1.75, bg = 'deepskyblue2')
lines(x_var, y_coral_beta, lwd = 7.5, col = 'black')
lines(x_var, y_coral_beta, lwd = 5, col = 'deepskyblue2')
lines(x_var, pred_intervals_coral[,1], lty = 2, col = 'deepskyblue2', lwd = 3)
lines(x_var, pred_intervals_coral[,2], lty = 2, col = 'deepskyblue2', lwd = 3)
mtext("Percent of reef habitat", side = 1, cex = 2.5, line = 3)
mtext("Proportional cover of coral", side = 2, cex = 2.5, line = 3)