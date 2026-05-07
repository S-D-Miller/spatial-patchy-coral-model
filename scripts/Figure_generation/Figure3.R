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

####
#Create filtered dataframes for each herbivory value
###
scatter_dv6 <- hysteresis_outputs %>%
  filter(dV == 6)

scatter_dv4 <- hysteresis_outputs %>%
  filter(dV == 4)

scatter_dv2 <- hysteresis_outputs %>%
  filter(dV == 2)

scatter_dv05 <- hysteresis_outputs %>%
  filter(dV == 0.05)

#Sets color scheme
coral_col <- "firebrick2"
macro_col <- "deepskyblue4"

coral_circle_col <- "tan3"
macro_circle_col <- "slateblue"

##########
#Figure 5a --> dV = 6
##########
dev.new()
par(oma = c(1,1,0,0), xpd = NA)
plot(scatter_dv6$p_struc[scatter_dv6$high_initial_var == "C"],
     jitter(scatter_dv6$prop_C[scatter_dv6$high_initial_var == "C"]),
     ylim = c(0,1.15), pch = 16, col = alpha(coral_col, 0.5),
     xlab = "Proportion of habitat cells in landscape",
     ylab = "Proportion of cells dominated by coral",
     cex.axis = 1.75, cex.lab = 2)
points(scatter_dv6$p_struc[scatter_dv6$high_initial_var == "Mi"],
       jitter(scatter_dv6$prop_C[scatter_dv6$high_initial_var == "Mi"]),
       pch = 17, col = alpha(macro_col, 0.5))
#Places circle around our example landscape
points(scatter_dv6$p_struc[scatter_dv6$unique_id_rep == "0.1_0.1_10"],
       jitter(scatter_dv6$prop_C[scatter_dv6$unique_id_rep == "0.1_0.1_10"]),
       pch = 21, col = c(coral_circle_col,macro_circle_col), cex = 5, lwd = 5)
text(0.01, 1.125, labels = "a)", cex = 2.5)


##########
#Figure 5b --> dV = 4
##########
dev.new()
par(oma = c(1,1,0,0), xpd = NA)
plot(scatter_dv4$p_struc[scatter_dv4$high_initial_var == "C"],
     jitter(scatter_dv4$prop_C[scatter_dv4$high_initial_var == "C"]),
     ylim = c(0,1.15), pch = 16, col = alpha(coral_col, 0.5),
     xlab = "Proportion of habitat cells in landscape",
     ylab = "Proportion of cells dominated by coral",
     cex.axis = 1.75, cex.lab = 2)
points(scatter_dv4$p_struc[scatter_dv4$high_initial_var == "Mi"],
       jitter(scatter_dv4$prop_C[scatter_dv4$high_initial_var == "Mi"], factor = 100),
       pch = 17, col = alpha(macro_col, 0.5))
#Places circle around our example landscape
points(scatter_dv4$p_struc[scatter_dv4$unique_id_rep == "0.1_0.1_10"],
       scatter_dv4$prop_C[scatter_dv4$unique_id_rep == "0.1_0.1_10"],
       pch = 21, col = c(coral_circle_col,macro_circle_col), cex = 5, lwd = 5)
text(0.01, 1.125, labels = "b)", cex = 2.5)

##########
#Figure 5c --> dV = 2
##########
dev.new()
par(oma = c(1,1,0,0), xpd = NA)
plot(scatter_dv2$p_struc[scatter_dv2$high_initial_var == "Mi"],
     jitter(scatter_dv2$prop_C[scatter_dv2$high_initial_var == "Mi"]),
     ylim = c(0,1.15), pch = 17, col = alpha(macro_col, 0.5),
     xlab = "Proportion of habitat cells in landscape",
     ylab = "Proportion of cells dominated by coral",
     cex.axis = 1.75, cex.lab = 2)
points(scatter_dv2$p_struc[scatter_dv2$high_initial_var == "C"],
       jitter(scatter_dv2$prop_C[scatter_dv2$high_initial_var == "C"], factor = 100),
       pch = 16, col = alpha(coral_col, 0.5))
#Places circle around our example landscape
points(scatter_dv2$p_struc[scatter_dv2$unique_id_rep == "0.1_0.1_10"],
       scatter_dv2$prop_C[scatter_dv2$unique_id_rep == "0.1_0.1_10"],
       pch = 21, col = c(coral_circle_col,macro_circle_col), cex = 5, lwd = 5)
text(0.01, 1.125, labels = "c)", cex = 2.5)

##########
#Figure 5d --> dV = 0.05
##########
dev.new()
par(oma = c(1,1,0,0), xpd = NA)
plot(scatter_dv05$p_struc[scatter_dv05$high_initial_var == "C"],
     jitter(scatter_dv05$prop_C[scatter_dv05$high_initial_var == "C"]),
     ylim = c(0,1.15), pch = 16, col = alpha(coral_col, 0.5),
     xlab = "Proportion of habitat cells in landscape",
     ylab = "Proportion of cells dominated by coral",
     cex.axis = 1.75, cex.lab = 2)
points(scatter_dv05$p_struc[scatter_dv05$high_initial_var == "Mi"],
       jitter(scatter_dv05$prop_C[scatter_dv05$high_initial_var == "Mi"]),
       pch = 17, col = alpha(macro_col, 0.5))
#Places circle around our example landscape
points(scatter_dv05$p_struc[scatter_dv05$unique_id_rep == "0.1_0.1_10"],
       jitter(scatter_dv05$prop_C[scatter_dv05$unique_id_rep == "0.1_0.1_10"]),
       pch = 21, col = c(coral_circle_col,macro_circle_col), cex = 5, lwd = 5)
text(0.01, 1.125, labels = "d)", cex = 2.5)