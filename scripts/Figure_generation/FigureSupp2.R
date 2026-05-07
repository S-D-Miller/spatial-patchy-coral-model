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
nonspatial_outputs <- read.csv("data\\nonspatial_outputs.csv")

uniform_outputs <- read.csv("data\\uniform_outputs_default_params.csv")

##############
##############
#MAKES PLOTS
##############
##############
coral_col <- "firebrick2"
macro_col <- "deepskyblue4"

###############
###############
#Figure 1a --> Nonspatial
###############
###############
dev.new()
par(oma = c(1,2,0,0), xpd = NA)
#Starts by plotting gray reference lines for unifrom landscape
plot(nonspatial_outputs$dV[nonspatial_outputs$high_initial_var == "C"],
     nonspatial_outputs$end_C[nonspatial_outputs$high_initial_var == "C"],
     type = 'l', col = coral_col, lwd = 2.5, ylim = c(0,0.5), xlab = "Herbivory (dV)",
     ylab = "Ending coral cover", main = "", cex.axis = 1.75, cex.lab = 2)
points(nonspatial_outputs$dV[nonspatial_outputs$high_initial_var == "Mi"],
       nonspatial_outputs$end_C[nonspatial_outputs$high_initial_var == "Mi"],
       type = 'l', col = macro_col, lty = 2, lwd = 2.5)
text(0.15, 0.49, labels = "a)", cex = 2.5)
legend('bottomright', col = c(coral_col, macro_col), legend = c("High C initial", "High M initial"),
       lty = c(1,2, 1, 2), lwd = 2.5, cex = 1)

###############
###############
#Figure 1b --> Uniform landscape
###############
###############
dev.new()
par(oma = c(1,2,0,0), xpd = NA)
#Starts with high initial coral
plot(uniform_outputs$dV[uniform_outputs$high_initial_var == "C"],
     uniform_outputs$end_C[uniform_outputs$high_initial_var == "C"],
     type = 'l', col = coral_col, lwd = 2.5, ylim = c(0,0.5), xlab = "Herbivory (dV)",
     ylab = "Mean ending coral cover", main = "", cex.axis = 1.75, cex.lab = 2)
#Then plots high initial macroalgae
points(uniform_outputs$dV[uniform_outputs$high_initial_var == "Mi"],
       uniform_outputs$end_C[uniform_outputs$high_initial_var == "Mi"],
       type = 'l', col = macro_col, lty = 2, lwd = 2.5)
text(0.15, 0.49, labels = "b)", cex = 2.5)
legend('bottomright', col = c(coral_col, macro_col), legend = c("High C initial", "High M initial"),
       lty = c(1,2, 1, 2), lwd = 2.5, cex = 1)