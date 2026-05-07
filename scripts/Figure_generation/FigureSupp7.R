rm(list=ls())

##############
##############
#Loads in necessary packages
##############
##############
library(tidyverse)
library(akima)

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

scatter_4 <- hysteresis_outputs %>%
  filter(dV == 4)
scatter_2 <- hysteresis_outputs %>%
  filter(dV == 2)

##############
##############
#Data processing
##############
##############
#Creates our color scheme
heat_cols <- viridis::plasma(10)

#########
#########
#Figure S4a
#########
#########
#Filters for just when Mi was initially high
plot_s4a <- scatter_4 %>%
  filter(high_initial_var == "Mi")

#Interpolates values for the heat map and sets x/y values
#First three arguments mean we're interpolating prop_C (third argument) by values in prop_structure and clump
interp_s4a <- with(plot_s4a, interp(prop_structure,clump,prop_C, duplicate = "mean", xo=seq(0.2,0.7, by = 0.1), yo = seq(0,0.8,by=0.1), extrap = T))

#########
#########
#Figure S4b
#########
#########
#Filters for just when C was initially high
plot_s4b <- scatter_2 %>%
  filter(high_initial_var == "C")
#Interpolates values for the heat map and sets x/y values
#First three arguments mean we're interpolating prop_C (third argument) by values in prop_structure and clump
interp_s4b <- with(plot_s4b, interp(prop_structure,clump,end_C, duplicate = "mean", xo=seq(0.2,0.7, by = 0.1), yo = seq(0,0.8,by=0.1), extrap = T))

##############
##############
#Make plots
##############
##############
########
########
#Figure S4a --> dv = 4, high Mi initial
########
########
#Makes the JPEG call for saving
jpeg("output\\heat_dv4_high_macro_figS4a.jpg", 
     width = 10, height = 6, units = "in", res = 300)

#Makes the layout for saving the legend, too
layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
#Creates the plot with image
with(interp_s4a,graphics::image(x,y,z, xlab = "Proportion of habitat cells", ylab = "Degree of clumpiness", col = heat_cols,
                          main = ""), breaks = seq(0,1, by = 0.1), cex.lab = 2, cex.axis = 1.75)
#Makes our legend and plots
legend_image <- grDevices::as.raster(matrix(rev(heat_cols), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
text(x=1.5, y = seq(0,1,l=5), labels = c(0, 0.2, 0.4, 0.6, 0.8, 1))
rasterImage(legend_image, 0, 0, 1,1)
dev.off() #Turn off the plotting device

########
########
#Figure S4b --> dv = 2, high C initial
########
########
#Makes the JPEG call for saving
jpeg("output\\heat_dv2_high_coral_figS4b.jpg", 
     width = 10, height = 6, units = "in", res = 300)
#Makes the layout for saving the legend, too
layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
#Creates the plot with image
with(interp_s4b,graphics::image(x,y,z, xlab = "Proportion of habitat cells", ylab = "Degree of clumpiness", col = heat_cols,
                          main = ""), breaks = seq(0,1, by = 0.1), cex.lab = 2, cex.axis = 1.75)
#Makes our legend and plots
legend_image <- grDevices::as.raster(matrix(rev(heat_cols), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
text(x=1.5, y = seq(0,1,l=6), labels = c(0, 0.2, 0.4, 0.6, 0.8, 1))
rasterImage(legend_image, 0, 0, 1,1)
dev.off() #Turn off the plotting device