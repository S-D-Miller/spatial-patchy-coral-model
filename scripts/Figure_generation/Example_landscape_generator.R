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
#Sets up size and resolution for the landscapes
##############
##############
nrow = 15
ncol = 15
cell.res = 0.5

###########
#Makes plots for hysteresis diagrams (Figure 2)
###########
#Set our random seed
set.seed(42)

#Generate example landscape for 0-25% habitat
landscape0_25 <- generate_inputs_from_random_landscape_no_cluster(nrow = nrow, ncol = ncol, resolution = cell.res,
                                                                  prop_structure = 0.125)
#Plot it
dev.new()
plot(landscape0_25, col = c("white", "black"), legend = NA, pax = list(lab = "", tick = F))

#Repeat for 25-50% habitat
set.seed(42)
landscape25_50 <- generate_inputs_from_random_landscape_no_cluster(nrow = nrow, ncol = ncol, resolution = cell.res,
                                                                   prop_structure = 0.35)
dev.new()
plot(landscape25_50, col = c("white", "black"), legend = NA, pax = list(lab = "", tick = F))

#Repeat for 50-75% habitat
set.seed(42)
landscape50_75 <- generate_inputs_from_random_landscape_no_cluster(nrow = nrow, ncol = ncol, resolution = cell.res,
                                                                   prop_structure = 0.625)
dev.new()
plot(landscape50_75, col = c("white", "black"), legend = NA, pax = list(lab = "", tick = F))

#Repeat for 75-100% habitat
set.seed(42)
landscape75_100 <- generate_inputs_from_random_landscape_no_cluster(nrow = nrow, ncol = ncol, resolution = cell.res,
                                                                    prop_structure = 0.875)
dev.new()
plot(landscape75_100, col = c("white", "black"), legend = NA, pax = list(lab = "", tick = F))


#################
#Makes plots for supplement example Figure S4
#################
#Same as before...set random seed and generate
#Start with no cluster
set.seed(42)
landscape1 <- generate_inputs_from_random_landscape_no_cluster(nrow = nrow, ncol = ncol, resolution = cell.res,
                                                               prop_structure = 0.1)
set.seed(42)
landscape2 <- generate_inputs_from_random_landscape_no_cluster(nrow = nrow, ncol = ncol, resolution = cell.res,
                                                               prop_structure = 0.4)
set.seed(42)
landscape3 <- generate_inputs_from_random_landscape_no_cluster(nrow = nrow, ncol = ncol, resolution = cell.res,
                                                               prop_structure = 0.7)

#Mild clustering
set.seed(100)
landscape4 <- generate_inputs_from_random_landscape(nrow = nrow, ncol = ncol, resolution = cell.res,
                                                    prop_structure = 0.1, prop_cluster = 0.25)
set.seed(42)
landscape5 <- generate_inputs_from_random_landscape(nrow = nrow, ncol = ncol, resolution = cell.res,
                                                    prop_structure = 0.3, prop_cluster = 0.25)
set.seed(100)
landscape6 <- generate_inputs_from_random_landscape(nrow = nrow, ncol = ncol, resolution = cell.res,
                                                    prop_structure = 0.6, prop_cluster = 0.25)

#High clustering
set.seed(42)
landscape7 <- generate_inputs_from_random_landscape(nrow = nrow, ncol = ncol, resolution = cell.res,
                                                    prop_structure = 0.1, prop_cluster = 0.5)
set.seed(42)
landscape8 <- generate_inputs_from_random_landscape(nrow = nrow, ncol = ncol, resolution = cell.res,
                                                    prop_structure = 0.3, prop_cluster = 0.5)
set.seed(42)
landscape9 <- generate_inputs_from_random_landscape(nrow = nrow, ncol = ncol, resolution = cell.res,
                                                    prop_structure = 0.7, prop_cluster = 0.4)

#Plot them all to save and combine later
dev.new()
plot(landscape1, col = c("white", "black"), legend = NA, pax = list(lab = "", tick = F))
dev.new()
plot(landscape2, col = c("white", "black"), legend = NA, pax = list(lab = "", tick = F))
dev.new()
plot(landscape3, col = c("white", "black"), legend = NA, pax = list(lab = "", tick = F))
dev.new()
plot(landscape4, col = c("white", "black"), legend = NA, pax = list(lab = "", tick = F))
dev.new()
plot(landscape5, col = c("white", "black"), legend = NA, pax = list(lab = "", tick = F))
dev.new()
plot(landscape6, col = c("white", "black"), legend = NA, pax = list(lab = "", tick = F))
dev.new()
plot(landscape7, col = c("white", "black"), legend = NA, pax = list(lab = "", tick = F))
dev.new()
plot(landscape8, col = c("white", "black"), legend = NA, pax = list(lab = "", tick = F))
dev.new()
plot(landscape9, col = c("white", "black"), legend = NA, pax = list(lab = "", tick = F))
