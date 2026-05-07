rm(list=ls())

##############
##############
#Loads in necessary packages
##############
##############
library(tidyverse)
library(data.table)
source("scripts\\model_functions.R")

##############
##############
#Loads in data
##############
##############

########
#Landscape level
########
hysteresis_outputs <- read.csv("data\\nocluster_outputs_default_params.csv")
hysteresis_outputs$unique_id_rep <- paste(hysteresis_outputs$prop_cluster, hysteresis_outputs$prop_structure, hysteresis_outputs$rep, sep = "_")
hysteresis_outputs_cluster <- read.csv("data\\cluster_outputs_default_params.csv")
hysteresis_outputs_cluster$unique_id_rep <- paste(hysteresis_outputs_cluster$prop_cluster, hysteresis_outputs_cluster$prop_structure, hysteresis_outputs_cluster$rep, sep = "_")
hysteresis_outputs <- rbind(hysteresis_outputs, hysteresis_outputs_cluster)

#########
#Cell level
#########
#Gets the high coral equilibrium for each landscape
to_merge <- hysteresis_outputs %>%
  group_by(unique_id_rep) %>%
  reframe(high_C_init_high_C = first(high_C_init_high_C))

###Loads in data
dat_cells <- readRDS("data\\cell_outputs_cluster.rds")
dat_cells2 <- readRDS("data\\cell_outputs_nocluster.rds")

#Converts list to a dataframe
#First element
cell_outputs <- rbindlist(dat_cells[[1]]) %>%
  mutate(unique_id_rep = paste(prop_cluster, prop_structure, rep, sep = "_"))
#Now loop through rest
for(i in 2:length(dat_cells)){
  
  temp <- rbindlist(dat_cells[[i]]) %>%
    mutate(unique_id_rep = paste(prop_cluster, prop_structure, rep, sep = "_"))
  
  cell_outputs <- rbind(cell_outputs, temp)
  
}

#Repeat for the non-clustered outputs
cell_outputs2 <- rbindlist(dat_cells2[[1]]) %>%
  mutate(unique_id_rep = paste(prop_cluster, prop_structure, rep, sep = "_"))

for(i in 2:length(dat_cells2)){
  
  temp <- rbindlist(dat_cells2[[i]]) %>%
    mutate(unique_id_rep = paste(prop_cluster, prop_structure, rep, sep = "_"))
  
  cell_outputs2 <- rbind(cell_outputs2, temp)
  
}

#Combine cluster and non-cluster into one dataframe and add equilibrium values of landscape
cell_outputs <- rbind(cell_outputs, cell_outputs2) %>%
  left_join(to_merge)

################
################
#Data processing
################
################
id_list <- unique(cell_outputs$unique_id_rep) #get unique landscape IDs

#Calculates transition speed for all the cells
#Just one landscape
speed_outputs <- calculate_transition_speed_cell(cell_outputs, rep_id = id_list[1])

#Then loop for the rest
for(i in 2:length(id_list)){
  
  speed_outputs <- rbind(speed_outputs, calculate_transition_speed_cell(cell_outputs, rep_id = id_list[i]))
  
}

#Saves amount of habitat and clumpiness to merge onto the speed_outputs
to_merge2 <- cell_outputs %>%
  group_by(unique_id_rep) %>%
  reframe(p_struc = first(p_struc),
          clump = first(clump))

cell_data <- speed_outputs %>%
  left_join(to_merge2) 

#Calculate transition speed for all the landscapes
#Start with one
landscape_outputs <- calculate_transition_speed_landscape(hysteresis_outputs, rep_id = id_list[1])

#Then loop the rest
for(i in 2:length(id_list)){
  
  landscape_outputs <- rbind(landscape_outputs, calculate_transition_speed_landscape(hysteresis_outputs, rep_id = id_list[i]))
  
}

landscape_data <- landscape_outputs %>%
  left_join(to_merge2)

#Finds overall max and min values for the range across scales
max_C_range <- max(c(cell_data$C_range, landscape_data$C_range)) #1.75
min_C_range <- min(c(cell_data$C_range, landscape_data$C_range)) #0.25
max_M_range <- max(c(cell_data$M_range, landscape_data$M_range)) #2
min_M_range <- min(c(cell_data$M_range, landscape_data$M_range)) #0.25

min_range <- min(min_C_range, min_M_range) #0.25
max_range <- max(max_C_range, max_M_range) #2

#Calculates abruptness for each state variable
#Normalizes so a value with min_range is 1 and a value with max_range is 0
#Cells
cell_data <- cell_data %>%
  mutate(C_abruptness = abs(1 - ((C_range - min_range) / (max_range - min_range))),
         M_abruptness = abs(1 - ((M_range - min_range) / (max_range - min_range))))
#Landscapes
landscape_data <- landscape_data %>%
  mutate(C_abruptness = abs(1 - ((C_range - min_range) / (max_range - min_range))),
         M_abruptness = abs(1 - ((M_range - min_range) / (max_range - min_range))))

#Filters data by amount of habitat for the next plots
#Cells
cells_25 <- cell_data %>%
  filter(p_struc <= .25)
cells_50 <- cell_data %>%
  filter(p_struc <= .5 & p_struc > 0.25)
cells_75 <- cell_data %>%
  filter(p_struc <= .75 & p_struc > 0.5)
cells_100 <- cell_data %>%
  filter(p_struc <= 1 & p_struc > 0.75)

#Landscapes
landscapes_25 <- landscape_data %>%
  filter(p_struc <= .25)
landscapes_50 <- landscape_data %>%
  filter(p_struc <= .5 & p_struc > 0.25)
landscapes_75 <- landscape_data %>%
  filter(p_struc <= .75 & p_struc > 0.5)
landscapes_100 <- landscape_data %>%
  filter(p_struc <= 1 & p_struc > 0.75)

#Calculates breaks and colors
macro_col <- rgb(0,0.40784313725490196,0.5450980392156862, 0.5)
coral_col <- rgb(0.9333333333333333, 0.17254901960784313, 0.17254901960784313, 0.5)
breaks = seq(0,1, by = 0.2)

#Makes plots
###########
#LANDSCAPES
###########
#Overlapping histogram for landscapes with 75-100% habitat
dev.new()
hist(landscapes_100$C_abruptness, freq = T, breaks = breaks, col = coral_col,
     xlab = "", ylab = "", cex.axis = 1.5, main = "")
hist(landscapes_100$M_abruptness, freq = T, breaks = breaks, col = macro_col, add = T)
text(0.01, 80, "a)", cex = 2.75)
legend(0.0, 25, c("High C Initial", "High M Initial"), fill = c(coral_col, macro_col), cex = 2)


#Overlapping histogram for landscapes with 50-75% habitat
dev.new()
par(xpd = NA)
hist(landscapes_75$C_abruptness, freq = T, breaks = breaks, col = coral_col,
     xlab = "", ylab = "", cex.axis = 1.5, main = "")
hist(landscapes_75$M_abruptness, freq = T, breaks = breaks, col = macro_col, add = T)
text(0.01, 140, "b)", cex = 2.75)

#Overlapping histogram for landscapes with 25-50% habitat
dev.new()
hist(landscapes_50$C_abruptness, freq = T, breaks = breaks, col = coral_col,
     xlab = "", ylab = "", cex.axis = 1.5, main = "", ylim = c(0, 70))
hist(landscapes_50$M_abruptness, freq = T, breaks = breaks, col = macro_col, add = T)
text(0.01, 70, "c)", cex = 2.75)

#Overlapping histogram for landscapes with 0-25% habitat
dev.new()
hist(landscapes_25$C_abruptness, freq = T, breaks = breaks, col = coral_col,
     xlab = "", ylab = "", cex.axis = 1.5, main = "", ylim = c(0,60))
hist(landscapes_25$M_abruptness, freq = T, breaks = breaks, col = macro_col, add = T)
text(0.01, 60, "d)", cex = 2.75)


###########
#CELLS
###########
#Overlapping histogram for cells in landscapes with 75-100% habitat
dev.new()
hist(cells_100$C_abruptness, freq = T, breaks = breaks, col = coral_col,
     xlab = "", ylab = "", cex.axis = 1.5, main = "")
hist(cells_100$M_abruptness, freq = T, breaks = breaks, col = macro_col, add = T)
text(0.01, 15000, "e)", cex = 2.75)

#Overlapping histogram for cells in landscapes with 50-75% habitat
dev.new()
hist(cells_75$C_abruptness, freq = T, breaks = breaks, col = coral_col,
     xlab = "", ylab = "", cex.axis = 1.5, main = "", ylim = c(0,20000))
hist(cells_75$M_abruptness, freq = T, breaks = breaks, col = macro_col, add = T)
text(0.01, 20000, "f)", cex = 2.75)

#Overlapping histogram for cells in landscapes with 25-50% habitat
dev.new()
hist(cells_50$C_abruptness, freq = T, breaks = breaks, col = coral_col,
     xlab = "", ylab = "", cex.axis = 1.5, main = "")
hist(cells_50$M_abruptness, freq = T, breaks = breaks, col = macro_col, add = T)
text(0.01, 10000, "g)", cex = 2.75)

#Overlapping histogram for cells in landscapes with 0-25% habitat
dev.new()
hist(cells_25$C_abruptness, freq = T, breaks = breaks, col = coral_col,
     xlab = "", ylab = "", cex.axis = 1.5, main = "", ylim = c(0,4000))
hist(cells_25$M_abruptness, freq = T, breaks = breaks, col = macro_col, add = T)
text(0.01, 4000, "h)", cex = 2.75)