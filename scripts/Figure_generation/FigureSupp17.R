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

#Gets data on cluster assignment
dat_clusters <- readRDS("data\\cluster_classification.rds")

#########
#Processing cell data
#########
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

#Saves amount of habitat and clumpiness to merge onto the speed_outputs
to_merge2 <- cell_outputs %>%
  group_by(unique_id_rep) %>%
  reframe(p_struc = first(p_struc),
          clump = first(clump))

#########
#Processing cluster metadata
#########
cluster_metadata <- rbindlist(dat_clusters[[1]]) %>%
  mutate(unique_id_rep = names(dat_clusters)[1])

for(i in 2:length(dat_clusters)){
  
  temp <- rbindlist(dat_clusters[[i]]) %>%
    mutate(unique_id_rep = names(dat_clusters)[i])
  
  cluster_metadata <- rbind(cluster_metadata, temp)
  
}

#Adds unique cell ID to the metadata
cluster_metadata <- cluster_metadata %>%
  mutate(unique_cell_id = paste(unique_id_rep, cell_list, sep = "_"))

#Groups cells by clusters and calculates mean equilibrium values for all cells within the cluster
cluster_outputs <- cell_outputs %>%
  mutate(unique_cell_id = paste(unique_id_rep, cell, sep = "_")) %>%
  left_join(cluster_metadata %>% select(id, n_cells, unique_cell_id), by = join_by(unique_cell_id)) %>%
  rename(cluster_id = id,
         cluster_size = n_cells.y) %>%
  mutate(unique_cluster_id = paste(unique_id_rep, cluster_id, sep = "_")) %>%
  group_by(unique_cluster_id, high_initial_var, dV) %>%
  summarize(unique_id_rep = first(unique_id_rep),
            p_struc = first(p_struc),
            clump = first(clump),
            high_C_init_high_C = first(high_C_init_high_C),
            cluster_id = first(cluster_id),
            cluster_size = first(cluster_size),
            end_C = mean(end_C),
            end_Mv = mean(end_Mv),
            end_Mi = mean(end_Mi),
            end_M = mean(end_M)
  )

#Gets list of unique landscape IDs
id_list <- unique(cell_outputs$unique_id_rep) #get unique landscape IDs

#Calculates cluster level transition speeds for the first landscape
cluster_data <- calculate_transition_speed_cluster(cluster_outputs, rep_id = id_list[1])
#Then loops across others
for(i in 2:length(id_list)){
  
  cluster_data <- rbind(cluster_data, calculate_transition_speed_cluster(cluster_outputs, rep_id = id_list[i]))
  
}

cluster_data <- cluster_data %>%
  left_join(to_merge2)

############
############
#Makes cluster-level plots
############
############

#Sets color scheme
coral_col_legend <- "firebrick2"
macro_col_legend <- "deepskyblue4"
macro_col <- rgb(0,0.40784313725490196,0.5450980392156862, 0.5)
coral_col <- rgb(0.9333333333333333, 0.17254901960784313, 0.17254901960784313, 0.5)

#Makes plot for lowest herbivory with high coral based on initial condition
#red = point high coral state is lost
#blue = point low coral state recovers
dev.new()
par(oma = c(1,1,0,0))
plot(jitter(cluster_data$cluster_size, 10),
     jitter(cluster_data$C_collapse, 1.5), ylim = c(0,5.5), pch = 16, col = coral_col, cex = 1.25,
     cex.axis = 1.5, ylab = "", xlab = "")
points(jitter(cluster_data$cluster_size, 10),
       jitter(cluster_data$M_collapse, 1.5), pch = 16, col = macro_col, cex = 1.25)
text(4,5.4, "a)", cex = 2.75)
mtext("Cells in patch", side = 1, cex = 2, line = 3.5)
mtext("Lowest herbivory with high coral", side = 2, cex = 2, line = 3)
legend("bottomright", c("High C initial", "High M initial"), fill = c(coral_col_legend, macro_col_legend), cex = 1.5)

#Makes plot for highest herbivory with low coral based on initial condition
#red = point high coral state flips to low coral
#blue = point low coral state begins recovery to high coral
dev.new()
par(oma = c(1,1,0,0))
plot(jitter(cluster_data$cluster_size, 10),
     jitter(cluster_data$C_recover, 1.5), ylim = c(0,5), pch = 16, col = coral_col, cex = 1.25,
     cex.axis = 1.5, ylab = "", xlab = "")
points(jitter(cluster_data$cluster_size, 10),
       jitter(cluster_data$M_recover, 1.5), pch = 16, col = macro_col, cex = 1.25)
text(5,4.9, "b)", cex = 2.75)
mtext("Cells in patch", side = 1, cex = 2, line = 3.5)
mtext("Highest herbivory with low coral", side = 2, cex = 2, line = 3)
legend("bottomright", c("High C initial", "High M initial"), fill = c(coral_col_legend, macro_col_legend), cex = 1.5)