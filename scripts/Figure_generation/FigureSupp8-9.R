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
#DEFAULT VALUES
hysteresis_outputs <- read.csv("data\\nocluster_outputs_default_params.csv")
hysteresis_outputs$unique_id_rep <- paste(hysteresis_outputs$prop_cluster, hysteresis_outputs$prop_structure, hysteresis_outputs$rep, sep = "_")
hysteresis_outputs_cluster <- read.csv("data\\cluster_outputs_default_params.csv")
hysteresis_outputs_cluster$unique_id_rep <- paste(hysteresis_outputs_cluster$prop_cluster, hysteresis_outputs_cluster$prop_structure, hysteresis_outputs_cluster$rep, sep = "_")
hysteresis_outputs <- rbind(hysteresis_outputs, hysteresis_outputs_cluster)

#DOUBLE a 
hysteresis_outputs_a15 <- read.csv("data\\nocluster_outputs_a15.csv")
hysteresis_outputs_a15$unique_id_rep <- paste(hysteresis_outputs_a15$prop_cluster, hysteresis_outputs_a15$prop_structure, hysteresis_outputs_a15$rep, sep = "_")
hysteresis_outputs_cluster_a15 <- read.csv("data\\cluster_outputs_a15.csv")
hysteresis_outputs_cluster_a15$unique_id_rep <- paste(hysteresis_outputs_cluster_a15$prop_cluster, hysteresis_outputs_cluster_a15$prop_structure, hysteresis_outputs_cluster_a15$rep, sep = "_")
hysteresis_outputs_a15 <- rbind(hysteresis_outputs_a15, hysteresis_outputs_cluster_a15)

#Uniform landscape
uniform_hysteresis <- read.csv("data\\uniform_outputs_default_params.csv")

##############
##############
#CALCULATE EQUILIBRIUM CHARACTERISTICS
##############
##############
#######
#######
#Default values
equilibrium_characteristics_default <- hysteresis_outputs %>%
  group_by(unique_id_rep) %>%
  reframe(p_struc = mean(p_struc),
          clump = mean(clump),
          high_C = mean(high_C_init_high_C),
          low_C = mean(low_C_init_low_C),
          high_M = mean(high_M_init_low_C),
          low_M = mean(low_M_init_high_C))
#######
#######
#Double dispersal --> a = 1.5
equilibrium_characteristics_a15 <- hysteresis_outputs_a15 %>%
  group_by(unique_id_rep) %>%
  reframe(p_struc = mean(p_struc),
          clump = mean(clump),
          high_C = mean(high_C_init_high_C),
          low_C = mean(low_C_init_low_C),
          high_M = mean(high_M_init_low_C),
          low_M = mean(low_M_init_high_C))

##############
##############
#MAKE PLOTS
##############
##############
########
########
#FIGURE S5 --> DEFAULT DISPERSAL
###
#LOW HERBIVORY
#Macroalgae equilibria
dev.new()
par(oma = c(0,1,0,0), xpd = NA)
plot(equilibrium_characteristics_default$p_struc, equilibrium_characteristics_default$high_M,
     pch = 16, xlab = "Amount of reef habitat", ylab = "Ending macroalgae cover under low herbivory",
     cex.lab = 1.75, cex.axis = 1.35)
text(0.01, 0.559, "a)", cex = 2.25)

#Coral equilibria
dev.new()
par(oma = c(0,1,0,0), xpd = NA)
plot(equilibrium_characteristics_default$p_struc, equilibrium_characteristics_default$low_C, pch = 16,
     xlab = "Amount of reef habitat", ylab = "Ending coral cover under low herbivory",
     cex.lab = 1.75, cex.axis = 1.35, ylim = c(0.005, 0.011))
text(0.01, 0.011, "b)", cex = 2.25)

###
#HIGH HERBIVORY
#Macroalgae equilibria
dev.new()
par(oma = c(0,1,0,0), xpd = NA)
plot(equilibrium_characteristics_default$p_struc, equilibrium_characteristics_default$low_M,
     pch = 16, xlab = "Amount of reef habitat", ylab = "Ending macroalgae cover under high herbivory",
     cex.lab = 1.75, cex.axis = 1.35, ylim = c(0,0.0001))
text(0.01, 0.0001, "c)", cex = 2.25)

#Coral equilibria
dev.new()
par(oma = c(0,1,0,0), xpd = NA)
plot(equilibrium_characteristics_default$p_struc, equilibrium_characteristics_default$high_C,
     pch = 16, xlab = "Amount of reef habitat", ylab = "Ending coral cover under high herbivory",
     cex.lab = 1.75, cex.axis = 1.35, ylim = c(0.509, 0.51))
text(0.01, 0.51, "d)", cex = 2.25)

########
########
#FIGURE S6 --> DOUBLE a DISPERSAL
###
#LOW HERBIVORY
#Macroalgae equilibria
dev.new()
par(oma = c(0,1,0,0), xpd = NA)
plot(equilibrium_characteristics_a15$p_struc, equilibrium_characteristics_a15$high_M,
     pch = 16, xlab = "Amount of reef habitat", ylab = "Ending macroalgae cover under low herbivory",
     cex.lab = 1.75, cex.axis = 1.35)
text(0.01, 0.559, "a)", cex = 2.25)

#Coral equilibria
dev.new()
par(oma = c(0,1,0,0), xpd = NA)
plot(equilibrium_characteristics_a15$p_struc, equilibrium_characteristics_a15$low_C, pch = 16,
     xlab = "Amount of reef habitat", ylab = "Ending coral cover under low herbivory",
     cex.lab = 1.75, cex.axis = 1.35, ylim = c(0.005, 0.013))
text(0.01, 0.013, "b)", cex = 2.25)

###
#HIGH HERBIVORY
#Macroalgae equilibria
dev.new()
par(oma = c(0,1,0,0), xpd = NA)
plot(equilibrium_characteristics_a15$p_struc, equilibrium_characteristics_a15$low_M,
     pch = 16, xlab = "Amount of reef habitat", ylab = "Ending macroalgae cover under high herbivory",
     cex.lab = 1.75, cex.axis = 1.35, ylim = c(0,0.0001))
text(0.01, 0.0001, "c)", cex = 2.25)

#Coral equilibria
dev.new()
par(oma = c(0,1,0,0), xpd = NA)
plot(equilibrium_characteristics_a15$p_struc, equilibrium_characteristics_a15$high_C, pch = 16,
     xlab = "Amount of reef habitat", ylab = "Ending coral cover under high herbivory",
     cex.lab = 1.75, cex.axis = 1.35, ylim = c(0.509, 0.510))
text(0.01, 0.51, "d)", cex = 2.25)