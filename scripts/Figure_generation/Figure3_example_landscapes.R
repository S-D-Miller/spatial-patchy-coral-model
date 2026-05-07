rm(list=ls())

##############
##############
#Loads in necessary packages
##############
##############
library(tidyverse)
source("scripts\\model_functions.R")

library(landscapemetrics)
library(JuliaCall)
library(data.table)
julia_setup()
julia_library("DifferentialEquations, RecursiveArrayTools, LinearAlgebra")
julia_command("using DifferentialEquations, RecursiveArrayTools, LinearAlgebra")
de <- diffeqr::diffeq_setup()

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

#############################
#############################
#RUNS MODEL
#############################
#############################

#Sets up our function to send to Julia
f <- julia_eval("function f(du,u,p,t)
    Rc, Rm, pm, Gtc, Gtv, Gti, y, Dc, Dv, Di, w, neighbor_mat_m, neighbor_mat_c, dispersal_mat, C_growth, Mi_growth, Mi_disp, turf = p
    C = @view u[1,:]
    Mv = @view u[2,:]
    Mi = @view u[3,:]
    dC = @view du[1,:]
    dMv = @view du[2,:]
    dMi = @view du[3,:]

    @. turf = 1 - (C+Mv+Mi)
    LinearAlgebra.mul!(C_growth, neighbor_mat_c, C)
    LinearAlgebra.mul!(Mi_growth, neighbor_mat_m, Mi)
    LinearAlgebra.mul!(Mi_disp, dispersal_mat, Mi)
    
    @. dC = (Rc * turf) + (Gtc * turf * C_growth) - (y * Gti * Mi_growth * C) - (Dc * C)
    @. dMv = (Rm * turf) + (pm * turf * Mi_disp) + (Gtv * turf * Mv) - (Dv * Mv) - (w * Mv)
    @. dMi = (w * Mv) + (Gti * turf * Mi_growth) + (y * Gti * C * Mi_growth) - (Di * Mi)
    
end")

###########
#Defines default parameters
###########
Rc <- 0.001
Rm <- 0.0001
pm <- 0.5
Gtc <- 0.1
Gtv <- 0.2
Gti <- 0.4
y <- 0.4
Dc <- 0.05
Dv <- 4
Di <- 0.4
w <- 2

disp_a <- 0.75
disp_b <- 3
macro_r <- 0.013
coral_r <- 0.0165

julia_assign("Rc", Rc)
julia_assign("Rm", Rm)
julia_assign("pm", pm)
julia_assign("Gtc", Gtc)
julia_assign("Gtv", Gtv)
julia_assign("Gti", Gti)
julia_assign("y", y)
julia_assign("Dc", Dc)
julia_assign("Di", Di)
julia_assign("w", w)

years_to_run <- 10000.0
cell_n <- 15
cell_res <- 0.5
save_int <- 10

tspan <- c(0.0,years_to_run)
julia_assign("tspan", tspan)

#Generates our example landscape 
set.seed(132)
input_landscape <- generate_inputs_from_random_landscape(nrow = cell_n, ncol = cell_n, resolution = cell_res, prop_cluster = 0.1,
                                                                    prop_structure = 0.1)
  


prob_func <- julia_eval("function prob_func(prob, i, repeat)
    remake(prob, p = VectorOfArray([Rc, Rm, pm, Gtc, Gtv, Gti, y, Dc, Dv_list[i], Di, w, neighbor_mat_m, neighbor_mat_c, dispersal_mat, C_growth, Mi_growth, Mi_disp, turf]))
end")


#Creates the initial grids and matrices to start the model
input_grids <- generate_grids_from_raster(input_landscape, macro_r = macro_r, coral_r = coral_r,
                                          dispersal_a = disp_a, dispersal_b = disp_b, as.torus = T,
                                          lat_rebound = T, disp_kernel = "exp_power",
                                          normalize_self = 0.5)

#Uses initial grids and initializes them based on high and low coral conditions
initial_grid_high_C <- initialize_julia_inputs(input_grids, C = 0.6, Mv = 0.1, Mi = 0.1)
initial_grid_low_C <- initialize_julia_inputs(input_grids, C = 0.1, Mv = 0.1, Mi = 0.6)

#Saves the high and low coral initial models in a list that gets accessed later
model_initializations <- list(initial_grid_high_C, initial_grid_low_C)

neighbor_mat_m <- initial_grid_high_C[[1]] #Neighbor dispersal for macroalgae
neighbor_mat_c <- initial_grid_high_C[[2]] #Neighbor dispersal for coral
dispersal_mat <- initial_grid_high_C[[3]] #Dispersal matrix for Mi -> Mv

julia_assign("neighbor_mat_m", neighbor_mat_m)
julia_assign("neighbor_mat_c", neighbor_mat_c)
julia_assign("dispersal_mat", dispersal_mat)

#####################
#####################
#SOLVES MODEL FOR EVERY EXAMPLE DV AND INITIAL CONDITION
#####################
#####################
########
#DV = 6
########
julia_assign("Dv", 6)

#Initializes with high coral
julia_assign("init_C", model_initializations[[1]][[4]][1,])
julia_assign("init_Mv", model_initializations[[1]][[4]][2,])
julia_assign("init_Mi", model_initializations[[1]][[4]][3,])

#Runs the model
julia_source("Scripts\\spatial_briggs_for_R.jl")
prob <- julia_eval("ODEProblem(f, u0, tspan, p)")
temp <- solve_julia_model(prob)
dimnames(temp)[[3]] <- which(input_grids[[4]] == 1)

#Saves output
dv6_high <- temp

#Initializes with low coral
julia_assign("init_C", model_initializations[[2]][[4]][1,])
julia_assign("init_Mv", model_initializations[[2]][[4]][2,])
julia_assign("init_Mi", model_initializations[[2]][[4]][3,])

#Runs the model
julia_source("Scripts\\spatial_briggs_for_R.jl")
prob <- julia_eval("ODEProblem(f, u0, tspan, p)")
temp <- solve_julia_model(prob)
dimnames(temp)[[3]] <- which(input_grids[[4]] == 1)

#Saves output
dv6_low <- temp

########
#DV = 4
########
julia_assign("Dv", 4)

#Initializes with high coral
julia_assign("init_C", model_initializations[[1]][[4]][1,])
julia_assign("init_Mv", model_initializations[[1]][[4]][2,])
julia_assign("init_Mi", model_initializations[[1]][[4]][3,])

#Runs the model
julia_source("Scripts\\spatial_briggs_for_R.jl")
prob <- julia_eval("ODEProblem(f, u0, tspan, p)")
temp <- solve_julia_model(prob)
dimnames(temp)[[3]] <- which(input_grids[[4]] == 1)

#Saves output
dv4_high <- temp

#Initializes with low coral
julia_assign("init_C", model_initializations[[2]][[4]][1,])
julia_assign("init_Mv", model_initializations[[2]][[4]][2,])
julia_assign("init_Mi", model_initializations[[2]][[4]][3,])

#Runs the model
julia_source("Scripts\\spatial_briggs_for_R.jl")
prob <- julia_eval("ODEProblem(f, u0, tspan, p)")
temp <- solve_julia_model(prob)
dimnames(temp)[[3]] <- which(input_grids[[4]] == 1)

#Saves output
dv4_low <- temp

########
#DV = 2
########
julia_assign("Dv", 2)

#Initializes with high coral
julia_assign("init_C", model_initializations[[1]][[4]][1,])
julia_assign("init_Mv", model_initializations[[1]][[4]][2,])
julia_assign("init_Mi", model_initializations[[1]][[4]][3,])

#Runs the model
julia_source("Scripts\\spatial_briggs_for_R.jl")
prob <- julia_eval("ODEProblem(f, u0, tspan, p)")
temp <- solve_julia_model(prob)
dimnames(temp)[[3]] <- which(input_grids[[4]] == 1)

#Saves output
dv2_high <- temp

#Initializes with low coral
julia_assign("init_C", model_initializations[[2]][[4]][1,])
julia_assign("init_Mv", model_initializations[[2]][[4]][2,])
julia_assign("init_Mi", model_initializations[[2]][[4]][3,])

#Runs the model
julia_source("Scripts\\spatial_briggs_for_R.jl")
prob <- julia_eval("ODEProblem(f, u0, tspan, p)")
temp <- solve_julia_model(prob)
dimnames(temp)[[3]] <- which(input_grids[[4]] == 1)

#Saves output
dv2_low <- temp

########
#DV = 0.05
########
julia_assign("Dv", 0.05)

#Initializes with high coral
julia_assign("init_C", model_initializations[[1]][[4]][1,])
julia_assign("init_Mv", model_initializations[[1]][[4]][2,])
julia_assign("init_Mi", model_initializations[[1]][[4]][3,])

#Runs the model
julia_source("Scripts\\spatial_briggs_for_R.jl")
prob <- julia_eval("ODEProblem(f, u0, tspan, p)")
temp <- solve_julia_model(prob)
dimnames(temp)[[3]] <- which(input_grids[[4]] == 1)

#Saves output
dv05_high <- temp

#Initializes with low coral
julia_assign("init_C", model_initializations[[2]][[4]][1,])
julia_assign("init_Mv", model_initializations[[2]][[4]][2,])
julia_assign("init_Mi", model_initializations[[2]][[4]][3,])

#Runs the model
julia_source("Scripts\\spatial_briggs_for_R.jl")
prob <- julia_eval("ODEProblem(f, u0, tspan, p)")
temp <- solve_julia_model(prob)
dimnames(temp)[[3]] <- which(input_grids[[4]] == 1)

#Saves output
dv05_low <- temp

#######################
#######################
#######################
#MAKES PLOTS
#######################
#######################
#######################

make_landscape_plot <- function(temp){
  #Function that takes cell-level model output for this landscape and sets the values on the grid to those values
  new_vals <- data.frame("val" = input_grids[[5]])
  new_vals$val[new_vals$val == 1] <- temp[nrow(temp),'C',]
  
  col_frame <- c("white", plasma(12))
  
  plot_landscape <- input_landscape %>%
    setValues(new_vals)
  
  #Finally, makes the plot
  br <- c(0, 0.0001, seq(0.05, 0.55, by = 0.05))
  dev.new()
  terra::plot(plot_landscape, main = "", col = col_frame, legend = F, pax = list(lab = "", tick = F), breaks = br)
  
}

###########
#Makes all the plots
###########
make_landscape_plot(dv6_high)
make_landscape_plot(dv6_low)
make_landscape_plot(dv4_high)
make_landscape_plot(dv4_low)
make_landscape_plot(dv2_high)
make_landscape_plot(dv2_low)
make_landscape_plot(dv05_high)
make_landscape_plot(dv05_low)