rm(list=ls())

library(tidyverse)
library(landscapemetrics)
library(JuliaCall)
library(data.table)

source("scripts\\model_functions.R")

##########
#SETS UP JULIA
##########
julia_setup()
julia_library("DifferentialEquations, RecursiveArrayTools, LinearAlgebra")
julia_command("using DifferentialEquations, RecursiveArrayTools, LinearAlgebra")
de <- diffeqr::diffeq_setup()

sensitivity_reps <- read.csv("data\\sensitivity_landscape_metadata.csv")

#Creates our differential equation function that will be solved
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

#Sets parameters
Rc <- 0.001
Rm <- 0.0001
pm <- 0.5
Gtc <- 0.1
Gtv <- 0.2
Gti <- 0.4
y <- 0.4
Dc <- 0.05
Di <- 0.4
w <- 2

#Assigns parameters to Julia
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

tspan <- c(0.0,years_to_run)
save_int <- 10
julia_assign("tspan", tspan)

###################
###################
#CORAL GROWTH
###################
###################
#Creates range of values to explore sensitivity
#The default of each spatial parameter up to +- 20%
rc_list <- seq(0.0165*.8, 0.0165*1.2, by = 0.0165*.05) #default: 0.0165 +-20%
disp_a <- 0.75
disp_b <- 3
macro_r <- 0.013

#Now loops through each value, running the analysis
for(rc in rc_list){
  
  #Sets the value we're testing
  coral_r <- rc
  print(paste("Progress ---> coral_r = ", coral_r, " of ", max(rc_list), sep = ""))
  
  #Loops through the different landscapes
  for(rep in unique(sensitivity_reps$unique_id_rep)){
    
    #Sets the generation values to that of this landscape
    temp_landscape <- sensitivity_reps %>%
      filter(unique_id_rep == rep)
    
    set.seed(temp_landscape$rand.seed)
    
    if(temp_landscape$prop_cluster == 0){
      
      input_landscape <- generate_inputs_from_random_landscape_no_cluster(nrow = sqrt(temp_landscape$size), 
                                                                          ncol = sqrt(temp_landscape$size), 
                                                                          resolution = temp_landscape$res,
                                                                          prop_structure = temp_landscape$prop_structure)
      
    } else{
      
      input_landscape <- generate_inputs_from_random_landscape(nrow = sqrt(temp_landscape$size), 
                                                               ncol = sqrt(temp_landscape$size), 
                                                               resolution = temp_landscape$res,
                                                               prop_cluster = temp_landscape$prop_cluster, 
                                                               prop_structure = temp_landscape$prop_structure)
      
    }
    
    #Sets our herbivory values
    Dv_list <- c(0.05, 0.1, 0.15, 0.2, seq(0.25, 8, by = 0.25))
    julia_assign("Dv_list", Dv_list)
    #Sets our function in Julia with all our parameters
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
    
    ###############
    #Loads in dispersal matrices for this landscape and saves to Julia
    ###############
    neighbor_mat_m <- initial_grid_high_C[[1]] #Neighbor dispersal for macroalgae
    neighbor_mat_c <- initial_grid_high_C[[2]] #Neighbor dispersal for coral
    dispersal_mat <- initial_grid_high_C[[3]] #Dispersal matrix for Mi -> Mv
    
    julia_assign("neighbor_mat_m", neighbor_mat_m)
    julia_assign("neighbor_mat_c", neighbor_mat_c)
    julia_assign("dispersal_mat", dispersal_mat)
    
    #####################
    #RUNS THE MODEL ON VERY HIGH AND LOW HERBIVORY TO LEARN DIFFERENT EQUILIBRIA
    #####################
    start.time <- Sys.time()
    
    #Sets dV high
    julia_assign("Dv", 12)
    
    #Initializes with high coral
    julia_assign("init_C", model_initializations[[1]][[4]][1,])
    julia_assign("init_Mv", model_initializations[[1]][[4]][2,])
    julia_assign("init_Mi", model_initializations[[1]][[4]][3,])
    
    #Runs the model
    julia_source("scripts\\spatial_briggs_for_R.jl")
    prob <- julia_eval("ODEProblem(f, u0, tspan, p)")
    temp <- solve_julia_model(prob)
    dimnames(temp)[[3]] <- which(input_grids[[4]] == 1)
    #Saves output
    high_c_init_high_c <- mean(temp[nrow(temp),'C',])
    low_m_init_high_c <- mean(temp[nrow(temp),'Mi',] + temp[nrow(temp),'Mv',])
    
    #Initializes with low coral
    julia_assign("init_C", model_initializations[[1]][[4]][1,])
    julia_assign("init_Mv", model_initializations[[1]][[4]][2,])
    julia_assign("init_Mi", model_initializations[[1]][[4]][3,])
    
    #Runs the model
    julia_source("scripts\\spatial_briggs_for_R.jl")
    prob <- julia_eval("ODEProblem(f, u0, tspan, p)")
    temp <- solve_julia_model(prob)
    dimnames(temp)[[3]] <- which(input_grids[[4]] == 1)
    #Saves output
    high_c_init_low_c <- mean(temp[nrow(temp),'C',])
    low_m_init_low_c <- mean(temp[nrow(temp),'Mi',] + temp[nrow(temp),'Mv',])
    
    #Runs the model on this landscape at very low herbivory under both starting conditions
    julia_assign("Dv", 0.01)
    #Initializes with high coral
    julia_assign("init_C", model_initializations[[1]][[4]][1,])
    julia_assign("init_Mv", model_initializations[[1]][[4]][2,])
    julia_assign("init_Mi", model_initializations[[1]][[4]][3,])
    
    #Runs the model
    julia_source("scripts\\spatial_briggs_for_R.jl")
    prob <- julia_eval("ODEProblem(f, u0, tspan, p)")
    temp <- solve_julia_model(prob)
    dimnames(temp)[[3]] <- which(input_grids[[4]] == 1)
    #Saves output
    low_c_init_high_c <- mean(temp[nrow(temp),'C',])
    high_m_init_high_c <- mean(temp[nrow(temp),'Mi',] + temp[nrow(temp),'Mv',])
    
    #Initializes with low coral
    julia_assign("init_C", model_initializations[[1]][[4]][1,])
    julia_assign("init_Mv", model_initializations[[1]][[4]][2,])
    julia_assign("init_Mi", model_initializations[[1]][[4]][3,])
    
    #Runs the model
    julia_source("scripts\\spatial_briggs_for_R.jl")
    prob <- julia_eval("ODEProblem(f, u0, tspan, p)")
    temp <- solve_julia_model(prob)
    dimnames(temp)[[3]] <- which(input_grids[[4]] == 1)
    #Saves output
    low_c_init_low_c <- mean(temp[nrow(temp),'C',])
    high_m_init_low_c <- mean(temp[nrow(temp),'Mi',] + temp[nrow(temp),'Mv',])
    
    output_list <- list()
    
    #Loops through the two initial conditions
    for(j in 1:2){
      hysteresis_outputs <- data.frame("unique_id_rep" = character(),
                                       "alpha" = numeric(),
                                       "beta" = numeric(),
                                       "coral_r" = numeric(),
                                       "macro_r" = numeric(),
                                       "high_initial_var" = character(),
                                       "dV" = numeric(),
                                       "years_ran" = numeric(),
                                       "end_C" = numeric(),
                                       "sd_C" = numeric(),
                                       "high_C_init_high_C" = numeric(),
                                       "low_C_init_high_C" = numeric(),
                                       "high_C_init_low_C" = numeric(),
                                       "low_C_init_low_C" = numeric(),
                                       "prop_C" = numeric(),
                                       "end_Mi" = numeric(),
                                       "sd_Mi" = numeric(),
                                       "end_Mv" = numeric(),
                                       "sd_Mv" = numeric(),
                                       "high_M_init_high_C" = numeric(),
                                       "low_M_init_high_C" = numeric(),
                                       "high_M_init_low_C" = numeric(),
                                       "low_M_init_low_C" = numeric(),
                                       "prop_M" = numeric(),
                                       "prop_ratio" = numeric(),
                                       "delta_C" = numeric(),
                                       "delta_Mi" = numeric())
      
      #j goes through the high and low initial conditions and runs the model
      julia_assign("init_C", model_initializations[[j]][[4]][1,])
      julia_assign("init_Mv", model_initializations[[j]][[4]][2,])
      julia_assign("init_Mi", model_initializations[[j]][[4]][3,])
      
      #Runs the model
      julia_source("scripts\\spatial_briggs_for_R.jl")
      prob <- julia_eval("ODEProblem(f, u0, tspan, p)")
      ensembleprob <- de$EnsembleProblem(prob, prob_func = prob_func)
      
      #Solves in Julia and brings outputs back into R
      sol <- de$solve(ensembleprob,de$Tsit5(),de$EnsembleSerial(),trajectories=length(Dv_list))
      solves <- sapply(sol$u,identity)
      
      #Goes through each solve (herbivory value) and extracts data
      for(i in 1:length(Dv_list)){
        
        mat <- solves[[i]]
        mat <- sapply(mat$u,identity)
        
        model_outputs <- array(dim = c(ncol(mat), 4, dim(mat)[1]/3))
        colnames(model_outputs) <- c("C", "Mv", "Mi", "t")
        model_outputs[,"C",] <- t(mat[seq(1, nrow(mat), by = 3),])
        model_outputs[,"Mv",] <- t(mat[seq(2, nrow(mat), by = 3),])
        model_outputs[,"Mi",] <- t(mat[seq(3, nrow(mat), by = 3),])
        model_outputs[,"t",] <- (1 - (model_outputs[,'C',] + model_outputs[,'Mv',] + model_outputs[,'Mi',]))
        
        temp <- model_outputs
        dimnames(temp)[[3]] <- which(input_grids[[4]] == 1)
        
        temp_frame <- data.frame("unique_id_rep" = temp_landscape$unique_id_rep,
                                 "alpha" = disp_a,
                                 "beta" = disp_b,
                                 "coral_r" = coral_r,
                                 "macro_r" = macro_r,
                                 "high_initial_var" = ifelse(j == 1, "C", "Mi"),
                                 "dV" = Dv_list[i],
                                 "years_ran" = years_to_run,
                                 "end_C" = ifelse(dim(temp)[3] == 1, temp[nrow(temp), "C", 1], rowMeans(temp[,"C",], na.rm=T)[length(rowMeans(temp[,"C",]))]),
                                 "sd_C" = ifelse(dim(temp)[3] == 1, 0, apply(temp[,"C",], 1, FUN = sd, na.rm=T)[length(apply(temp[,"C",], 1, FUN = sd))]),
                                 "high_C_init_high_C" = high_c_init_high_c,
                                 "low_C_init_high_C" = low_c_init_high_c,
                                 "high_C_init_low_C" = high_c_init_low_c,
                                 "low_C_init_low_C" = low_c_init_low_c,
                                 "prop_C" = ifelse(j == 1, sum(temp[nrow(temp),'C',] > (0.85 * high_c_init_high_c)) / dim(temp)[3],
                                                   sum(temp[nrow(temp),'C',] > (0.85 * high_c_init_low_c)) / dim(temp)[3]),
                                 "end_Mi" = ifelse(dim(temp)[3] == 1, temp[nrow(temp), "Mi", 1], rowMeans(temp[,"Mi",], na.rm=T)[length(rowMeans(temp[,"Mi",]))]),
                                 "sd_Mi" = ifelse(dim(temp)[3] == 1, 0, apply(temp[,"Mi",], 1, FUN = sd, na.rm=T)[length(apply(temp[,"Mi",], 1, FUN = sd))]),
                                 "end_Mv" = ifelse(dim(temp)[3] == 1, temp[nrow(temp), "Mv", 1], rowMeans(temp[,"Mv",], na.rm=T)[length(rowMeans(temp[,"Mv",]))]),
                                 "sd_Mv" = ifelse(dim(temp)[3] == 1, 0, apply(temp[,"Mv",], 1, FUN = sd, na.rm=T)[length(apply(temp[,"Mv",], 1, FUN = sd))]),
                                 "high_M_init_high_C" = high_m_init_high_c,
                                 "low_M_init_high_C" = low_m_init_high_c,
                                 "high_M_init_low_C" = high_m_init_low_c,
                                 "low_M_init_low_C" = low_m_init_low_c,
                                 "prop_M" = ifelse(j == 1, sum((temp[nrow(temp),'Mi',] + temp[nrow(temp), 'Mv',]) > (0.85 * high_m_init_high_c)) / dim(temp)[3],
                                                   sum((temp[nrow(temp),'Mi',] + temp[nrow(temp), 'Mv',]) > (0.85 * high_m_init_low_c)) / dim(temp)[3]),
                                 "prop_ratio_C" =  sum(temp[nrow(temp),'C',] > (temp[nrow(temp),'Mi',] + temp[nrow(temp), 'Mv',])) / dim(temp)[3],
                                 "prop_ratio_M" =  sum(temp[nrow(temp),'C',] < (temp[nrow(temp),'Mi',] + temp[nrow(temp), 'Mv',])) / dim(temp)[3],
                                 "delta_C" = mean(temp[dim(temp)[1],"C",], na.rm=T) - mean(temp[dim(temp)[1] - 1,"C",], na.rm=T),
                                 "delta_Mi" = mean(temp[dim(temp)[1],"Mi",], na.rm=T) - mean(temp[dim(temp)[1] - 1,"Mi",], na.rm=T)
        )
        
        hysteresis_outputs <- rbind(hysteresis_outputs, temp_frame)
        
      }
      
      output_list[[j]] <- hysteresis_outputs
      
    }
    
    end.time <- Sys.time()
    time.taken <- round(end.time - start.time,2)
    print(paste("This rep took ", time.taken, sep = ""))
    
    outer_output_list[[paste("a", disp_a, "_b", disp_b, "_rc", coral_r, "_rm", macro_r, "_rep", rep, sep = "")]] <- output_list
    
  }
}

#Converts our larger list into a dataframe starting with the first entry
hysteresis_outputs <- rbindlist(outer_output_list[[1]])

#...then looping through the rest
for(i in 2:length(outer_output_list)){
  
  hysteresis_outputs <- rbind(hysteresis_outputs, rbindlist(outer_output_list[[i]]))
  
}

#Writes sensitivity outputs
#write.csv(hysteresis_outputs, file = "outputs\\sensitivity_analysis_coral.csv", row.names = F)

###################
###################
#MACROALGAL GROWTH
###################
###################
#Creates range of values to explore sensitivity
#The default of each spatial parameter up to +- 20%
rm_list <- seq(0.013*.8, 0.013*1.20, by = 0.013*.05) #default: 0.013 +-20%

#Sets other parameters to their default
disp_a <- 0.75
disp_b <- 3
coral_r <- 0.0165

#Loops through each rm, analyzing the model
for(rm in rm_list){
  
  macro_r <- rm
  print(paste("Progress ---> macro_r = ", macro_r, " of ", max(rm_list), sep = ""))
  
  #Loops through the different landscapes
  for(rep in unique(sensitivity_reps$unique_id_rep)){
    
    #Sets the generation values to that of this landscape
    temp_landscape <- sensitivity_reps %>%
      filter(unique_id_rep == rep)
    
    set.seed(temp_landscape$rand.seed)
    
    if(temp_landscape$prop_cluster == 0){
      
      input_landscape <- generate_inputs_from_random_landscape(nrow = sqrt(temp_landscape$size), 
                                                               ncol = sqrt(temp_landscape$size), 
                                                               resolution = temp_landscape$res,
                                                               prop_cluster = temp_landscape$prop_cluster, 
                                                               prop_structure = temp_landscape$prop_structure)
      
    } else{
      
      input_landscape <- generate_inputs_from_random_landscape_no_cluster(nrow = sqrt(temp_landscape$size), 
                                                                          ncol = sqrt(temp_landscape$size), 
                                                                          resolution = temp_landscape$res,
                                                                          prop_structure = temp_landscape$prop_structure)
      
    }
    
    #Sets our herbivory values
    Dv_list <- c(0.05, 0.1, 0.15, 0.2, seq(0.25, 8, by = 0.25))
    julia_assign("Dv_list", Dv_list)
    #Sets our function in Julia with all our parameters
    prob_func <- julia_eval("function prob_func(prob, i, repeat)
      remake(prob, p = VectorOfArray([Rc, Rm, pm, Gtc, Gtv, Gti, y, Dc, Dv_list[i], Di, w, neighbor_mat, dispersal_mat, C_growth, Mv_growth, Mi_growth, Mi_disp, turf]))
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
    
    ###############
    #Loads in dispersal matrices for this landscape and saves to Julia
    ###############
    neighbor_mat_m <- initial_grid_high_C[[1]] #Neighbor dispersal for macroalgae
    neighbor_mat_c <- initial_grid_high_C[[2]] #Neighbor dispersal for coral
    dispersal_mat <- initial_grid_high_C[[3]] #Dispersal matrix for Mi -> Mv
    
    julia_assign("neighbor_mat_m", neighbor_mat_m)
    julia_assign("neighbor_mat_c", neighbor_mat_c)
    julia_assign("dispersal_mat", dispersal_mat)
    
    #####################
    #RUNS THE MODEL ON VERY HIGH AND LOW HERBIVORY TO LEARN DIFFERENT EQUILIBRIA
    #####################
    start.time <- Sys.time()
    
    #Sets dV high
    julia_assign("Dv", 12)
    
    #Initializes with high coral
    julia_assign("init_C", model_initializations[[1]][[4]][1,])
    julia_assign("init_Mv", model_initializations[[1]][[4]][2,])
    julia_assign("init_Mi", model_initializations[[1]][[4]][3,])
    
    #Runs the model
    julia_source("scripts\\spatial_briggs_for_R.jl")
    prob <- julia_eval("ODEProblem(f, u0, tspan, p)")
    temp <- solve_julia_model(prob)
    dimnames(temp)[[3]] <- which(input_grids[[4]] == 1)
    #Saves output
    high_c_init_high_c <- mean(temp[nrow(temp),'C',])
    low_m_init_high_c <- mean(temp[nrow(temp),'Mi',] + temp[nrow(temp),'Mv',])
    
    #Initializes with low coral
    julia_assign("init_C", model_initializations[[1]][[4]][1,])
    julia_assign("init_Mv", model_initializations[[1]][[4]][2,])
    julia_assign("init_Mi", model_initializations[[1]][[4]][3,])
    
    #Runs the model
    julia_source("scripts\\spatial_briggs_for_R.jl")
    prob <- julia_eval("ODEProblem(f, u0, tspan, p)")
    temp <- solve_julia_model(prob)
    dimnames(temp)[[3]] <- which(input_grids[[4]] == 1)
    #Saves output
    high_c_init_low_c <- mean(temp[nrow(temp),'C',])
    low_m_init_low_c <- mean(temp[nrow(temp),'Mi',] + temp[nrow(temp),'Mv',])
    
    #Runs the model on this landscape at very low herbivory under both starting conditions
    julia_assign("Dv", 0.01)
    #Initializes with high coral
    julia_assign("init_C", model_initializations[[1]][[4]][1,])
    julia_assign("init_Mv", model_initializations[[1]][[4]][2,])
    julia_assign("init_Mi", model_initializations[[1]][[4]][3,])
    
    #Runs the model
    julia_source("scripts\\spatial_briggs_for_R.jl")
    prob <- julia_eval("ODEProblem(f, u0, tspan, p)")
    temp <- solve_julia_model(prob)
    dimnames(temp)[[3]] <- which(input_grids[[4]] == 1)
    #Saves output
    low_c_init_high_c <- mean(temp[nrow(temp),'C',])
    high_m_init_high_c <- mean(temp[nrow(temp),'Mi',] + temp[nrow(temp),'Mv',])
    
    #Initializes with low coral
    julia_assign("init_C", model_initializations[[1]][[4]][1,])
    julia_assign("init_Mv", model_initializations[[1]][[4]][2,])
    julia_assign("init_Mi", model_initializations[[1]][[4]][3,])
    
    #Runs the model
    julia_source("scripts\\spatial_briggs_for_R.jl")
    prob <- julia_eval("ODEProblem(f, u0, tspan, p)")
    temp <- solve_julia_model(prob)
    dimnames(temp)[[3]] <- which(input_grids[[4]] == 1)
    #Saves output
    low_c_init_low_c <- mean(temp[nrow(temp),'C',])
    high_m_init_low_c <- mean(temp[nrow(temp),'Mi',] + temp[nrow(temp),'Mv',])
    
    output_list <- list()
    
    #Loops through the two initial conditions
    for(j in 1:2){
      hysteresis_outputs <- data.frame("unique_id_rep" = character(),
                                       "alpha" = numeric(),
                                       "beta" = numeric(),
                                       "coral_r" = numeric(),
                                       "macro_r" = numeric(),
                                       "high_initial_var" = character(),
                                       "dV" = numeric(),
                                       "years_ran" = numeric(),
                                       "end_C" = numeric(),
                                       "sd_C" = numeric(),
                                       "high_C_init_high_C" = numeric(),
                                       "low_C_init_high_C" = numeric(),
                                       "high_C_init_low_C" = numeric(),
                                       "low_C_init_low_C" = numeric(),
                                       "prop_C" = numeric(),
                                       "end_Mi" = numeric(),
                                       "sd_Mi" = numeric(),
                                       "end_Mv" = numeric(),
                                       "sd_Mv" = numeric(),
                                       "high_M_init_high_C" = numeric(),
                                       "low_M_init_high_C" = numeric(),
                                       "high_M_init_low_C" = numeric(),
                                       "low_M_init_low_C" = numeric(),
                                       "prop_M" = numeric(),
                                       "prop_ratio" = numeric(),
                                       "delta_C" = numeric(),
                                       "delta_Mi" = numeric())
      
      #j goes through the high and low initial conditions and runs the model
      julia_assign("init_C", model_initializations[[j]][[4]][1,])
      julia_assign("init_Mv", model_initializations[[j]][[4]][2,])
      julia_assign("init_Mi", model_initializations[[j]][[4]][3,])
      
      #Runs the model
      julia_source("scripts\\spatial_briggs_for_R.jl")
      prob <- julia_eval("ODEProblem(f, u0, tspan, p)")
      ensembleprob <- de$EnsembleProblem(prob, prob_func = prob_func)
      
      #Solves in Julia and brings outputs back into R
      sol <- de$solve(ensembleprob,de$Tsit5(),de$EnsembleSerial(),trajectories=length(Dv_list))
      solves <- sapply(sol$u,identity)
      
      #Goes through each solve (herbivory value) and extracts data
      for(i in 1:length(Dv_list)){
        
        mat <- solves[[i]]
        mat <- sapply(mat$u,identity)
        
        model_outputs <- array(dim = c(ncol(mat), 4, dim(mat)[1]/3))
        colnames(model_outputs) <- c("C", "Mv", "Mi", "t")
        model_outputs[,"C",] <- t(mat[seq(1, nrow(mat), by = 3),])
        model_outputs[,"Mv",] <- t(mat[seq(2, nrow(mat), by = 3),])
        model_outputs[,"Mi",] <- t(mat[seq(3, nrow(mat), by = 3),])
        model_outputs[,"t",] <- (1 - (model_outputs[,'C',] + model_outputs[,'Mv',] + model_outputs[,'Mi',]))
        
        temp <- model_outputs
        dimnames(temp)[[3]] <- which(input_grids[[4]] == 1)
        
        temp_frame <- data.frame("unique_id_rep" = temp_landscape$unique_id_rep,
                                 "alpha" = disp_a,
                                 "beta" = disp_b,
                                 "coral_r" = coral_r,
                                 "macro_r" = macro_r,
                                 "high_initial_var" = ifelse(j == 1, "C", "Mi"),
                                 "dV" = Dv_list[i],
                                 "years_ran" = years_to_run,
                                 "end_C" = ifelse(dim(temp)[3] == 1, temp[nrow(temp), "C", 1], rowMeans(temp[,"C",], na.rm=T)[length(rowMeans(temp[,"C",]))]),
                                 "sd_C" = ifelse(dim(temp)[3] == 1, 0, apply(temp[,"C",], 1, FUN = sd, na.rm=T)[length(apply(temp[,"C",], 1, FUN = sd))]),
                                 "high_C_init_high_C" = high_c_init_high_c,
                                 "low_C_init_high_C" = low_c_init_high_c,
                                 "high_C_init_low_C" = high_c_init_low_c,
                                 "low_C_init_low_C" = low_c_init_low_c,
                                 "prop_C" = ifelse(j == 1, sum(temp[nrow(temp),'C',] > (0.85 * high_c_init_high_c)) / dim(temp)[3],
                                                   sum(temp[nrow(temp),'C',] > (0.85 * high_c_init_low_c)) / dim(temp)[3]),
                                 "end_Mi" = ifelse(dim(temp)[3] == 1, temp[nrow(temp), "Mi", 1], rowMeans(temp[,"Mi",], na.rm=T)[length(rowMeans(temp[,"Mi",]))]),
                                 "sd_Mi" = ifelse(dim(temp)[3] == 1, 0, apply(temp[,"Mi",], 1, FUN = sd, na.rm=T)[length(apply(temp[,"Mi",], 1, FUN = sd))]),
                                 "end_Mv" = ifelse(dim(temp)[3] == 1, temp[nrow(temp), "Mv", 1], rowMeans(temp[,"Mv",], na.rm=T)[length(rowMeans(temp[,"Mv",]))]),
                                 "sd_Mv" = ifelse(dim(temp)[3] == 1, 0, apply(temp[,"Mv",], 1, FUN = sd, na.rm=T)[length(apply(temp[,"Mv",], 1, FUN = sd))]),
                                 "high_M_init_high_C" = high_m_init_high_c,
                                 "low_M_init_high_C" = low_m_init_high_c,
                                 "high_M_init_low_C" = high_m_init_low_c,
                                 "low_M_init_low_C" = low_m_init_low_c,
                                 "prop_M" = ifelse(j == 1, sum((temp[nrow(temp),'Mi',] + temp[nrow(temp), 'Mv',]) > (0.85 * high_m_init_high_c)) / dim(temp)[3],
                                                   sum((temp[nrow(temp),'Mi',] + temp[nrow(temp), 'Mv',]) > (0.85 * high_m_init_low_c)) / dim(temp)[3]),
                                 "prop_ratio_C" =  sum(temp[nrow(temp),'C',] > (temp[nrow(temp),'Mi',] + temp[nrow(temp), 'Mv',])) / dim(temp)[3],
                                 "prop_ratio_M" =  sum(temp[nrow(temp),'C',] < (temp[nrow(temp),'Mi',] + temp[nrow(temp), 'Mv',])) / dim(temp)[3],
                                 "delta_C" = mean(temp[dim(temp)[1],"C",], na.rm=T) - mean(temp[dim(temp)[1] - 1,"C",], na.rm=T),
                                 "delta_Mi" = mean(temp[dim(temp)[1],"Mi",], na.rm=T) - mean(temp[dim(temp)[1] - 1,"Mi",], na.rm=T)
        )
        
        hysteresis_outputs <- rbind(hysteresis_outputs, temp_frame)
        
      }
      
      output_list[[j]] <- hysteresis_outputs
      
    }
    
    end.time <- Sys.time()
    time.taken <- round(end.time - start.time,2)
    print(paste("This rep took ", time.taken, sep = ""))
    
    outer_output_list[[paste("a", disp_a, "_b", disp_b, "_rc", coral_r, "_rm", macro_r, "_rep", rep, sep = "")]] <- output_list
    
  }
}

#Converts our larger list into a dataframe starting with the first entry
hysteresis_outputs <- rbindlist(outer_output_list[[1]])

#...then looping through the rest
for(i in 2:length(outer_output_list)){
  
  hysteresis_outputs <- rbind(hysteresis_outputs, rbindlist(outer_output_list[[i]]))
  
}

#Writes sensitivity outputs
#write.csv(hysteresis_outputs, file = "outputs\\sensitivity_analysis_macroalgae.csv", row.names = F)

###################
###################
#ALPHA
###################
###################
#Creates range of values to explore sensitivity
#The default of each spatial parameter up to +- 20%
alpha_list <- seq(0.75*.8, 0.75*1.20, by = 0.75*.05) #default: 0.75 +-20%

#Sets other parameters to their default
disp_b <- 3
macro_r <- 0.013
coral_r <- 0.0165

#Loops through each alpha, analyzing the model
for(alpha_val in alpha_list){
  
  disp_a <- alpha_val
  print(paste("Progress ---> alpha = ", disp_a, " of ", max(alpha_list), sep = ""))
  
  #Loops through the different landscapes
  for(rep in unique(sensitivity_reps$unique_id_rep)){
    
    #Sets the generation values to that of this landscape
    temp_landscape <- sensitivity_reps %>%
      filter(unique_id_rep == rep)
    
    set.seed(temp_landscape$rand.seed)
    
    if(temp_landscape$prop_cluster == 0){
      
      input_landscape <- generate_inputs_from_random_landscape(nrow = sqrt(temp_landscape$size), 
                                                               ncol = sqrt(temp_landscape$size), 
                                                               resolution = temp_landscape$res,
                                                               prop_cluster = temp_landscape$prop_cluster, 
                                                               prop_structure = temp_landscape$prop_structure)
      
    } else{
      
      input_landscape <- generate_inputs_from_random_landscape_no_cluster(nrow = sqrt(temp_landscape$size), 
                                                                          ncol = sqrt(temp_landscape$size), 
                                                                          resolution = temp_landscape$res,
                                                                          prop_structure = temp_landscape$prop_structure)
      
    }
    
    #Sets our herbivory values
    Dv_list <- c(0.05, 0.1, 0.15, 0.2, seq(0.25, 8, by = 0.25))
    julia_assign("Dv_list", Dv_list)
    #Sets our function in Julia with all our parameters
    prob_func <- julia_eval("function prob_func(prob, i, repeat)
      remake(prob, p = VectorOfArray([Rc, Rm, pm, Gtc, Gtv, Gti, y, Dc, Dv_list[i], Di, w, neighbor_mat, dispersal_mat, C_growth, Mv_growth, Mi_growth, Mi_disp, turf]))
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
    
    ###############
    #Loads in dispersal matrices for this landscape and saves to Julia
    ###############
    neighbor_mat_m <- initial_grid_high_C[[1]] #Neighbor dispersal for macroalgae
    neighbor_mat_c <- initial_grid_high_C[[2]] #Neighbor dispersal for coral
    dispersal_mat <- initial_grid_high_C[[3]] #Dispersal matrix for Mi -> Mv
    
    julia_assign("neighbor_mat_m", neighbor_mat_m)
    julia_assign("neighbor_mat_c", neighbor_mat_c)
    julia_assign("dispersal_mat", dispersal_mat)
    
    #####################
    #RUNS THE MODEL ON VERY HIGH AND LOW HERBIVORY TO LEARN DIFFERENT EQUILIBRIA
    #####################
    start.time <- Sys.time()
    
    #Sets dV high
    julia_assign("Dv", 12)
    
    #Initializes with high coral
    julia_assign("init_C", model_initializations[[1]][[4]][1,])
    julia_assign("init_Mv", model_initializations[[1]][[4]][2,])
    julia_assign("init_Mi", model_initializations[[1]][[4]][3,])
    
    #Runs the model
    julia_source("scripts\\spatial_briggs_for_R.jl")
    prob <- julia_eval("ODEProblem(f, u0, tspan, p)")
    temp <- solve_julia_model(prob)
    dimnames(temp)[[3]] <- which(input_grids[[4]] == 1)
    #Saves output
    high_c_init_high_c <- mean(temp[nrow(temp),'C',])
    low_m_init_high_c <- mean(temp[nrow(temp),'Mi',] + temp[nrow(temp),'Mv',])
    
    #Initializes with low coral
    julia_assign("init_C", model_initializations[[1]][[4]][1,])
    julia_assign("init_Mv", model_initializations[[1]][[4]][2,])
    julia_assign("init_Mi", model_initializations[[1]][[4]][3,])
    
    #Runs the model
    julia_source("scripts\\spatial_briggs_for_R.jl")
    prob <- julia_eval("ODEProblem(f, u0, tspan, p)")
    temp <- solve_julia_model(prob)
    dimnames(temp)[[3]] <- which(input_grids[[4]] == 1)
    #Saves output
    high_c_init_low_c <- mean(temp[nrow(temp),'C',])
    low_m_init_low_c <- mean(temp[nrow(temp),'Mi',] + temp[nrow(temp),'Mv',])
    
    #Runs the model on this landscape at very low herbivory under both starting conditions
    julia_assign("Dv", 0.01)
    #Initializes with high coral
    julia_assign("init_C", model_initializations[[1]][[4]][1,])
    julia_assign("init_Mv", model_initializations[[1]][[4]][2,])
    julia_assign("init_Mi", model_initializations[[1]][[4]][3,])
    
    #Runs the model
    julia_source("scripts\\spatial_briggs_for_R.jl")
    prob <- julia_eval("ODEProblem(f, u0, tspan, p)")
    temp <- solve_julia_model(prob)
    dimnames(temp)[[3]] <- which(input_grids[[4]] == 1)
    #Saves output
    low_c_init_high_c <- mean(temp[nrow(temp),'C',])
    high_m_init_high_c <- mean(temp[nrow(temp),'Mi',] + temp[nrow(temp),'Mv',])
    
    #Initializes with low coral
    julia_assign("init_C", model_initializations[[1]][[4]][1,])
    julia_assign("init_Mv", model_initializations[[1]][[4]][2,])
    julia_assign("init_Mi", model_initializations[[1]][[4]][3,])
    
    #Runs the model
    julia_source("scripts\\spatial_briggs_for_R.jl")
    prob <- julia_eval("ODEProblem(f, u0, tspan, p)")
    temp <- solve_julia_model(prob)
    dimnames(temp)[[3]] <- which(input_grids[[4]] == 1)
    #Saves output
    low_c_init_low_c <- mean(temp[nrow(temp),'C',])
    high_m_init_low_c <- mean(temp[nrow(temp),'Mi',] + temp[nrow(temp),'Mv',])
    
    output_list <- list()
    
    #Loops through the two initial conditions
    for(j in 1:2){
      hysteresis_outputs <- data.frame("unique_id_rep" = character(),
                                       "alpha" = numeric(),
                                       "beta" = numeric(),
                                       "coral_r" = numeric(),
                                       "macro_r" = numeric(),
                                       "high_initial_var" = character(),
                                       "dV" = numeric(),
                                       "years_ran" = numeric(),
                                       "end_C" = numeric(),
                                       "sd_C" = numeric(),
                                       "high_C_init_high_C" = numeric(),
                                       "low_C_init_high_C" = numeric(),
                                       "high_C_init_low_C" = numeric(),
                                       "low_C_init_low_C" = numeric(),
                                       "prop_C" = numeric(),
                                       "end_Mi" = numeric(),
                                       "sd_Mi" = numeric(),
                                       "end_Mv" = numeric(),
                                       "sd_Mv" = numeric(),
                                       "high_M_init_high_C" = numeric(),
                                       "low_M_init_high_C" = numeric(),
                                       "high_M_init_low_C" = numeric(),
                                       "low_M_init_low_C" = numeric(),
                                       "prop_M" = numeric(),
                                       "prop_ratio" = numeric(),
                                       "delta_C" = numeric(),
                                       "delta_Mi" = numeric())
      
      #j goes through the high and low initial conditions and runs the model
      julia_assign("init_C", model_initializations[[j]][[4]][1,])
      julia_assign("init_Mv", model_initializations[[j]][[4]][2,])
      julia_assign("init_Mi", model_initializations[[j]][[4]][3,])
      
      #Runs the model
      julia_source("scripts\\spatial_briggs_for_R.jl")
      prob <- julia_eval("ODEProblem(f, u0, tspan, p)")
      ensembleprob <- de$EnsembleProblem(prob, prob_func = prob_func)
      
      #Solves in Julia and brings outputs back into R
      sol <- de$solve(ensembleprob,de$Tsit5(),de$EnsembleSerial(),trajectories=length(Dv_list))
      solves <- sapply(sol$u,identity)
      
      #Goes through each solve (herbivory value) and extracts data
      for(i in 1:length(Dv_list)){
        
        mat <- solves[[i]]
        mat <- sapply(mat$u,identity)
        
        model_outputs <- array(dim = c(ncol(mat), 4, dim(mat)[1]/3))
        colnames(model_outputs) <- c("C", "Mv", "Mi", "t")
        model_outputs[,"C",] <- t(mat[seq(1, nrow(mat), by = 3),])
        model_outputs[,"Mv",] <- t(mat[seq(2, nrow(mat), by = 3),])
        model_outputs[,"Mi",] <- t(mat[seq(3, nrow(mat), by = 3),])
        model_outputs[,"t",] <- (1 - (model_outputs[,'C',] + model_outputs[,'Mv',] + model_outputs[,'Mi',]))
        
        temp <- model_outputs
        dimnames(temp)[[3]] <- which(input_grids[[4]] == 1)
        
        temp_frame <- data.frame("unique_id_rep" = temp_landscape$unique_id_rep,
                                 "alpha" = disp_a,
                                 "beta" = disp_b,
                                 "coral_r" = coral_r,
                                 "macro_r" = macro_r,
                                 "high_initial_var" = ifelse(j == 1, "C", "Mi"),
                                 "dV" = Dv_list[i],
                                 "years_ran" = years_to_run,
                                 "end_C" = ifelse(dim(temp)[3] == 1, temp[nrow(temp), "C", 1], rowMeans(temp[,"C",], na.rm=T)[length(rowMeans(temp[,"C",]))]),
                                 "sd_C" = ifelse(dim(temp)[3] == 1, 0, apply(temp[,"C",], 1, FUN = sd, na.rm=T)[length(apply(temp[,"C",], 1, FUN = sd))]),
                                 "high_C_init_high_C" = high_c_init_high_c,
                                 "low_C_init_high_C" = low_c_init_high_c,
                                 "high_C_init_low_C" = high_c_init_low_c,
                                 "low_C_init_low_C" = low_c_init_low_c,
                                 "prop_C" = ifelse(j == 1, sum(temp[nrow(temp),'C',] > (0.85 * high_c_init_high_c)) / dim(temp)[3],
                                                   sum(temp[nrow(temp),'C',] > (0.85 * high_c_init_low_c)) / dim(temp)[3]),
                                 "end_Mi" = ifelse(dim(temp)[3] == 1, temp[nrow(temp), "Mi", 1], rowMeans(temp[,"Mi",], na.rm=T)[length(rowMeans(temp[,"Mi",]))]),
                                 "sd_Mi" = ifelse(dim(temp)[3] == 1, 0, apply(temp[,"Mi",], 1, FUN = sd, na.rm=T)[length(apply(temp[,"Mi",], 1, FUN = sd))]),
                                 "end_Mv" = ifelse(dim(temp)[3] == 1, temp[nrow(temp), "Mv", 1], rowMeans(temp[,"Mv",], na.rm=T)[length(rowMeans(temp[,"Mv",]))]),
                                 "sd_Mv" = ifelse(dim(temp)[3] == 1, 0, apply(temp[,"Mv",], 1, FUN = sd, na.rm=T)[length(apply(temp[,"Mv",], 1, FUN = sd))]),
                                 "high_M_init_high_C" = high_m_init_high_c,
                                 "low_M_init_high_C" = low_m_init_high_c,
                                 "high_M_init_low_C" = high_m_init_low_c,
                                 "low_M_init_low_C" = low_m_init_low_c,
                                 "prop_M" = ifelse(j == 1, sum((temp[nrow(temp),'Mi',] + temp[nrow(temp), 'Mv',]) > (0.85 * high_m_init_high_c)) / dim(temp)[3],
                                                   sum((temp[nrow(temp),'Mi',] + temp[nrow(temp), 'Mv',]) > (0.85 * high_m_init_low_c)) / dim(temp)[3]),
                                 "prop_ratio_C" =  sum(temp[nrow(temp),'C',] > (temp[nrow(temp),'Mi',] + temp[nrow(temp), 'Mv',])) / dim(temp)[3],
                                 "prop_ratio_M" =  sum(temp[nrow(temp),'C',] < (temp[nrow(temp),'Mi',] + temp[nrow(temp), 'Mv',])) / dim(temp)[3],
                                 "delta_C" = mean(temp[dim(temp)[1],"C",], na.rm=T) - mean(temp[dim(temp)[1] - 1,"C",], na.rm=T),
                                 "delta_Mi" = mean(temp[dim(temp)[1],"Mi",], na.rm=T) - mean(temp[dim(temp)[1] - 1,"Mi",], na.rm=T)
        )
        
        hysteresis_outputs <- rbind(hysteresis_outputs, temp_frame)
        
      }
      
      output_list[[j]] <- hysteresis_outputs
      
    }
    
    end.time <- Sys.time()
    time.taken <- round(end.time - start.time,2)
    print(paste("This rep took ", time.taken, sep = ""))
    
    outer_output_list[[paste("a", disp_a, "_b", disp_b, "_rc", coral_r, "_rm", macro_r, "_rep", rep, sep = "")]] <- output_list
    
  }
}

#Converts our larger list into a dataframe starting with the first entry
hysteresis_outputs <- rbindlist(outer_output_list[[1]])

#...then looping through the rest
for(i in 2:length(outer_output_list)){
  
  hysteresis_outputs <- rbind(hysteresis_outputs, rbindlist(outer_output_list[[i]]))
  
}

#Writes sensitivity outputs
#write.csv(hysteresis_outputs, file = "outputs\\sensitivity_analysis_alpha.csv", row.names = F)

###################
###################
#BETA
###################
###################
#Creates range of values to explore sensitivity
#The default of each spatial parameter up to +- 20%
beta_list <- seq(3*.8, 3*1.20, by = 3*.05) #default: 3 +-20%

#Sets other parameters to their default
disp_a <- 0.75
macro_r <- 0.013
coral_r <- 0.0165

#Loops through each beta, analyzing the model
for(beta_val in beta_list){
  
  disp_b <- beta_val
  print(paste("Progress ---> beta = ", disp_b, " of ", max(beta_list), sep = ""))
  
  #Loops through the different landscapes
  for(rep in unique(sensitivity_reps$unique_id_rep)){
    
    #Sets the generation values to that of this landscape
    temp_landscape <- sensitivity_reps %>%
      filter(unique_id_rep == rep)
    
    set.seed(temp_landscape$rand.seed)
    
    if(temp_landscape$prop_cluster == 0){
      
      input_landscape <- generate_inputs_from_random_landscape(nrow = sqrt(temp_landscape$size), 
                                                               ncol = sqrt(temp_landscape$size), 
                                                               resolution = temp_landscape$res,
                                                               prop_cluster = temp_landscape$prop_cluster, 
                                                               prop_structure = temp_landscape$prop_structure)
      
    } else{
      
      input_landscape <- generate_inputs_from_random_landscape_no_cluster(nrow = sqrt(temp_landscape$size), 
                                                                          ncol = sqrt(temp_landscape$size), 
                                                                          resolution = temp_landscape$res,
                                                                          prop_structure = temp_landscape$prop_structure)
      
    }
    
    #Sets our herbivory values
    Dv_list <- c(0.05, 0.1, 0.15, 0.2, seq(0.25, 8, by = 0.25))
    julia_assign("Dv_list", Dv_list)
    #Sets our function in Julia with all our parameters
    prob_func <- julia_eval("function prob_func(prob, i, repeat)
      remake(prob, p = VectorOfArray([Rc, Rm, pm, Gtc, Gtv, Gti, y, Dc, Dv_list[i], Di, w, neighbor_mat, dispersal_mat, C_growth, Mv_growth, Mi_growth, Mi_disp, turf]))
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
    
    ###############
    #Loads in dispersal matrices for this landscape and saves to Julia
    ###############
    neighbor_mat_m <- initial_grid_high_C[[1]] #Neighbor dispersal for macroalgae
    neighbor_mat_c <- initial_grid_high_C[[2]] #Neighbor dispersal for coral
    dispersal_mat <- initial_grid_high_C[[3]] #Dispersal matrix for Mi -> Mv
    
    julia_assign("neighbor_mat_m", neighbor_mat_m)
    julia_assign("neighbor_mat_c", neighbor_mat_c)
    julia_assign("dispersal_mat", dispersal_mat)
    
    #####################
    #RUNS THE MODEL ON VERY HIGH AND LOW HERBIVORY TO LEARN DIFFERENT EQUILIBRIA
    #####################
    start.time <- Sys.time()
    
    #Sets dV high
    julia_assign("Dv", 12)
    
    #Initializes with high coral
    julia_assign("init_C", model_initializations[[1]][[4]][1,])
    julia_assign("init_Mv", model_initializations[[1]][[4]][2,])
    julia_assign("init_Mi", model_initializations[[1]][[4]][3,])
    
    #Runs the model
    julia_source("scripts\\spatial_briggs_for_R.jl")
    prob <- julia_eval("ODEProblem(f, u0, tspan, p)")
    temp <- solve_julia_model(prob)
    dimnames(temp)[[3]] <- which(input_grids[[4]] == 1)
    #Saves output
    high_c_init_high_c <- mean(temp[nrow(temp),'C',])
    low_m_init_high_c <- mean(temp[nrow(temp),'Mi',] + temp[nrow(temp),'Mv',])
    
    #Initializes with low coral
    julia_assign("init_C", model_initializations[[1]][[4]][1,])
    julia_assign("init_Mv", model_initializations[[1]][[4]][2,])
    julia_assign("init_Mi", model_initializations[[1]][[4]][3,])
    
    #Runs the model
    julia_source("scripts\\spatial_briggs_for_R.jl")
    prob <- julia_eval("ODEProblem(f, u0, tspan, p)")
    temp <- solve_julia_model(prob)
    dimnames(temp)[[3]] <- which(input_grids[[4]] == 1)
    #Saves output
    high_c_init_low_c <- mean(temp[nrow(temp),'C',])
    low_m_init_low_c <- mean(temp[nrow(temp),'Mi',] + temp[nrow(temp),'Mv',])
    
    #Runs the model on this landscape at very low herbivory under both starting conditions
    julia_assign("Dv", 0.01)
    #Initializes with high coral
    julia_assign("init_C", model_initializations[[1]][[4]][1,])
    julia_assign("init_Mv", model_initializations[[1]][[4]][2,])
    julia_assign("init_Mi", model_initializations[[1]][[4]][3,])
    
    #Runs the model
    julia_source("scripts\\spatial_briggs_for_R.jl")
    prob <- julia_eval("ODEProblem(f, u0, tspan, p)")
    temp <- solve_julia_model(prob)
    dimnames(temp)[[3]] <- which(input_grids[[4]] == 1)
    #Saves output
    low_c_init_high_c <- mean(temp[nrow(temp),'C',])
    high_m_init_high_c <- mean(temp[nrow(temp),'Mi',] + temp[nrow(temp),'Mv',])
    
    #Initializes with low coral
    julia_assign("init_C", model_initializations[[1]][[4]][1,])
    julia_assign("init_Mv", model_initializations[[1]][[4]][2,])
    julia_assign("init_Mi", model_initializations[[1]][[4]][3,])
    
    #Runs the model
    julia_source("scripts\\spatial_briggs_for_R.jl")
    prob <- julia_eval("ODEProblem(f, u0, tspan, p)")
    temp <- solve_julia_model(prob)
    dimnames(temp)[[3]] <- which(input_grids[[4]] == 1)
    #Saves output
    low_c_init_low_c <- mean(temp[nrow(temp),'C',])
    high_m_init_low_c <- mean(temp[nrow(temp),'Mi',] + temp[nrow(temp),'Mv',])
    
    output_list <- list()
    
    #Loops through the two initial conditions
    for(j in 1:2){
      hysteresis_outputs <- data.frame("unique_id_rep" = character(),
                                       "alpha" = numeric(),
                                       "beta" = numeric(),
                                       "coral_r" = numeric(),
                                       "macro_r" = numeric(),
                                       "high_initial_var" = character(),
                                       "dV" = numeric(),
                                       "years_ran" = numeric(),
                                       "end_C" = numeric(),
                                       "sd_C" = numeric(),
                                       "high_C_init_high_C" = numeric(),
                                       "low_C_init_high_C" = numeric(),
                                       "high_C_init_low_C" = numeric(),
                                       "low_C_init_low_C" = numeric(),
                                       "prop_C" = numeric(),
                                       "end_Mi" = numeric(),
                                       "sd_Mi" = numeric(),
                                       "end_Mv" = numeric(),
                                       "sd_Mv" = numeric(),
                                       "high_M_init_high_C" = numeric(),
                                       "low_M_init_high_C" = numeric(),
                                       "high_M_init_low_C" = numeric(),
                                       "low_M_init_low_C" = numeric(),
                                       "prop_M" = numeric(),
                                       "prop_ratio" = numeric(),
                                       "delta_C" = numeric(),
                                       "delta_Mi" = numeric())
      
      #j goes through the high and low initial conditions and runs the model
      julia_assign("init_C", model_initializations[[j]][[4]][1,])
      julia_assign("init_Mv", model_initializations[[j]][[4]][2,])
      julia_assign("init_Mi", model_initializations[[j]][[4]][3,])
      
      #Runs the model
      julia_source("scripts\\spatial_briggs_for_R.jl")
      prob <- julia_eval("ODEProblem(f, u0, tspan, p)")
      ensembleprob <- de$EnsembleProblem(prob, prob_func = prob_func)
      
      #Solves in Julia and brings outputs back into R
      sol <- de$solve(ensembleprob,de$Tsit5(),de$EnsembleSerial(),trajectories=length(Dv_list))
      solves <- sapply(sol$u,identity)
      
      #Goes through each solve (herbivory value) and extracts data
      for(i in 1:length(Dv_list)){
        
        mat <- solves[[i]]
        mat <- sapply(mat$u,identity)
        
        model_outputs <- array(dim = c(ncol(mat), 4, dim(mat)[1]/3))
        colnames(model_outputs) <- c("C", "Mv", "Mi", "t")
        model_outputs[,"C",] <- t(mat[seq(1, nrow(mat), by = 3),])
        model_outputs[,"Mv",] <- t(mat[seq(2, nrow(mat), by = 3),])
        model_outputs[,"Mi",] <- t(mat[seq(3, nrow(mat), by = 3),])
        model_outputs[,"t",] <- (1 - (model_outputs[,'C',] + model_outputs[,'Mv',] + model_outputs[,'Mi',]))
        
        temp <- model_outputs
        dimnames(temp)[[3]] <- which(input_grids[[4]] == 1)
        
        temp_frame <- data.frame("unique_id_rep" = temp_landscape$unique_id_rep,
                                 "alpha" = disp_a,
                                 "beta" = disp_b,
                                 "coral_r" = coral_r,
                                 "macro_r" = macro_r,
                                 "high_initial_var" = ifelse(j == 1, "C", "Mi"),
                                 "dV" = Dv_list[i],
                                 "years_ran" = years_to_run,
                                 "end_C" = ifelse(dim(temp)[3] == 1, temp[nrow(temp), "C", 1], rowMeans(temp[,"C",], na.rm=T)[length(rowMeans(temp[,"C",]))]),
                                 "sd_C" = ifelse(dim(temp)[3] == 1, 0, apply(temp[,"C",], 1, FUN = sd, na.rm=T)[length(apply(temp[,"C",], 1, FUN = sd))]),
                                 "high_C_init_high_C" = high_c_init_high_c,
                                 "low_C_init_high_C" = low_c_init_high_c,
                                 "high_C_init_low_C" = high_c_init_low_c,
                                 "low_C_init_low_C" = low_c_init_low_c,
                                 "prop_C" = ifelse(j == 1, sum(temp[nrow(temp),'C',] > (0.85 * high_c_init_high_c)) / dim(temp)[3],
                                                   sum(temp[nrow(temp),'C',] > (0.85 * high_c_init_low_c)) / dim(temp)[3]),
                                 "end_Mi" = ifelse(dim(temp)[3] == 1, temp[nrow(temp), "Mi", 1], rowMeans(temp[,"Mi",], na.rm=T)[length(rowMeans(temp[,"Mi",]))]),
                                 "sd_Mi" = ifelse(dim(temp)[3] == 1, 0, apply(temp[,"Mi",], 1, FUN = sd, na.rm=T)[length(apply(temp[,"Mi",], 1, FUN = sd))]),
                                 "end_Mv" = ifelse(dim(temp)[3] == 1, temp[nrow(temp), "Mv", 1], rowMeans(temp[,"Mv",], na.rm=T)[length(rowMeans(temp[,"Mv",]))]),
                                 "sd_Mv" = ifelse(dim(temp)[3] == 1, 0, apply(temp[,"Mv",], 1, FUN = sd, na.rm=T)[length(apply(temp[,"Mv",], 1, FUN = sd))]),
                                 "high_M_init_high_C" = high_m_init_high_c,
                                 "low_M_init_high_C" = low_m_init_high_c,
                                 "high_M_init_low_C" = high_m_init_low_c,
                                 "low_M_init_low_C" = low_m_init_low_c,
                                 "prop_M" = ifelse(j == 1, sum((temp[nrow(temp),'Mi',] + temp[nrow(temp), 'Mv',]) > (0.85 * high_m_init_high_c)) / dim(temp)[3],
                                                   sum((temp[nrow(temp),'Mi',] + temp[nrow(temp), 'Mv',]) > (0.85 * high_m_init_low_c)) / dim(temp)[3]),
                                 "prop_ratio_C" =  sum(temp[nrow(temp),'C',] > (temp[nrow(temp),'Mi',] + temp[nrow(temp), 'Mv',])) / dim(temp)[3],
                                 "prop_ratio_M" =  sum(temp[nrow(temp),'C',] < (temp[nrow(temp),'Mi',] + temp[nrow(temp), 'Mv',])) / dim(temp)[3],
                                 "delta_C" = mean(temp[dim(temp)[1],"C",], na.rm=T) - mean(temp[dim(temp)[1] - 1,"C",], na.rm=T),
                                 "delta_Mi" = mean(temp[dim(temp)[1],"Mi",], na.rm=T) - mean(temp[dim(temp)[1] - 1,"Mi",], na.rm=T)
        )
        
        hysteresis_outputs <- rbind(hysteresis_outputs, temp_frame)
        
      }
      
      output_list[[j]] <- hysteresis_outputs
      
    }
    
    end.time <- Sys.time()
    time.taken <- round(end.time - start.time,2)
    print(paste("This rep took ", time.taken, sep = ""))
    
    outer_output_list[[paste("a", disp_a, "_b", disp_b, "_rc", coral_r, "_rm", macro_r, "_rep", rep, sep = "")]] <- output_list
    
  }
}

#Converts our larger list into a dataframe starting with the first entry
hysteresis_outputs <- rbindlist(outer_output_list[[1]])

#...then looping through the rest
for(i in 2:length(outer_output_list)){
  
  hysteresis_outputs <- rbind(hysteresis_outputs, rbindlist(outer_output_list[[i]]))
  
}

#Writes sensitivity outputs
#write.csv(hysteresis_outputs, file = "outputs\\sensitivity_analysis_beta.csv", row.names = F)

