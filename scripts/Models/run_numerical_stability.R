rm(list=ls())

##################
#LOADS IN LIBRARIES
##################
library(landscapemetrics)
library(JuliaCall)
library(data.table)
source("scripts\\model_functions.R")

#################
#SETS UP JULIA
#################
julia_setup()
julia_library("DifferentialEquations, RecursiveArrayTools, LinearAlgebra")
julia_command("using DifferentialEquations, RecursiveArrayTools, LinearAlgebra")
de <- diffeqr::diffeq_setup()

sensitivity_reps <- read.csv("data\\sensitivity_landscape_metadata.csv")

#Defines our system of differential equations in Julia
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
#Nonspatial
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

#Spatial parameters
disp_a <- 0.75
disp_b <- 3
macro_r <- 0.013
coral_r <- 0.0165

#Now assigns parameters to Julia
julia_assign("Rc", Rc)
julia_assign("Rm", Rm)
julia_assign("pm", pm)
julia_assign("Gtc", Gtc)
julia_assign("Gtv", Gtv)
julia_assign("Gti", Gti)
julia_assign("y", y)
julia_assign("Dc", Dc)
julia_assign("Di", Di)
julia_assign("Dv", 12)
julia_assign("w", w)

years_to_run <- 10000.0 #Maximum time steps to run
cell_n <- 15 #Number of cells in each dimension
cell_res <- 0.5 #Resolution of cells (in meters)

tspan <- c(0.0,years_to_run) #Range of potential t
save_int <- 10 #How frequently to save outputs
julia_assign("tspan", tspan)

init_frame <- expand_grid(c_init = seq(0.1, 0.6, by = 0.05),
                          m_init = seq(0.1, 0.6, by = 0.05)) %>%
  filter(c_init + m_init <= 0.9)

outer_output_list <- list() #saves landscape results
outer_time_series <- list()

for(rep in unique(sensitivity_reps$unique_id_rep)){
  print(paste("Progress ---> landscape_ID = ", rep, sep = ""))
  outer_output_list <- list() #saves landscape results
  outer_time_series <- list()
  
  for(trial in 1:nrow(init_frame)){
    print(paste("Progress ---> trial = ", trial, " of ", nrow(init_frame), sep = ""))
    
    #Loads in the landscape from the 
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
    
    #####################
    #CREATES MODEL INITIALIZATIONS AND DV LIST
    #####################
    Dv_list <- c(0.05, 0.1, 0.15, 0.2, seq(0.25, 8, by = 0.25)) #Creates list of dV values to run over
    #Dv_list <- c(2.25)
    julia_assign("Dv_list", Dv_list)
    #Sets up our ensemble problem for Julia
    prob_func <- julia_eval("function prob_func(prob, i, repeat)
    remake(prob, p = VectorOfArray([Rc, Rm, pm, Gtc, Gtv, Gti, y, Dc, Dv_list[i], Di, w, neighbor_mat_m, neighbor_mat_c, dispersal_mat, C_growth, Mi_growth, Mi_disp, turf]))
end")
    
    
    #Creates the initial grids and matrices to start the model
    input_grids <- generate_grids_from_raster(input_landscape, macro_r = macro_r, coral_r = coral_r,
                                              dispersal_a = disp_a, dispersal_b = disp_b, as.torus = T,
                                              lat_rebound = T, disp_kernel = "exp_power",
                                              normalize_self = 0.5)
    cells_with_structure <- input_grids[[5]]
    
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
    
    output_list <- list() #for landscape metrics
    cell_outputs_list <- list() #for cell metrics
    output_time_series <- list()
    
    #Uses initial grids and initializes them based on high and low coral conditions
    initial_grid <- initialize_julia_inputs(input_grids, C = init_frame$c_init[trial], Mv = 0.1, Mi = init_frame$m_init[trial])
    
    #Saves the high and low coral initial models in a list that gets accessed later
    model_initializations <- list(initial_grid)
    
    for(j in 1:1){
      hysteresis_outputs <- data.frame("unique_id_rep" = character(),
                                       "init_C" = numeric(),
                                       "init_M" = numeric(),
                                       "res" = numeric(),
                                       "size" = numeric(),
                                       "trial" = integer(),
                                       "p_struc" = numeric(),
                                       "clump" = numeric(),
                                       "dV" = numeric(),
                                       "years_ran" = numeric(),
                                       "end_C" = numeric(),
                                       "sd_C" = numeric(),
                                       "prop_C" = numeric(),
                                       "end_Mi" = numeric(),
                                       "sd_Mi" = numeric(),
                                       "end_Mv" = numeric(),
                                       "sd_Mv" = numeric(),
                                       "prop_M" = numeric(),
                                       "prop_ratio" = numeric())
      
      cell_outputs <- data.frame("unique_id_rep" = character(),
                                 "init_C" = numeric(),
                                 "init_M" = numeric(),
                                 "res" = numeric(),
                                 "size" = numeric(),
                                 "trial" = integer(),
                                 "p_struc" = numeric(),
                                 "clump" = numeric(),
                                 "dV" = numeric(),
                                 "cell" = integer(),
                                 "end_C" = numeric(),
                                 "end_Mv" = numeric(),
                                 "end_Mi" = numeric(),
                                 "n_cells" = integer())
      
      #j goes through the high and low initial conditions and runs the model
      julia_assign("init_C", model_initializations[[j]][[4]][1,])
      julia_assign("init_Mv", model_initializations[[j]][[4]][2,])
      julia_assign("init_Mi", model_initializations[[j]][[4]][3,])
      
      #Runs the model
      julia_source("scripts\\spatial_briggs_for_R.jl")
      prob <- julia_eval("ODEProblem(f, u0, tspan, p)")
      ensembleprob <- de$EnsembleProblem(prob, prob_func = prob_func)
      
      #Solves the ensemble problem (all outputs for every value of Dv in Dv_list) 
      sol <- de$solve(ensembleprob,de$Tsit5(),de$EnsembleSerial(),trajectories=length(Dv_list))
      solves <- sapply(sol$u,identity)
      
      #Pulls out information we want from each model run
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
                                 "init_C" = init_frame$c_init[trial],
                                 "init_M" = init_frame$m_init[trial],
                                 "res" = cell_res,
                                 "size" = temp_landscape$size,
                                 "trial" = trial,
                                 "p_struc" = temp_landscape$p_struc,
                                 "clump" = temp_landscape$clump,
                                 "dV" = Dv_list[i],
                                 "years_ran" = years_to_run,
                                 "end_C" = ifelse(dim(temp)[3] == 1, temp[nrow(temp), "C", 1], rowMeans(temp[,"C",], na.rm=T)[length(rowMeans(temp[,"C",]))]),
                                 "sd_C" = ifelse(dim(temp)[3] == 1, 0, apply(temp[,"C",], 1, FUN = sd, na.rm=T)[length(apply(temp[,"C",], 1, FUN = sd))]),
                                 "prop_C" = ifelse(j == 1, sum(temp[nrow(temp),'C',] > (0.85 * high_c_init_high_c)) / dim(temp)[3],
                                                   sum(temp[nrow(temp),'C',] > (0.85 * high_c_init_low_c)) / dim(temp)[3]),
                                 "end_Mi" = ifelse(dim(temp)[3] == 1, temp[nrow(temp), "Mi", 1], rowMeans(temp[,"Mi",], na.rm=T)[length(rowMeans(temp[,"Mi",]))]),
                                 "sd_Mi" = ifelse(dim(temp)[3] == 1, 0, apply(temp[,"Mi",], 1, FUN = sd, na.rm=T)[length(apply(temp[,"Mi",], 1, FUN = sd))]),
                                 "end_Mv" = ifelse(dim(temp)[3] == 1, temp[nrow(temp), "Mv", 1], rowMeans(temp[,"Mv",], na.rm=T)[length(rowMeans(temp[,"Mv",]))]),
                                 "sd_Mv" = ifelse(dim(temp)[3] == 1, 0, apply(temp[,"Mv",], 1, FUN = sd, na.rm=T)[length(apply(temp[,"Mv",], 1, FUN = sd))]),
                                 "prop_M" = ifelse(j == 1, sum((temp[nrow(temp),'Mi',] + temp[nrow(temp), 'Mv',]) > (0.85 * high_m_init_high_c)) / dim(temp)[3],
                                                   sum((temp[nrow(temp),'Mi',] + temp[nrow(temp), 'Mv',]) > (0.85 * high_m_init_low_c)) / dim(temp)[3]),
                                 "prop_ratio_C" =  sum(temp[nrow(temp),'C',] > (temp[nrow(temp),'Mi',] + temp[nrow(temp), 'Mv',])) / dim(temp)[3],
                                 "prop_ratio_M" =  sum(temp[nrow(temp),'C',] < (temp[nrow(temp),'Mi',] + temp[nrow(temp), 'Mv',])) / dim(temp)[3]
        )
        
        n_cells <- dim(temp)[3]
        temp_cell_output <- data.frame("unique_id_rep" = temp_landscape$unique_id_rep,
                                       "init_C" = rep(init_frame$c_init[trial], times = n_cells),
                                       "init_M" = rep(init_frame$m_init[trial], times = n_cells),
                                       "res" = rep(cell_res, times = n_cells),
                                       "size" = rep(temp_landscape$size, times = n_cells),
                                       "trial" = rep(trial, times = n_cells),
                                       "p_struc" = rep(temp_landscape$p_struc, times = n_cells),
                                       "clump" = rep(temp_landscape$clump, times = n_cells),
                                       "dV" = rep(Dv_list[i], times = n_cells),
                                       "cell" = which(cells_with_structure == 1),
                                       "end_C" = temp[dim(temp)[1],"C",],
                                       "end_Mv" = temp[dim(temp)[1],"Mv",],
                                       "end_Mi" = temp[dim(temp)[1],"Mi",],
                                       "end_M" = temp[dim(temp)[1],"Mv",] + temp[dim(temp)[1],"Mi",],
                                       "n_cells" = n_cells)
        
        hysteresis_outputs <- rbind(hysteresis_outputs, temp_frame)
        cell_outputs <- rbind(cell_outputs, temp_cell_output)
        
        #Cuts down the time series to every 100 time steps for size reduction
        temp <- temp[seq(1, nrow(temp), by = 50), ,]
        temp_series <- data.frame("time" = 1:nrow(temp),
                                  "C" = rowMeans(temp[,'C',]),
                                  "Mi" = rowMeans(temp[,'Mi',]),
                                  "Mv" = rowMeans(temp[,'Mv',]),
                                  "")
        output_time_series[[paste("dV:", Dv_list[i], sep = "")]] <- temp
        
      }
      
      output_list[[j]] <- hysteresis_outputs
      cell_outputs_list[[j]] <- cell_outputs
      
      
    }

    end.time <- Sys.time()
    time.taken <- round(end.time - start.time,2)
    print(paste("trial = ", trial, "; landscape_ID = ", rep, " took ", time.taken, sep = ""))
    
    #####################
    #Adds to our larger list
    #####################
    
    outer_output_list[[paste("trial", trial, "; landscape_ID = ", rep, sep = "")]] <- output_list
    outer_time_series[[paste("trial", trial, "; landscape_ID = ", rep, sep = "")]] <- output_time_series
    
  }
  #Save the output for this landscape
  #saveRDS(outer_time_series, file = paste("data\\output_time_series_", rep, ".rds", sep = ""))
  #saveRDS(outer_output_list, file = paste("data\\output_list_", rep, ".rds", sep = ""))
  
}