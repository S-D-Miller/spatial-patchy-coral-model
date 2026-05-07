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
Rc <- 0.001 #External recruitment of coral (phi_C in text)
Rm <- 0.0001 #External recruitment of macroalgae (phi_M in text)
pm <- 0.5 #Local production of vulnerable from invulnerable macroalgae (r_M in text)
Gtc <- 0.1 #Growth of coral over turf
Gtv <- 0.2 #Growth of vulnerable macroalgae over turf
Gti <- 0.4 #Growth of invulnerable macroalgae over turf
y <- 0.4 #gamma...Growth of macroalgae over coral rather than bare space
Dc <- 0.05 #mortality of coral
Di <- 0.4 #mortality of adult macroalgae
w <- 2 #omega...maturation of vulnerable to invulnerable macroalgae

#Spatial parameters
disp_a <- 0.75 #alpha..dispersal kernel scale parameter
disp_b <- 3 #beta...dispersal kernel shape parameter
macro_r <- 0.013 #s_M in text...lateral growth of macroalgae
coral_r <- 0.0165 #s_C in text...lateral growth of coral
rms <- 0.5 #How much local production is kept inside focal cell

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
julia_assign("w", w)

years_to_run <- 10000.0 #Maximum time steps to run
cell_n <- 15 #Number of cells in each dimension
cell_res <- 0.5 #Resolution of cells (in meters)

tspan <- c(0.0,years_to_run) #Range of potential t
save_int <- 10 #How frequently to save outputs
julia_assign("tspan", tspan)

ps_list <- seq(0.1, 0.8, by = 0.1)
n_reps <- 10
outer_output_list <- list() #saves landscape results
outer_output_list_cells <- list() #saves cell results
z <- 0

for(ps in ps_list){
  print(paste("Progress ---> ps = ", ps, " of ", max(ps_list), sep = ""))
  for(k in 1:n_reps){
    print(paste("Progress ---> rep = ", k, " of ", n_reps, sep = ""))
    
    set.seed(42 + z) #Start at 42 for the answer to the Ultimate Question of Life, the Universe, and Everything
    #Generates a random landscape without clusters
    input_landscape <- generate_inputs_from_random_landscape_no_cluster(nrow = cell_n, ncol = cell_n, resolution = cell_res,
                                                                        prop_structure = ps)
    
    #####################
    #Calculates landscape metrics on the random landscape
    #####################
    landscape_dat <- as.data.frame(input_landscape)
    
    #Detect if the random landscape has 0 structure cells
    if(sum(landscape_dat$layer) == 0){
      
      continue_searching <- T
      
    } else{
      
      continue_searching <- F
      
    }
    
    #And if it does have 0 cells, have it continue searching
    while(continue_searching){
      
      z <- z + 1
      set.seed(42 + z)
      input_landscape <- generate_inputs_from_random_landscape_no_cluster(nrow = cell_n, ncol = cell_n, resolution = cell_res,
                                                                          prop_structure = ps)
      
      landscape_dat <- as.data.frame(input_landscape)
      
      if(sum(landscape_dat$layer) > 0){
        
        continue_searching <- F
        
      }
      
    }
    
    #Calculates landscape metrics on the landscape
    temp_class <- calculate_lsm(input_landscape, level = "class") %>%
      dplyr::filter(class == 1)
    
    #####################
    #CREATES MODEL INITIALIZATIONS AND DV LIST
    #####################
    Dv_list <- c(0.05, 0.1, 0.15, 0.2, seq(0.25, 8, by = 0.25)) #Creates list of dV values to run over
    julia_assign("Dv_list", Dv_list)
    #Sets up our ensemble problem for Julia
    prob_func <- julia_eval("function prob_func(prob, i, repeat)
    remake(prob, p = VectorOfArray([Rc, Rm, pm, Gtc, Gtv, Gti, y, Dc, Dv_list[i], Di, w, neighbor_mat_m, neighbor_mat_c, dispersal_mat, C_growth, Mi_growth, Mi_disp, turf]))
end")
    
    
    #Creates the initial grids and matrices to start the model
    input_grids <- generate_grids_from_raster(input_landscape, macro_r = macro_r, coral_r = coral_r,
                                              dispersal_a = disp_a, dispersal_b = disp_b, as.torus = T,
                                              lat_rebound = T, disp_kernel = "exp_power",
                                              normalize_self = rms)
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
    
    for(j in 1:2){
      hysteresis_outputs <- data.frame("rand.seed" = numeric(), 
                                       "prop_structure" = numeric(),
                                       "prop_cluster" = numeric(),
                                       "res" = numeric(),
                                       "size" = numeric(),
                                       "rep" = integer(),
                                       "p_struc" = numeric(),
                                       "clump" = numeric(),
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
                                       "prop_ratio" = numeric())
      
      cell_outputs <- data.frame("rand.seed" = numeric(),
                                 "prop_cluster" = numeric(),
                                 "prop_structure" = numeric(),
                                 "res" = numeric(),
                                 "size" = numeric(),
                                 "rep" = integer(),
                                 "p_struc" = numeric(),
                                 "clump" = numeric(),
                                 "high_initial_var" = character(),
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
        
        temp_frame <- data.frame("rand.seed" = z+42,
                                 "prop_cluster" = 0,
                                 "prop_structure" = ps,
                                 "res" = cell_res,
                                 "size" = length(landscape_dat$layer),
                                 "rep" = k,
                                 "p_struc" = sum(landscape_dat$layer == 1) / length(landscape_dat$layer),
                                 "clump" = temp_class$value[temp_class$metric == "clumpy"],
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
                                 "prop_ratio_M" =  sum(temp[nrow(temp),'C',] < (temp[nrow(temp),'Mi',] + temp[nrow(temp), 'Mv',])) / dim(temp)[3]
        )
        
        n_cells <- dim(temp)[3]
        temp_cell_output <- data.frame("rand.seed" = rep(z+42, times = n_cells),
                                       "prop_cluster" = rep(0, times = n_cells),
                                       "prop_structure" = rep(ps, times = n_cells),
                                       "res" = rep(cell_res, times = n_cells),
                                       "size" = rep(length(landscape_dat[,1]), times = n_cells),
                                       "rep" = rep(k, times = n_cells),
                                       "p_struc" = rep(sum(landscape_dat[,1] == 1) / length(landscape_dat[,1]), times = n_cells),
                                       "clump" = rep(temp_class$value[temp_class$metric == "clumpy"], times = n_cells),
                                       "high_initial_var" = rep(ifelse(j == 1, "C", "Mi"), times = n_cells),
                                       "dV" = rep(Dv_list[i], times = n_cells),
                                       "cell" = which(cells_with_structure == 1),
                                       "end_C" = temp[dim(temp)[1],"C",],
                                       "end_Mv" = temp[dim(temp)[1],"Mv",],
                                       "end_Mi" = temp[dim(temp)[1],"Mi",],
                                       "end_M" = temp[dim(temp)[1],"Mv",] + temp[dim(temp)[1],"Mi",],
                                       "n_cells" = n_cells)
        
        hysteresis_outputs <- rbind(hysteresis_outputs, temp_frame)
        cell_outputs <- rbind(cell_outputs, temp_cell_output)
        
      }
      
      output_list[[j]] <- hysteresis_outputs
      cell_outputs_list[[j]] <- cell_outputs
      
    }
    
    end.time <- Sys.time()
    time.taken <- round(end.time - start.time,2)
    print(paste("ps = ", ps, "; rep = ", k, " took ", time.taken, sep = ""))
    
    #####################
    #Adds to our larger list
    #####################
    
    outer_output_list[[paste("ps", ps, "_rep", k, sep = "")]] <- output_list
    outer_output_list_cells[[paste("p0", "_ps", ps, "_rep", k, sep = "")]] <- cell_outputs_list
    
    #Save the current state of outer_output_list in case something happens
    #saveRDS(outer_output_list_cells, file = "outputs\\temp_cell_outputs.rds")
    
    z <- z + 10 #Increase our seed
    
  }
}

#In case you want to preserve all your data for later use
#saveRDS(outer_output_list_cells, file = "outputs\\cell_outputs_list_nocluster.rds")

#Converts our list of outputs into a dataframe
hysteresis_outputs <- rbindlist(outer_output_list[[1]])

for(i in 2:length(outer_output_list)){
  
  hysteresis_outputs <- rbind(hysteresis_outputs, rbindlist(outer_output_list[[i]]))
  
}

hysteresis_outputs$unique_id_rep <- paste(hysteresis_outputs$prop_structure, 
                                          hysteresis_outputs$rep, sep = "_")

outputs_by_cell <- rbindlist(outer_output_list_cells[[1]])

for(i in 2:length(outer_output_list_cells)){
  
  outputs_by_cell <- rbind(outputs_by_cell, rbindlist(outer_output_list_cells[[i]]))
  
}

outputs_by_cell$unique_id_rep <- paste(outputs_by_cell$prop_structure,
                                       outputs_by_cell$prop_cluster,
                                       outputs_by_cell$rep, sep = "_")

#Saves outputs as csv
#write.csv(hysteresis_outputs, file = "outputs\\noncluster_outputs.csv", row.names = F)