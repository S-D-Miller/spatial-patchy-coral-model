library(tidyverse)
library(terra)
library(viridis)
library(som.nn)
library(NLMR)

calc_torus_dist <- function(coors, grid_res){
  #Calculates distances on a torus
  #Arguments are coordinates in x/y space and the resolution of the grid
  x <- coors[, 1]
  dx <- stats::dist(x, diag = TRUE, upper = TRUE)
  y <- coors[, 2]
  dy <- stats::dist(y, diag = TRUE, upper = TRUE)
  max.x <- max(dx) + grid_res
  max.y <- max(dy) + grid_res
  mdx <- max.x - dx
  mdy <- max.y - dy
  dx <- pmin(dx, mdx)
  dy <- pmin(dy, mdy)
  d <- sqrt(dx^2 + dy^2)
  
  #Returns the distance matrix
  return(d)
}

solve_julia_model <- function(prob, saveat = 10){
  #Runs the model through the julia solver
  #prob is the julia ODEproblem (R output from julia_eval)
  #saveat is the interval at which you want to save outputs for time series
  
  #Solves the model
  sol <- de$solve(prob, de$Tsit5(), saveat = save_int)
  #Converts to a matrix output
  mat <- sapply(sol$u,identity)
  #Saves the time series for the state variables
  model_outputs <- array(dim = c(ncol(mat), 5, dim(mat)[1]/3))
  colnames(model_outputs) <- c("time", "C", "Mv", "Mi", "t")
  model_outputs[,"time",] <- matrix(rep(seq(0,max(tspan), by = save_int), times = dim(model_outputs)[3]), 
                                    nrow = dim(model_outputs)[1], ncol = dim(model_outputs)[3])
  model_outputs[,"C",] <- t(mat[seq(1, nrow(mat), by = 3),])
  model_outputs[,"Mv",] <- t(mat[seq(2, nrow(mat), by = 3),])
  model_outputs[,"Mi",] <- t(mat[seq(3, nrow(mat), by = 3),])
  model_outputs[,"t",] <- (1 - (model_outputs[,'C',] + model_outputs[,'Mv',] + model_outputs[,'Mi',]))
  
  temp <- model_outputs
  
  return(temp)
}

initialize_julia_inputs <- function(input_grids,
                                    C = 0.6,
                                    Mi = 0.1,
                                    Mv = 0.1){
  #Initializes the model for Julia
  #input_grids is the output from generate_grids_from_raster()
  #C, Mi, and Mv are initial values for state variables
  
  #Pulls in the spatial grid to initialize
  spat_grid <- input_grids[[1]]
  neighbor_mat_m <- input_grids[[2]] #Neighbor dispersal for macro
  neighbor_mat_c <- input_grids[[3]] #Neighbor dispersal for coral
  dispersal_mat <- input_grids[[4]] #Dispersal matrix for Mi -> Mv
  cells_with_structure <- input_grids[[5]] #Cells that are used in model
  
  #Creates the array used to do calculations
  model_inputs <- array(dim = c(3, dim(spat_grid)[1]*dim(spat_grid)[2]))
  
  #Removes patches without structure to increase efficiency
  removed_cells <- which(cells_with_structure == 0) #Flags entries to remove (b/c on sand)
  
  #Only attempts to remove cells if there are actually cells to remove
  if(length(removed_cells)>0){
    
    neighbor_mat_m <- neighbor_mat_m[-removed_cells, -removed_cells]
    neighbor_mat_c <- neighbor_mat_c[-removed_cells, -removed_cells]
    dispersal_mat <- dispersal_mat[-removed_cells, -removed_cells]
    model_inputs <- model_inputs[, -removed_cells]
    
  } 
  
  #Initializes all cells the same based on function arguments
  #Can change later to add more interesting initilizations
  model_inputs[1,] <- C
  model_inputs[2,] <- Mv
  model_inputs[3,] <- Mi
  
  return(list(neighbor_mat_m, neighbor_mat_c, dispersal_mat, model_inputs))
  
}

generate_inputs_from_random_landscape <- function(nrow, ncol, resolution,
                                                  prop_cluster, prop_structure){
  #Creates a random landscape with clustering
  #nrow/ncol are number of rows/columns
  #resolution is the cell resolution
  #prop_cluster is proportion of cells chosen to cluster
  #prop_structure is the proportion of cells to be reef habitat
  
  random_cluster <- NLMR::nlm_randomcluster(nrow = nrow,
                                            ncol = ncol,
                                            resolution = resolution,
                                            p    = prop_cluster,
                                            ai   = c(1 - prop_structure, prop_structure),
                                            rescale = FALSE)
  
  output_rast <- rast(random_cluster)
  output_rast <- terra::classify(output_rast, matrix(c(1,0,
                                                       2,1), ncol = 2, byrow = T))
  
  return(output_rast)
  
}

generate_inputs_from_random_landscape_no_cluster <- function(nrow, ncol, resolution,
                                                             prop_structure){
  #Creates a random landscape without clustering
  #nrow/ncol are number of rows/columns
  #resolution is the cell resolution
  #prop_structure is the proportion of cells to be reef habitat
  
  random_cluster <- NLMR::nlm_random(nrow = nrow,
                                     ncol = ncol,
                                     resolution = resolution,
                                     rescale = T)
  
  output_rast <- rast(random_cluster)
  output_rast <- terra::classify(output_rast, matrix(c(0,prop_structure,1,
                                                       prop_structure,1,0), ncol = 3, byrow = T))
  
  return(output_rast)
  
}

generate_grids_from_raster <- function(input_rast, macro_r = 0.05,
                                       coral_r = 0.05,
                                       dispersal_a = 0.75, dispersal_b = 3, 
                                       as.torus = T, lat_rebound = T,
                                       disp_kernel = "exp_power",
                                       normalize_self = NULL){
  #Generates landscape grid and connectivity matrices from an input raster
  #Input raster is output raster from the random landscape generation
  #coral_r and macro_r are the spatial growth terms
  #dispersal_a and dispersal_b are the shape parameters for the dispersal kernel
  #as.torus and lat_rebound turn on torus and lateral grow redirecting back if it would be in a sand cell
  #disp_kernel can be "exp_power" or "gaussian" for the kernel type
  #normalize_self is how much production stays within a focal cell
  
  #Creates a test matrix and calculates the total number of cells
  tot_cells <- dim(input_rast)[1] * dim(input_rast)[2]
  
  #Calculates cell spacing for calculating neighborhood interactions
  cell_space <- terra::res(input_rast)[1]
  
  #Saves vector of structure values for each cell
  cells_with_structure <- as.data.frame(input_rast)[,1]
  
  if(as.torus == T){
    
    #Creates distance matrix from centroid of input assuming a torus
    dist_grid <- as.matrix(zapsmall(calc_torus_dist(xyFromCell(input_rast, 1:tot_cells), cell_space)))
    
  }else{
    
    #Creates distance matrix from centroid of input
    dist_grid <- as.matrix(zapsmall(dist(xyFromCell(input_rast, 1:tot_cells), diag = T, upper = T)))
    
  }
  
  ########
  #Neighborhood interactions
  ########
  
  ########
  #INVULERNABLE MACROALGAE
  ########
  #Finds distance of diagonal cell
  diag_dist <- dist_grid[1,ncol(input_rast) + 2]
  
  #Determines the number of adjacent and diagonal cells for each cell
  dist_with_struc <- dist_grid * cells_with_structure
  adj_cells <- apply(dist_with_struc, 2, function(x) sum(x == cell_space))
  diag_cells <- apply(dist_with_struc, 2, function(x) sum(x == diag_dist))
  
  #Calculates probabilities for diagonal and orthogonal neighborhood interactions
  ortho_prob <- (((4*(macro_r/(cell_space * pi))) - (2*(macro_r**2/(cell_space**2*pi))))) / 4
  diag_prob <- ((macro_r**2)/(cell_space**2*pi)) / 4
  self_prob <- 1 - (diag_prob*4) - (ortho_prob*4)
  
  #Creates the neighborhood matrix using distance grid as starting point
  neighbor_mat_m <- dist_grid
  
  #Flags values from cells that are not in neighborhood
  neighbor_mat_m <- ifelse(neighbor_mat_m > diag_dist, 999, neighbor_mat_m)
  #Calculates probability for orthogonal cells
  neighbor_mat_m <- ifelse(neighbor_mat_m == cell_space, ortho_prob, neighbor_mat_m)
  #...and for diagonal
  neighbor_mat_m <- ifelse(neighbor_mat_m == diag_dist, diag_prob, neighbor_mat_m)
  #...and self
  neighbor_mat_m <- ifelse(neighbor_mat_m == 0, self_prob, neighbor_mat_m)
  #...then resets others to 0 probability
  neighbor_mat_m <- ifelse(neighbor_mat_m == 999, 0, neighbor_mat_m)
  
  #Calculates cell lateral dispersal into sand to rebound back into the focal cell if desired
  if(lat_rebound == T){
    
    #Defines a function to find out the new proportion into the focal cell if rebounded
    rebound_dispersal <- function(x){
      
      #Finds the amount of the total production that goes into empty cells and adds them back to the original self probability
      x[x == self_prob] <- x[x == self_prob] + (1 - sum(x * cells_with_structure))
      
    }
    
    #Applies the function over the columns of the neighbor_mat_m
    new_self_probs <- apply(neighbor_mat_m, 2, FUN = rebound_dispersal)
    
    #Sets the diagonal to the new self probabilities and removes dispersal into cells without structure
    diag(neighbor_mat_m) <- new_self_probs
    neighbor_mat_m <- neighbor_mat_m * cells_with_structure
    
    #Note that this allows cells to have < 1 total proportional dispersal (columns don't all sum to 1)
    #This is OK because the only cells affected are those without structure in them and they are removed later in the model
    
  }
  
  ########
  #CORAL
  ########
  #Finds distance of diagonal cell
  diag_dist <- dist_grid[1,ncol(input_rast) + 2]
  
  #Determines the number of adjacent and diagonal cells for each cell
  dist_with_struc <- dist_grid * cells_with_structure
  adj_cells <- apply(dist_with_struc, 2, function(x) sum(x == cell_space))
  diag_cells <- apply(dist_with_struc, 2, function(x) sum(x == diag_dist))
  
  #Calculates probabilities for diagonal and orthogonal neighborhood interactions
  ortho_prob <- (((4*(coral_r/(cell_space * pi))) - (2*(coral_r**2/(cell_space**2*pi))))) / 4
  diag_prob <- ((coral_r**2)/(cell_space**2*pi)) / 4
  self_prob <- 1 - (diag_prob*4) - (ortho_prob*4)
  
  #Creates the neighborhood matrix using distance grid as starting point
  neighbor_mat_c <- dist_grid
  
  #Flags values from cells that are not in neighborhood
  neighbor_mat_c <- ifelse(neighbor_mat_c > diag_dist, 999, neighbor_mat_c)
  #Calculates probability for orthogonal cells
  neighbor_mat_c <- ifelse(neighbor_mat_c == cell_space, ortho_prob, neighbor_mat_c)
  #...and for diagonal
  neighbor_mat_c <- ifelse(neighbor_mat_c == diag_dist, diag_prob, neighbor_mat_c)
  #...and self
  neighbor_mat_c <- ifelse(neighbor_mat_c == 0, self_prob, neighbor_mat_c)
  #...then resets others to 0 probability
  neighbor_mat_c <- ifelse(neighbor_mat_c == 999, 0, neighbor_mat_c)
  
  #Calculates cell lateral dispersal into sand to rebound back into the focal cell if desired
  if(lat_rebound == T){
    
    #Defines a function to find out the new proportion into the focal cell if rebounded
    rebound_dispersal <- function(x){
      
      #Finds the amount of the total production that goes into empty cells and adds them back to the original self probability
      x[x == self_prob] <- x[x == self_prob] + (1 - sum(x * cells_with_structure))
      
    }
    
    #Applies the function over the columns of the neighbor_mat_c
    new_self_probs <- apply(neighbor_mat_c, 2, FUN = rebound_dispersal)
    
    #Sets the diagonal to the new self probabilities and removes dispersal into cells without structure
    diag(neighbor_mat_c) <- new_self_probs
    neighbor_mat_c <- neighbor_mat_c * cells_with_structure
    
    #Note that this allows cells to have < 1 total proportional dispersal (columns don't all sum to 1)
    #This is OK because the only cells affected are those without structure in them and they are removed later in the model
    
  }
  
  ########
  #Longer-distance dispersal
  ########
  #Calculates the number of cells for each distance from a center cell to scale properly
  cent_vect <- vect(matrix(c(jitter(mean(c(xmin(input_rast), xmax(input_rast))), amount = res(input_rast)[1] * .5), 
                             jitter(mean(c(ymin(input_rast), ymax(input_rast))), amount = res(input_rast)[1] * .5)), 
                           nrow = 1, ncol = 2), type = "points") #Puts point on one of center cells
  cent_int <- terra::extract(input_rast, cent_vect, cells = T) #Extracts cell information for intersected cell
  
  #Calculates probability kernel using deisgnated kernel
  if(disp_kernel == "gaussian"){
    
    dispersal_mat <- (1/(2*pi*dispersal_a**2))*exp(-1*(dist_grid**2/(2*dispersal_a**2)))
    
  }else if(disp_kernel == "exp_power"){
    
    dispersal_mat <- (dispersal_b / (2*pi*gamma(2/dispersal_b)))*exp(-(dist_grid**dispersal_b / dispersal_a**dispersal_b))
    
  }
  
  #Calculates the denominator from a central cell
  standard_mat <- dispersal_mat[,cent_int$cell]
  tot_denom <- sum(standard_mat)
  
  #If there's a value in normalize_self, the code puts that value as the diagonal in dispersal_mat and standardizes off that
  if(is.null(normalize_self)){
    
    #Standardizes the dispersal probabilities based on a central cell denominator
    dispersal_mat <- dispersal_mat / tot_denom
    
  }else{
    
    #Re-normalize the rest of the dispersal to still add up to a total of 1
    dispersal_mat <- (dispersal_mat / tot_denom) * (1 - normalize_self)
    
    #Set diagonal to be normalize_self (each cell has a set probability for itself)
    diag(dispersal_mat) <- normalize_self + diag(dispersal_mat)
    
    
  }
  
  return(list(input_rast, neighbor_mat_m, neighbor_mat_c, dispersal_mat, cells_with_structure))
  
}

calculate_transition_speed_landscape <- function(dat, rep_id, tail_prop = 0.1){
  #Calculates metrics of transition abruptness on landscape-level data
  #dat is output from the model at the landscape level
  #rep_id is the name of the replicate
  #tail_prop is the range you consider to be "transitioning"
  #e.g., tail_prop = 0.1 means values between 10 and 90% of equilibrium are considered "intermediate" and thus part of a gradual transition
  
  #Filters data for ID and starting conditions
  temp <- dat %>%
    filter(unique_id_rep == rep_id)
  temp_C <- dat %>%
    filter(unique_id_rep == rep_id) %>%
    filter(high_initial_var == "C")
  temp_M <- dat %>%
    filter(unique_id_rep == rep_id) %>%
    filter(high_initial_var == "Mi")
  
  #Finds cutoffs for transition values
  high_val <- temp$high_C_init_high_C[1] * (1-tail_prop)
  low_val <- temp$high_C_init_high_C[1] * tail_prop
  
  #Value at which the coral collapses (leaves high coral state) and recovers (leaves low coral state)
  c_collapse <- min(temp_C$dV[temp_C$end_C > high_val])
  c_recover <- max(temp_C$dV[temp_C$end_C < low_val])
  
  #Value at which the macroalgae collapses (leaves high macro state) and recovers (leaves low macro state)
  m_collapse <- min(temp_M$dV[temp_M$end_C > high_val])
  m_recover <- max(temp_M$dV[temp_M$end_C < low_val])
  
  #Saves outputs
  outputs <- suppressWarnings(data.frame("unique_id_rep" = rep_id,
                                         "high_val" = high_val,
                                         "low_val" = low_val,
                                         "C_collapse" = c_collapse,
                                         "C_recover" = c_recover,
                                         "M_collapse" = m_collapse,
                                         "M_recover" = m_recover,
                                         "C_range" = c_collapse - c_recover,
                                         "M_range" = m_collapse - m_recover,
                                         "collapse_diff" = m_collapse - c_collapse,
                                         "recover_diff" = m_recover - c_recover))
  #Puts NA for infinite values
  for(i in 1:ncol(outputs)){
    
    outputs[,i][is.infinite(outputs[,i])] <- NA
    
  }
  
  return(outputs)
  
}

calculate_transition_speed_cell <- function(dat, rep_id, tail_prop = 0.1){
  #Calculates metrics of transition abruptness on cell-level data
  #dat is output from the model at the landscape level
  #rep_id is the name of the replicate
  #tail_prop is the range you consider to be "transitioning"
  #e.g., tail_prop = 0.1 means values between 10 and 90% of equilibrium are considered "intermediate" and thus part of a gradual transition
  
  #Filters data for ID and starting conditions
  temp <- dat %>%
    filter(unique_id_rep == rep_id)
  temp_C <- dat %>%
    filter(unique_id_rep == rep_id) %>%
    filter(high_initial_var == "C")
  temp_M <- dat %>%
    filter(unique_id_rep == rep_id) %>%
    filter(high_initial_var == "Mi")
  
  #Finds cutoffs for transition values
  high_val <- temp$high_C_init_high_C[1] * (1-tail_prop)
  low_val <- temp$high_C_init_high_C[1] * tail_prop
  
  #Gets list of unique cell IDs within the landscape
  cell_list <- unique(temp$cell)
  
  #Sets up our dataframe
  ind_output <- suppressWarnings(data.frame("unique_id_rep" = rep(rep_id, times = length(cell_list)),
                                            "cell_id" = integer(length = length(cell_list)),
                                            "high_val" = rep(high_val, times = length(cell_list)),
                                            "low_val" = rep(low_val, times = length(cell_list)),
                                            "C_collapse" = numeric(length = length(cell_list)),
                                            "C_recover" = numeric(length = length(cell_list)),
                                            "M_collapse" = numeric(length = length(cell_list)),
                                            "M_recover" = numeric(length = length(cell_list)),
                                            "C_range" = numeric(length = length(cell_list)),
                                            "M_range" = numeric(length = length(cell_list)),
                                            "collapse_diff" = numeric(length = length(cell_list)),
                                            "recover_diff" = numeric(length = length(cell_list))))
  
  #Loops through cell_list and makes calculation for each cell, saving in ind_output
  for(i in 1:length(cell_list)){
    
    #Filters for data just for this cell
    temp_cell_C <- temp_C %>%
      filter(cell == cell_list[i])
    temp_cell_M <- temp_M %>%
      filter(cell == cell_list[i])
    
    #Value at which the coral collapses (leaves high coral state) and recovers (leaves low coral state)
    c_collapse <- min(temp_cell_C$dV[temp_cell_C$end_C > high_val])
    c_recover <- max(temp_cell_C$dV[temp_cell_C$end_C < low_val])
    
    #Value at which the macroalgae collapses (leaves high macro state) and recovers (leaves low macro state)
    m_collapse <- min(temp_cell_M$dV[temp_cell_M$end_C > high_val])
    m_recover <- max(temp_cell_M$dV[temp_cell_M$end_C < low_val])
    
    #Saves data in the output dataframe
    ind_output$cell_id[i] <- cell_list[i]
    ind_output$C_collapse[i] <- c_collapse
    ind_output$C_recover[i] <- c_recover
    ind_output$M_collapse[i] <- m_collapse
    ind_output$M_recover[i] <- m_recover
    ind_output$C_range[i] <- c_collapse - c_recover
    ind_output$M_range[i] <- m_collapse - m_recover
    ind_output$collapse_diff[i] <- m_collapse - c_collapse
    ind_output$recover_diff[i] <- m_recover - c_recover
    
  }
  
  #Sets infinite values as NA
  for(i in 1:ncol(ind_output)){
    
    ind_output[,i][is.infinite(ind_output[,i])] <- NA
    
  }
  
  return(ind_output)
  
}

calculate_transition_speed_cluster <- function(dat, rep_id, tail_prop = 0.1){
  #Calculates metrics of transition abruptness on cluster-level data
  #dat is output from the model at the cluster level
  #rep_id is the unique ID of the replicate
  #tail_prop is the range you consider to be "transitioning"
  #e.g., tail_prop = 0.1 means values between 10 and 90% of equilibrium are considered "intermediate" and thus part of a gradual transition
  
  temp <- dat %>%
    filter(unique_id_rep == rep_id)
  temp_C <- dat %>%
    filter(unique_id_rep == rep_id) %>%
    filter(high_initial_var == "C")
  temp_M <- dat %>%
    filter(unique_id_rep == rep_id) %>%
    filter(high_initial_var == "Mi")
  
  high_val <- temp$high_C_init_high_C[1] * (1-tail_prop)
  low_val <- temp$high_C_init_high_C[1] * tail_prop
  
  cluster_list <- unique(temp$cluster_id)
  
  ind_output <- suppressWarnings(data.frame("unique_id_rep" = rep(rep_id, times = length(cluster_list)),
                                            "cluster_id" = integer(length = length(cluster_list)),
                                            "cluster_size" = integer(length = length(cluster_list)),
                                            "high_val" = rep(high_val, times = length(cluster_list)),
                                            "low_val" = rep(low_val, times = length(cluster_list)),
                                            "C_collapse" = numeric(length = length(cluster_list)),
                                            "C_recover" = numeric(length = length(cluster_list)),
                                            "M_collapse" = numeric(length = length(cluster_list)),
                                            "M_recover" = numeric(length = length(cluster_list)),
                                            "C_range" = numeric(length = length(cluster_list)),
                                            "M_range" = numeric(length = length(cluster_list)),
                                            "collapse_diff" = numeric(length = length(cluster_list)),
                                            "recover_diff" = numeric(length = length(cluster_list))))
  
  for(i in 1:length(cluster_list)){
    
    temp_cluster_C <- temp_C %>%
      filter(cluster_id == cluster_list[i])
    
    temp_cluster_M <- temp_M %>%
      filter(cluster_id == cluster_list[i])
    
    c_collapse <- min(temp_cluster_C$dV[temp_cluster_C$end_C > high_val])
    c_recover <- max(temp_cluster_C$dV[temp_cluster_C$end_C < low_val])
    
    m_collapse <- min(temp_cluster_M$dV[temp_cluster_M$end_C > high_val])
    m_recover <- max(temp_cluster_M$dV[temp_cluster_M$end_C < low_val])
    
    ind_output$cluster_id[i] <- cluster_list[i]
    ind_output$cluster_size[i] <- first(temp_cluster_C$cluster_size)
    ind_output$C_collapse[i] <- c_collapse
    ind_output$C_recover[i] <- c_recover
    ind_output$M_collapse[i] <- m_collapse
    ind_output$M_recover[i] <- m_recover
    ind_output$C_range[i] <- c_collapse - c_recover
    ind_output$M_range[i] <- m_collapse - m_recover
    ind_output$collapse_diff[i] <- m_collapse - c_collapse
    ind_output$recover_diff[i] <- m_recover - c_recover
    
  }
  
  for(i in 1:ncol(ind_output)){
    
    ind_output[,i][is.infinite(ind_output[,i])] <- NA
    
  }
  
  return(ind_output)
  
}

calculate_bistable_region <- function(dat, rep_id, cutoff = 0.01){
  #Calculates characteristics of the bistable region for a given replicate landscape
  #Determines the minimum, maximum, and range of herbivory values that the equilibrium coral cover differs (by at least as much as "cutoff") between initial conditions
  #Within this range, it then calculates the mean, median, and maximum difference in equilibrium coral cover between the initial conditions
  #dat is the hysteresis_outputs
  #rep_id is the identifier for the individual rep
  #cutoff is how far the initial conditions are to be considered as "bistable"
  
  #Filters out rep and C initially high
  temp_C <- dat %>%
    filter(unique_id_rep == rep_id) %>%
    filter(high_initial_var == "C") %>%
    arrange(dV)
  
  #Filters out rep and M initially high
  temp_M <- dat %>%
    filter(unique_id_rep == rep_id) %>%
    filter(high_initial_var == "Mi") %>%
    arrange(dV)
  
  #Does the calculations to quantify the bistable range
  bistable_range <- suppressWarnings(data.frame("unique_id_rep" = rep_id,
                                                "p_struc" = dat$p_struc[dat$unique_id_rep==rep_id][1],
                                                "clump" = dat$clump[dat$unique_id_rep==rep_id][1],
                                                "bistable_min" = min(temp_C$dV[temp_C$end_C - temp_M$end_C > cutoff]),
                                                "bistable_max" = max(temp_C$dV[temp_C$end_C - temp_M$end_C > cutoff]),
                                                "bistable_range" = max(temp_C$dV[temp_C$end_C - temp_M$end_C > cutoff]) - 
                                                  min(temp_C$dV[temp_C$end_C - temp_M$end_C > cutoff]),
                                                "max_diff" = max(temp_C$end_C - temp_M$end_C),
                                                "mean_diff" = mean((temp_C$end_C - temp_M$end_C)[temp_C$end_C - temp_M$end_C > cutoff]),
                                                "median_diff" = median((temp_C$end_C - temp_M$end_C)[temp_C$end_C - temp_M$end_C > cutoff]),
                                                "max_dV" = temp_C$dV[(temp_C$end_C - temp_M$end_C) == max(temp_C$end_C - temp_M$end_C)]))
  
  diff_values <- suppressWarnings(data.frame("unique_id_rep" = rep(rep_id, times = nrow(temp_C)),
                                             "p_struc" = rep(dat$p_struc[dat$unique_id_rep==rep_id][1], times = nrow(temp_C)),
                                             "clump" = rep(dat$clump[dat$unique_id_rep==rep_id][1], times = nrow(temp_C)),
                                             "dV" = temp_C$dV,
                                             "diffs" = temp_C$end_C - temp_M$end_C))
  
  return(list(bistable_range, diff_values))
  
}
