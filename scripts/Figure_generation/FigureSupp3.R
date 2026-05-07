rm(list=ls())

library(tidyverse)

#Sets up our initial values and herbivory values
init_frame <- expand_grid(c_init = seq(0.1, 0.6, by = 0.05),
                          m_init = seq(0.1, 0.6, by = 0.05)) %>%
  filter(c_init + m_init <= 0.9)
dv_list <- c(0.05, 2, 4, 6)

################
#LEFT COLUMN -- LANDSCAPE 0_0.1_4
################
time_series <- read.csv("data\\stability_figure_plotting.csv") %>%
  filter(unique_rep_id == "0_0.1_4") %>%
  mutate("init" = paste(c_init, m_init, sep = "_"))

#Goes through each herbivory value to make the four panels
for (j in 1:length(dv_list)){
  
  #Creates plot with initial values
  dev.new()
  plot(init_frame$m_init, init_frame$c_init, pch = 16, ylim = c(0,0.7), xlim = c(0,0.7), cex = 0.5,
       xlab = "Macroalgae cover", ylab = "Coral cover", main = paste("Herbivory = ", dv_list[j], sep = ""), cex.lab = 1.4, cex.axis = 1.3)

  #Plots the time series lines
  for(i in 1:length(unique(time_series$init))){
    
    temp <- time_series %>%
      filter(dv == dv_list[j],
             init == unique(time_series$init)[i])
    coords_x <- c(init_frame$m_init[i], temp$mean_m)
    coords_y <- c(init_frame$c_init[i], temp$mean_c)
    lines(coords_x, coords_y, col = 'gray80', lwd = 0.85)

  }
  #Replots the initial conditions so they are over the lines
  points(init_frame$m_init, init_frame$c_init, pch = 16, cex = 0.5)
  
  #Plots all the equilibria from all initial conditions
  for(i in 1:length(unique(time_series$init))){
    temp <- time_series %>%
      filter(dv == dv_list[j],
             init == unique(time_series$init)[i])
    points(temp$mean_m[nrow(temp)],
           temp$mean_c[nrow(temp)],
           pch = 16, col = 'red', cex = 2.5)
    
  }
  
  #Plots the trajectories and equilibria of the initial conditions used in other analyses
  #Start with high coral low macro (0.6 coral, 0.1 mi)
  temp <- time_series %>%
    filter(dv == dv_list[j],
           init == "0.6_0.1")
  coords_x <- c(temp$mean_m[1], temp$mean_m)
  coords_y <- c(temp$mean_c[1], temp$mean_c)
  lines(coords_x, coords_y, col = 'deepskyblue3', lwd = 2)
  
  #Now low coral high macro
  temp <- time_series %>%
    filter(dv == dv_list[j],
           init == "0.1_0.6")
  coords_x <- c(temp$mean_m[1], temp$mean_m)
  coords_y <- c(temp$mean_c[1], temp$mean_c)
  lines(coords_x, coords_y, col = 'deepskyblue3', lwd = 2)
  
  #Now plot the equilibria
  temp <- time_series %>%
    filter(dv == dv_list[j],
           init == "0.6_0.1")
  points(temp$mean_m[nrow(temp)],
         temp$mean_c[nrow(temp)],
         pch = 16, col = 'deepskyblue4', cex = 2.5)
  
  temp <- time_series %>%
    filter(dv == dv_list[j],
           init == "0.1_0.6")
  points(temp$mean_m[nrow(temp)],
         temp$mean_c[nrow(temp)],
         pch = 16, col = 'deepskyblue4', cex = 2.5)

}

####################
#RIGHT COLUMN -- LANDSCAPE 0_0.8_5
####################

time_series <- read.csv("data\\stability_figure_plotting.csv") %>%
  filter(unique_rep_id == "0_0.8_5") %>%
  mutate("init" = paste(c_init, m_init, sep = "_"))

#Goes through each herbivory value to make the four panels
for (j in 1:length(dv_list)){
  
  #Creates plot with initial values
  dev.new()
  plot(init_frame$m_init, init_frame$c_init, pch = 16, ylim = c(0,0.7), xlim = c(0,0.7), cex = 0.5,
       xlab = "Macroalgae cover", ylab = "Coral cover", main = paste("Herbivory = ", dv_list[j], sep = ""), cex.lab = 1.4, cex.axis = 1.3)
  
  #Plots the time series lines
  for(i in 1:length(unique(time_series$init))){
    
    temp <- time_series %>%
      filter(dv == dv_list[j],
             init == unique(time_series$init)[i])
    coords_x <- c(init_frame$m_init[i], temp$mean_m)
    coords_y <- c(init_frame$c_init[i], temp$mean_c)
    lines(coords_x, coords_y, col = 'gray80', lwd = 0.85)
    
  }
  #Replots the initial conditions so they are over the lines
  points(init_frame$m_init, init_frame$c_init, pch = 16, cex = 0.5)
  
  #Plots all the equilibria from all initial conditions
  for(i in 1:length(unique(time_series$init))){
    temp <- time_series %>%
      filter(dv == dv_list[j],
             init == unique(time_series$init)[i])
    points(temp$mean_m[nrow(temp)],
           temp$mean_c[nrow(temp)],
           pch = 16, col = 'red', cex = 2.5)
    
  }
  
  #Plots the trajectories and equilibria of the initial conditions used in other analyses
  #Start with high coral low macro (0.6 coral, 0.1 mi)
  temp <- time_series %>%
    filter(dv == dv_list[j],
           init == "0.6_0.1")
  coords_x <- c(temp$mean_m[1], temp$mean_m)
  coords_y <- c(temp$mean_c[1], temp$mean_c)
  lines(coords_x, coords_y, col = 'deepskyblue3', lwd = 2)
  
  #Now low coral high macro
  temp <- time_series %>%
    filter(dv == dv_list[j],
           init == "0.1_0.6")
  coords_x <- c(temp$mean_m[1], temp$mean_m)
  coords_y <- c(temp$mean_c[1], temp$mean_c)
  lines(coords_x, coords_y, col = 'deepskyblue3', lwd = 2)
  
  #Now plot the equilibria
  temp <- time_series %>%
    filter(dv == dv_list[j],
           init == "0.6_0.1")
  points(temp$mean_m[nrow(temp)],
         temp$mean_c[nrow(temp)],
         pch = 16, col = 'deepskyblue4', cex = 2.5)
  
  temp <- time_series %>%
    filter(dv == dv_list[j],
           init == "0.1_0.6")
  points(temp$mean_m[nrow(temp)],
         temp$mean_c[nrow(temp)],
         pch = 16, col = 'deepskyblue4', cex = 2.5)
  
}
