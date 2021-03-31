# Coeliac disease Markov model
# Edna Keeney
library(tictoc)
library(BCEA)
library(SimDesign)
library(BCEA)



tic()
rm(list=ls())
set.seed(14143)
  
  # Define the number and names of treatments

  n_treatments <- 3
  t_names <- c("IgAEMA", "IgATTGplusIgAEMA", "Double test")
  
  # Define the number and names of states of the model
  n_states <- 9
  state_names <- c("CD GFD no complications",
                  "CD GFD subfertility",
                  "CD GFD osteoporosis",
                  "CD GFD NHL",
                  "Undiagnosed CD no complications",
                  "Undiagnosed CD subfertility",
                  "Undiagnosed CD osteoporosis",
                  "Undiagnosed CD NHL",
                  "Death")
 
  
  
   # Define the number of cycles
  n_cycles <- 50
  
  # Define simulation parameters
  # This is the number of PSA samples to use
  n_samples <- 100
  
  #pre-test probability of coeliac disease 
  p_cd <- 0.15 
  
  starting_age <- 30 #Max is 50 with 50 cycles
 
  
  perspective <- "NHS+OOP" #Options are "NHS" or "NHS+OOP" if out-of-pocket costs for iron supplements and gluten free products are to be included
  
  source("generate_state_costs.R")
  source("generate_state_qalys.R")
  source("generate_model_parameters.R")
  source("generate_transition_matrices.R")
  source("generate_net_benefit.R")
  
  
  #generate input parameters
  input_parameters <- generate_model_parameters(starting_age)
  

  #generate results
  output <- generate_net_benefit(input_parameters)
  output
  
  
  # Now use the BCEA package to analyse the results
 # pkgs <- c("MASS","Rtools","devtools")
  #repos <- c("https://cran.rstudio.com", "https://www.math.ntnu.no/inla/R/stable") 
  #install.packages(pkgs,repos=repos,dependencies = "Depends")
  #devtools::install_github("giabaio/BCEA")
  
  m <- bcea(e = t(output$total_qalys), c = t(output$total_costs), ref = 1, interventions = t_names)
summary(m)
eib.plot(m, comparison = NULL, pos =
           c(1, 0), size = NULL, plot.cri = NULL, graph
         = c("base", "ggplot2", "plotly"))
evi.plot(m, graph = c("base", "ggplot2",
                       "plotly"))
ceac_plot(m, comparison = NULL,
          pos = c(1, 0), graph = c("base",
                                   "ggplot2", "plotly"))

ceplane_plot(m, comparison =
               NULL, pos = c(1, 0), graph = c("base",
                                              "ggplot2", "plotly"))
sim.table(m)
toc()

