# Coeliac disease Markov model
# Edna Keeney
library(tictoc)
tic()
rm(list=ls())
set.seed(14143)
  
  # Define the number and names of treatments

  n_treatments <- 3
  t_names <- c("Test", "Test + biopsy", "Double test")
  
  # Define the number and names of states of the model
  n_states <- 13
  state_names <- c("CD GFD no complications",
                  "CD GFD subfertility",
                  "CD GFD osteoporosis",
                  "CD GFD NHL",
                  "CD no GFD no complications",
                  "CD no GFD subfertility",
                  "CD no GFD osteoporosis",
                  "CD no GFD NHL",
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
  
  #prevalence of coeliac disease
  p_cd <- 0.5 
  
  starting_age <- 30 #Max is 50 with 50 cycles
  starting_age_columnandrow <- read.csv("starting_age_column.csv")
  starting_age_columnandrow$starting_age_row <- c(1:11)
  starting_age_column <- starting_age_columnandrow[starting_age_columnandrow$Starting.age == starting_age, 2]
  starting_age_row <- starting_age_columnandrow[starting_age_columnandrow$Starting.age == starting_age, 2]
  
  perspective <- "NHS" #Other option is "NHS+OOP" for when out-of-pocket costs for iron supplements and gluten free products are included
  
  source("generate_cohort_vectors.R")
  source("generate_state_costs.R")
  source("generate_state_qalys.R")
  source("generate_model_parameters.R")
  source("generate_transition_matrices.R")
  source("generate_net_benefit.R")
  #############################################################################
  ## Input parameters #########################################################
  #############################################################################
  
  input_parameters <- generate_model_parameters(starting_age)
  attach(input_parameters)
  # Osteoporosis probabilities On GFD
 
  osteoporosis_probability_GFD_all <- data.frame(osteoporosis_probability_GFD_0, osteoporosis_probability_GFD_10, osteoporosis_probability_GFD_20, osteoporosis_probability_GFD_30,
                                                 osteoporosis_probability_GFD_40, osteoporosis_probability_GFD_50, osteoporosis_probability_GFD_60,
                                                 osteoporosis_probability_GFD_70, osteoporosis_probability_GFD_80, osteoporosis_probability_GFD_90)
  
  # Subfertility probabilities on GFD
 
  subfertility_probability_GFD_all <- data.frame(subfertility_probability_GFD_0, subfertility_probability_GFD_10, subfertility_probability_GFD_20, subfertility_probability_GFD_30,
                                                 subfertility_probability_GFD_40, subfertility_probability_GFD_50, subfertility_probability_GFD_60,
                                                 subfertility_probability_GFD_70, subfertility_probability_GFD_80, subfertility_probability_GFD_90)
  
  # Corresponding probabilities not on GFD

  osteoporosis_probability_noGFD_all <- data.frame(osteoporosis_probability_noGFD_0, osteoporosis_probability_noGFD_10, osteoporosis_probability_noGFD_20, osteoporosis_probability_noGFD_30,
                                                   osteoporosis_probability_noGFD_40, osteoporosis_probability_noGFD_50, osteoporosis_probability_noGFD_60,
                                                   osteoporosis_probability_noGFD_70, osteoporosis_probability_noGFD_80, osteoporosis_probability_noGFD_90)
  
  
  subfertility_probability_noGFD_all <- data.frame(subfertility_probability_noGFD_0, subfertility_probability_noGFD_10, subfertility_probability_noGFD_20, subfertility_probability_noGFD_30,
                                                   subfertility_probability_noGFD_40, subfertility_probability_noGFD_50, subfertility_probability_noGFD_60,
                                                   subfertility_probability_noGFD_70, subfertility_probability_noGFD_80, subfertility_probability_noGFD_90)
  lifetables <- read.csv("lifetables.csv")
  percentage_male <- 0.5
  lifetables$Overall <- (percentage_male * lifetables$Males) + ((1-percentage_male) * lifetables$Females)
  death_probability_nocomplications	<- data.frame(lifetables$Age, lifetables$Overall)
  
  death_probability_osteoporosis <-	lifetables
  death_probability_osteoporosis$Males <- death_probability_osteoporosis$Males*3.5   #Currently  not probabilistic as lifetables are not probabilistic
  death_probability_osteoporosis$Females <- death_probability_osteoporosis$Females*2.4
  death_probability_osteoporosis$Overall <- (percentage_male * death_probability_osteoporosis$Males) + ((1-percentage_male) * death_probability_osteoporosis$Females)
  death_probability_osteoporosis	<- data.frame(death_probability_osteoporosis$Age, death_probability_osteoporosis$Overall)
  
   # HT: The code that follows should go in a a generate_transition_matrices() function which
   # takes input_parameters as input and is called within the (main model running) generate_net_benefit() function
  
  transition_matrices <- generate_transition_matrices(input_parameters)
   #the issue here is that 'subfertility_probability_GFD_all' is it's own dataframe that is used when generating the transition matrices. When I generate the 
  # data frame of input parameters however it's columns are just merged with the larger dataframe when it needs to stay as it's own dataframe. Same issue with osteoporosis and death

    # HT: Obviously remove this check, or include it in the generate_transition_matrices() function
    # and have it output a warning() if the sums are not 1.
    rowSums (transition_matrices[1, 43, , ], na.rm = FALSE, dims = 1)
    
    
   
    # HT: What follows will go in a generate_state_qalys() function
    # which takes input_parameters as input and is called from within generate_net_benefit()
    eq5d_norms <- read.csv("eq5d_norms.csv")
    eq5d_norms$age <- c(0, 10, 20, 30, 40, 50, 60, 70, 80)
    
    state_qalys <- generate_state_qalys(input_parameters)

    probability_IDA <- data.frame(probability_IDA_0 , probability_IDA_10 , probability_IDA_20 , probability_IDA_30 , probability_IDA_40
                                  ,probability_IDA_50 , probability_IDA_60 , probability_IDA_70 , probability_IDA_80 , probability_IDA_90)
    
    state_costs <- generate_state_costs(input_parameters)
    #same issue here with probability_IDA
    
  
#Need to include cost of GFD - Penny contacting Coeliac UK
  
  # Define the treatment costs
  # One for each PSA sample and each treatment

  treatment_costs<-array(dim=c(n_treatments,n_samples),dimnames=list(t_names,NULL))
  
  treatment_costs["Test", ] <- 10
  probability_biopsy <- runif(n = n_samples, min = 0.6, max = 0.8)
  treatment_costs["Test + biopsy", ] <- 10 + (probability_biopsy * 530)
  treatment_costs["Double test", ] <- 20
  
  #############################################################################
  ## Accuracy of tests ########################################################
  #############################################################################

  
  tp <- fn <- fp <- tn <- matrix(nrow=n_samples, ncol=n_treatments)
  
  #IgA EMA sensitivity in adults: 88.0 (75.2, 94.7)
  sens_test <- 0.88
  sens_se <- (0.947 - 0.752)/3.92
  sens_alpha <- (sens_test ^ 2 * (1-sens_test)/sens_se ^ 2) - sens_test
  sens_beta <- (sens_alpha / sens_test) - sens_alpha
  #IgA EMA specificity in adults: 99.6 (92.3, 100.0)
  spec_test <- 0.996
  spec_se <- (1 - 0.923)/3.92 
  spec_alpha <- (spec_test ^ 2 * (1 - spec_test)/spec_se ^ 2)- spec_test
  spec_beta <- (spec_alpha / spec_test) - spec_alpha
  
  sens_testbiopsy <- 1
  spec_testbiopsy <- 1
  sens_doubletest <- 1
  spec_doubletest <- 1
  
  # Probabilities for test 
  tp[,1] <- (n_samples * p_cd * rbeta(n=n_samples, shape1 = sens_alpha, shape2 = sens_beta))/n_samples
  fn[,1] <- 1 - tp[,1]  
  tn[,1] <- (n_samples * p_cd * rbeta(n=n_samples, shape1 = spec_alpha, shape2 = spec_beta))/n_samples
  fp[,1] <- 1 - tn[,1]
  
  # Probabilities for test + biopsy
  tp[,2] <- (n_samples * p_cd * sens_testbiopsy)/n_samples
  fn[,2] <- 1 - tp[,2] 
  tn[,2] <- (n_samples * p_cd * spec_testbiopsy)/n_samples
  fp[,2] <- 1 - tn[,2]
  
  # Probabilities for Double test
  tp[,3] <- (n_samples * p_cd * sens_doubletest)/n_samples
  fn[,3] <- 1 - tp[,3] 
  tn[,3] <- (n_samples * p_cd * spec_doubletest)/n_samples
  fp[,3] <- 1 - tn[,3]
  
  
  
  #############################################################################
  ## Simulation ###############################################################
  #############################################################################
  
  cohort_vectors <- generate_cohort_vectors(input_parameters)
 
  rowSums (cohort_vectors[, 2, , ], na.rm = FALSE, dims = 1)
 
   # Build an array to store the costs and QALYs accrued per cycle
  # One for each treatment, for each PSA sample, for each cycle
  # These will be filled in below in the main model code
  # Then discounted and summed to contribute to total costs and total QALYs
  cycle_costs <- array(dim = c(n_treatments, n_samples, n_cycles),
                     dimnames = list(t_names, NULL, NULL))
  cycle_qalys <- array(dim = c(n_treatments, n_samples, n_cycles),
                     dimnames = list(t_names, NULL, NULL))
  
  # Build arrays to store the total costs and total QALYs
  # There is one for each treatment and each PSA sample
  # These are filled in below using cycle_costs, 
  # treatment_costs, and cycle_qalys
  total_costs <- array(dim = c(n_treatments, n_samples),
                     dimnames = list(t_names, NULL))
  total_qalys <- array(dim = c(n_treatments, n_samples),
                     dimnames = list(t_names, NULL))
  
  disc_vec <- (1/1.035) ^ rep(c(0:(n_cycles/2-1)), each = 2)

  
  # Now use the BCEA package to analyse the results___
  output <- generate_net_benefit(transition_matrices, state_costs, state_qalys, cohort_vectors)
  output
  
  pkgs <- c("MASS","Rtools","devtools")
  repos <- c("https://cran.rstudio.com", "https://www.math.ntnu.no/inla/R/stable") 
  install.packages(pkgs,repos=repos,dependencies = "Depends")
  devtools::install_github("giabaio/BCEA")
   library(BCEA)
  m <- bcea(e = t(total_qalys), c = t(total_costs), ref = 1, interventions = t_names)
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

