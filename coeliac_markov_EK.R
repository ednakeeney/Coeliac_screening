# Coeliac disease Markov model
# Edna Keeney

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
  starting_age_column <- read.csv("starting_age_column.csv")
  starting_age_column <- starting_age_column[starting_age_column$Starting.age == starting_age, 2]
  #############################################################################
  ## Input parameters #########################################################
  #############################################################################
  
  adherence	<- 0.75      #will come from targeted review
  duration_of_symptoms <- rnorm(n = n_samples, mean = 10.93, sd = 14.88)     #calculated from Violato, what about SD? Should some of these be negative?
  rate_of_symptoms <- 1/duration_of_symptoms
  #probability_late_diagnosis <- 1 - exp(-rate_of_symptoms)
  probability_late_diagnosis <- 0.06

    prevalence <- read.csv("CPRD prevalence.csv")
  
  #Initial cohort at diagnosis - depends on age at diagnosis
 probability_subfertility <- rbeta(n=n_samples, shape1 = prevalence[prevalence$Age.categories == starting_age, 5], shape2 = prevalence[prevalence$Age.categories == starting_age, 2] - prevalence[prevalence$Age.categories == starting_age, 5])
 probability_osteoporosis <- rbeta(n=n_samples, shape1 = prevalence[prevalence$Age.categories == starting_age, 4], shape2 = prevalence[prevalence$Age.categories == starting_age, 2] - prevalence[prevalence$Age.categories == starting_age, 4])
 probability_NHL <- rbeta(n=n_samples, shape1 = prevalence[prevalence$Age.categories == starting_age, 3], shape2 = prevalence[prevalence$Age.categories == starting_age, 2] - prevalence[prevalence$Age.categories == starting_age, 3])
 probability_nocomplications <- rbeta(n=n_samples, shape1 = prevalence[prevalence$Age.categories == starting_age, 7], shape2 = prevalence[prevalence$Age.categories == starting_age, 2] - prevalence[prevalence$Age.categories == starting_age, 7])
  
  #On GFD
 osteoporosis_probability <- read.csv("osteoporosis_rate_nice.csv")
 osteoporosis_probability_GFD_0	<- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_alpha[1]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_beta[1]))
 osteoporosis_probability_GFD_10	<- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_alpha[2]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_beta[2]))
 osteoporosis_probability_GFD_20	<- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_alpha[3]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_beta[3]))
 osteoporosis_probability_GFD_30	<- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_alpha[4]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_beta[4]))
  osteoporosis_probability_GFD_40	<- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_alpha[5]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_beta[5]))
  osteoporosis_probability_GFD_50 <- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_alpha[6]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_beta[6]))
  osteoporosis_probability_GFD_60	<- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_alpha[7]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_beta[7]))
  osteoporosis_probability_GFD_70	<- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_alpha[8]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_beta[8]))
  osteoporosis_probability_GFD_80 <- 	rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_alpha[9]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_beta[9]))
  osteoporosis_probability_GFD_90 <- 	rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_alpha[10]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_beta[10]))
 osteoporosis_probability_GFD_all <- data.frame(osteoporosis_probability_GFD_0, osteoporosis_probability_GFD_10, osteoporosis_probability_GFD_20, osteoporosis_probability_GFD_30,
                                            osteoporosis_probability_GFD_40, osteoporosis_probability_GFD_50, osteoporosis_probability_GFD_60,
                                            osteoporosis_probability_GFD_70, osteoporosis_probability_GFD_80, osteoporosis_probability_GFD_90)

 
 subfertility_probability <- read.csv("subfertility_rate_nice.csv")
 subfertility_probability_GFD_0	<- rbeta(n=n_samples, shape1 = as.numeric(subfertility_probability$subfertility_GFD_alpha[1]), shape2 = as.numeric(subfertility_probability$subfertility_GFD_beta[1]))
 subfertility_probability_GFD_10	<- rbeta(n=n_samples, shape1 = as.numeric(subfertility_probability$subfertility_GFD_alpha[2]), shape2 = as.numeric(subfertility_probability$subfertility_GFD_beta[2]))
 subfertility_probability_GFD_20	<- rbeta(n=n_samples, shape1 = as.numeric(subfertility_probability$subfertility_GFD_alpha[3]), shape2 = as.numeric(subfertility_probability$subfertility_GFD_beta[3]))
  subfertility_probability_GFD_30 <-	rbeta(n=n_samples, shape1 = as.numeric(subfertility_probability$subfertility_GFD_alpha[4]), shape2 = as.numeric(subfertility_probability$subfertility_GFD_beta[4]))
  subfertility_probability_GFD_40 <- rbeta(n=n_samples, shape1 = as.numeric(subfertility_probability$subfertility_GFD_alpha[5]), shape2 = as.numeric(subfertility_probability$subfertility_GFD_beta[5]))
  subfertility_probability_GFD_50 <- 0
  subfertility_probability_GFD_60	<- 0
  subfertility_probability_GFD_70	<- 0
  subfertility_probability_GFD_80 <- 0
  subfertility_probability_GFD_90 <- 0 
  subfertility_probability_GFD_all <- data.frame(subfertility_probability_GFD_0, subfertility_probability_GFD_10, subfertility_probability_GFD_20, subfertility_probability_GFD_30,
                                                 subfertility_probability_GFD_40, subfertility_probability_GFD_50, subfertility_probability_GFD_60,
                                                 subfertility_probability_GFD_70, subfertility_probability_GFD_80, subfertility_probability_GFD_90)
  
  
  NHL_probability <- read.csv("NHL_rate_nice.csv")
  NHL_probability_GFD_18orless	<- rbeta(n=n_samples, shape1 = as.numeric(NHL_probability$NHL_GFD_alpha[1]), shape2 = as.numeric(NHL_probability$NHL_GFD_beta[1]))
  NHL_probability_GFD_18plus	<- rbeta(n=n_samples, shape1 = as.numeric(NHL_probability$NHL_GFD_alpha[2]), shape2 = as.numeric(NHL_probability$NHL_GFD_beta[2]))
  NHL_probability_GFD <- NHL_probability_GFD_18plus
  
  #Not on GFD
  osteoporosis_probability_noGFD_0	<- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_alpha[1]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_beta[1]))
  osteoporosis_probability_noGFD_10	<- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_alpha[2]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_beta[2]))
  osteoporosis_probability_noGFD_20	<- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_alpha[3]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_beta[3]))
  osteoporosis_probability_noGFD_30	<- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_alpha[4]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_beta[4]))
  osteoporosis_probability_noGFD_40	<- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_alpha[5]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_beta[5]))
  osteoporosis_probability_noGFD_50 <- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_alpha[6]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_beta[6]))
  osteoporosis_probability_noGFD_60	<- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_alpha[7]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_beta[7]))
  osteoporosis_probability_noGFD_70	<- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_alpha[8]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_beta[8]))
  osteoporosis_probability_noGFD_80 <- 	rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_alpha[9]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_beta[9]))
  osteoporosis_probability_noGFD_90 <- 	rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_alpha[10]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_beta[10]))
  osteoporosis_probability_noGFD_all <- data.frame(osteoporosis_probability_noGFD_0, osteoporosis_probability_noGFD_10, osteoporosis_probability_noGFD_20, osteoporosis_probability_noGFD_30,
                                                 osteoporosis_probability_noGFD_40, osteoporosis_probability_noGFD_50, osteoporosis_probability_noGFD_60,
                                                 osteoporosis_probability_noGFD_70, osteoporosis_probability_noGFD_80, osteoporosis_probability_noGFD_90)
 
  
  
  subfertility_probability_noGFD_0	<- rbeta(n=n_samples, shape1 = as.numeric(subfertility_probability$subfertility_noGFD_alpha[1]), shape2 = as.numeric(subfertility_probability$subfertility_noGFD_beta[1]))
  subfertility_probability_noGFD_10	<- rbeta(n=n_samples, shape1 = as.numeric(subfertility_probability$subfertility_noGFD_alpha[2]), shape2 = as.numeric(subfertility_probability$subfertility_noGFD_beta[2]))
  subfertility_probability_noGFD_20	<- rbeta(n=n_samples, shape1 = as.numeric(subfertility_probability$subfertility_noGFD_alpha[3]), shape2 = as.numeric(subfertility_probability$subfertility_noGFD_beta[3]))
  subfertility_probability_noGFD_30 <-	rbeta(n=n_samples, shape1 = as.numeric(subfertility_probability$subfertility_noGFD_alpha[4]), shape2 = as.numeric(subfertility_probability$subfertility_noGFD_beta[4]))
  subfertility_probability_noGFD_40 <- rbeta(n=n_samples, shape1 = as.numeric(subfertility_probability$subfertility_noGFD_alpha[5]), shape2 = as.numeric(subfertility_probability$subfertility_noGFD_beta[5]))
  subfertility_probability_noGFD_50 <- 0
  subfertility_probability_noGFD_60	<- 0
  subfertility_probability_noGFD_70	<- 0
  subfertility_probability_noGFD_80 <- 0
  subfertility_probability_noGFD_90 <- 0 
  subfertility_probability_noGFD_all <- data.frame(subfertility_probability_noGFD_0, subfertility_probability_noGFD_10, subfertility_probability_noGFD_20, subfertility_probability_noGFD_30,
                                                 subfertility_probability_noGFD_40, subfertility_probability_noGFD_50, subfertility_probability_noGFD_60,
                                                 subfertility_probability_noGFD_70, subfertility_probability_noGFD_80, subfertility_probability_noGFD_90)
 
  
  NHL_probability_noGFD_18orless	<- rbeta(n=n_samples, shape1 = as.numeric(NHL_probability$NHL_noGFD_alpha[1]), shape2 = as.numeric(NHL_probability$NHL_noGFD_beta[1]))
  NHL_probability_noGFD_18plus	<- rbeta(n=n_samples, shape1 = as.numeric(NHL_probability$NHL_noGFD_alpha[2]), shape2 = as.numeric(NHL_probability$NHL_noGFD_beta[2]))
  
  lifetables <- read.csv("lifetables.csv")
  percentage_male <- 0.5
  lifetables$Overall <- (percentage_male * lifetables$Males) + ((1-percentage_male) * lifetables$Females)
 # lifetables30plus <- subset(lifetables, lifetables$Age > 29)
  death_probability_nocomplications	<- data.frame(lifetables$Age, lifetables$Overall)
 
  mortality_NHL <- read.csv("NHL mortality.csv")
  
   death_probability_NHL <-	0.006 #will relate to table above
 
   #During the total study period, the hazard ratio for one-year all-cause mortality was 3.5 times (95% CI: 3.28–3.74) 
   #greater for male hip fracture patients than control subjects and 2.4 times (95% CI: 2.31–2.50) greater than 
   #controls for females. 
    
   death_probability_osteoporosis <-	lifetables
   death_probability_osteoporosis$Males <- death_probability_osteoporosis$Males*3.5   #Currently  not probabilistic as lifetables are not probabilistic
   death_probability_osteoporosis$Females <- death_probability_osteoporosis$Females*2.4
   death_probability_osteoporosis$Overall <- (percentage_male * death_probability_osteoporosis$Males) + ((1-percentage_male) * death_probability_osteoporosis$Females)
   death_probability_osteoporosis	<- data.frame(death_probability_osteoporosis$Age, death_probability_osteoporosis$Overall)
  
   
   # There is one transition matrix for each cycle and each PSA sample
  # Store them in an array with (before filling in below) NA entries
  transition_matrices <- array(dim=c(n_samples, n_cycles, n_states, n_states),
                             dimnames=list(NULL, NULL, state_names, state_names))

  
  # Now define the QALYS associated with the states per cycle
  # There is one for each PSA sample and each state
  # Store in an NA array and then fill in below
  state_qalys <- array(dim=c(n_samples, n_states), dimnames = list(NULL, state_names))
  
  # And finally define the state costs
  # There is one for each PSA sample and each state
  # Store in an NA array and then fill in below
  state_costs <- array(dim=c(n_samples, n_states), dimnames = list(NULL, state_names))
  
  
  # Define transition matrices, state utilities and costs 

  #CD GFD 
    transition_matrices[, c(1:10), "CD GFD no complications", "CD GFD subfertility"] <- subfertility_probability_GFD_all[, starting_age_column]
    transition_matrices[, c(11:20), "CD GFD no complications", "CD GFD subfertility"] <- subfertility_probability_GFD_all[, starting_age_column + 1]
    transition_matrices[, c(21:30), "CD GFD no complications", "CD GFD subfertility"] <- subfertility_probability_GFD_all[, starting_age_column + 2]
    transition_matrices[, c(31:40), "CD GFD no complications", "CD GFD subfertility"] <- subfertility_probability_GFD_all[, starting_age_column + 3]
    transition_matrices[, c(41:50), "CD GFD no complications", "CD GFD subfertility"] <- subfertility_probability_GFD_all[, starting_age_column + 4]
  
    transition_matrices[, c(1:10), "CD GFD no complications", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_all[, starting_age_column]
    transition_matrices[, c(11:20), "CD GFD no complications", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_all[, starting_age_column + 1]
    transition_matrices[, c(21:30), "CD GFD no complications", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_all[, starting_age_column + 2]
    transition_matrices[, c(31:40), "CD GFD no complications", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_all[, starting_age_column + 3]
    transition_matrices[, c(41:50), "CD GFD no complications", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_all[, starting_age_column + 4]
    transition_matrices[, , "CD GFD no complications", "CD GFD NHL"] <- NHL_probability_GFD
    for (i_age in 1:50){
    transition_matrices[, i_age, "CD GFD no complications", "Death"] <- death_probability_nocomplications[i_age, 2]
    }
   
    transition_matrices[, , "CD GFD no complications", "CD GFD no complications"] <- 1 - transition_matrices[, , "CD GFD no complications", "CD GFD subfertility"] - 
      transition_matrices[, , "CD GFD no complications", "CD GFD osteoporosis"] - transition_matrices[, , "CD GFD no complications", "CD GFD NHL"] -  transition_matrices[ , , "CD GFD no complications", "Death"]


    transition_matrices[, c(1:10),"CD GFD subfertility", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_all[, starting_age_column]
    transition_matrices[, c(11:20), "CD GFD subfertility", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_all[, starting_age_column + 1]
    transition_matrices[, c(21:30), "CD GFD subfertility", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_all[, starting_age_column + 2]
    transition_matrices[, c(31:40), "CD GFD subfertility", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_all[, starting_age_column + 3]
    transition_matrices[, c(41:50), "CD GFD subfertility", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_all[, starting_age_column + 4]
    transition_matrices[, ,"CD GFD subfertility", "CD GFD NHL"] <- NHL_probability_GFD
    for (i_age in 1:50){
      transition_matrices[, i_age, "CD GFD subfertility", "Death"] <- death_probability_nocomplications[i_age, 2]
    }
      transition_matrices[, , "CD GFD subfertility", "CD GFD subfertility"] <- 1 -  transition_matrices[, ,"CD GFD subfertility", "CD GFD osteoporosis"] - 
        transition_matrices[, ,"CD GFD subfertility", "CD GFD NHL"] - transition_matrices[, , "CD GFD subfertility", "Death"]
    
    
    transition_matrices[, ,"CD GFD osteoporosis", "CD GFD NHL"] <- NHL_probability_GFD
    for (i_age in 1:50){
    transition_matrices[, i_age, "CD GFD osteoporosis", "Death"] <- death_probability_osteoporosis[i_age, 2]
    }
      transition_matrices[, , "CD GFD osteoporosis", "CD GFD osteoporosis"] <- 1 - 
        transition_matrices[, ,"CD GFD osteoporosis", "CD GFD NHL"] - transition_matrices[, , "CD GFD osteoporosis", "Death"]
   
    transition_matrices[, , "CD GFD NHL", "Death"] <- death_probability_NHL
    transition_matrices[, ,"CD GFD NHL", "CD GFD NHL"] <- 1 - death_probability_NHL
   

#CD no GFD    
    transition_matrices[, c(1:10), "CD no GFD no complications", "CD no GFD subfertility"] <- subfertility_probability_noGFD_all[, starting_age_column]
    transition_matrices[, c(11:20), "CD no GFD no complications", "CD no GFD subfertility"] <- subfertility_probability_noGFD_all[, starting_age_column + 1]
    transition_matrices[, c(21:30), "CD no GFD no complications", "CD no GFD subfertility"] <- subfertility_probability_noGFD_all[, starting_age_column + 2]
    transition_matrices[, c(31:40), "CD no GFD no complications", "CD no GFD subfertility"] <- subfertility_probability_noGFD_all[, starting_age_column + 3]
    transition_matrices[, c(41:50), "CD no GFD no complications", "CD no GFD subfertility"] <- subfertility_probability_noGFD_all[, starting_age_column + 4]
    
    transition_matrices[, c(1:10), "CD no GFD no complications", "CD no GFD osteoporosis"] <- osteoporosis_probability_noGFD_all[, starting_age_column]
    transition_matrices[, c(11:20), "CD no GFD no complications", "CD no GFD osteoporosis"] <- osteoporosis_probability_noGFD_all[, starting_age_column + 1]
    transition_matrices[, c(21:30), "CD no GFD no complications", "CD no GFD osteoporosis"] <- osteoporosis_probability_noGFD_all[, starting_age_column + 2]
    transition_matrices[, c(31:40), "CD no GFD no complications", "CD no GFD osteoporosis"] <- osteoporosis_probability_noGFD_all[, starting_age_column + 3]
    transition_matrices[, c(41:50), "CD no GFD no complications", "CD no GFD osteoporosis"] <- osteoporosis_probability_noGFD_all[, starting_age_column + 4]
    transition_matrices[, , "CD no GFD no complications", "CD no GFD NHL"] <- NHL_probability_noGFD
    for (i_age in 1:50){
      transition_matrices[, i_age, "CD no GFD no complications", "Death"] <- death_probability_nocomplications[i_age, 2]
    }
    
    transition_matrices[, , "CD no GFD no complications", "CD no GFD no complications"] <- 1 - transition_matrices[, , "CD no GFD no complications", "CD no GFD subfertility"] - 
      transition_matrices[, , "CD no GFD no complications", "CD no GFD osteoporosis"] - transition_matrices[, , "CD no GFD no complications", "CD no GFD NHL"] -  transition_matrices[ , , "CD no GFD no complications", "Death"]
    
    
    transition_matrices[, c(1:10),"CD no GFD subfertility", "CD no GFD osteoporosis"] <- osteoporosis_probability_noGFD_all[, starting_age_column]
    transition_matrices[, c(11:20), "CD no GFD subfertility", "CD no GFD osteoporosis"] <- osteoporosis_probability_noGFD_all[, starting_age_column + 1]
    transition_matrices[, c(21:30), "CD no GFD subfertility", "CD no GFD osteoporosis"] <- osteoporosis_probability_noGFD_all[, starting_age_column + 2]
    transition_matrices[, c(31:40), "CD no GFD subfertility", "CD no GFD osteoporosis"] <- osteoporosis_probability_noGFD_all[, starting_age_column + 3]
    transition_matrices[, c(41:50), "CD no GFD subfertility", "CD no GFD osteoporosis"] <- osteoporosis_probability_noGFD_all[, starting_age_column + 4]
    transition_matrices[, ,"CD no GFD subfertility", "CD no GFD NHL"] <- NHL_probability_noGFD
    for (i_age in 1:50){
      transition_matrices[, i_age, "CD no GFD subfertility", "Death"] <- death_probability_nocomplications[i_age, 2]
    }
    transition_matrices[, , "CD no GFD subfertility", "CD no GFD subfertility"] <- 1 -  transition_matrices[, ,"CD no GFD subfertility", "CD no GFD osteoporosis"] - 
      transition_matrices[, ,"CD no GFD subfertility", "CD no GFD NHL"] - transition_matrices[, , "CD no GFD subfertility", "Death"]
  
    
    transition_matrices[, ,"CD no GFD osteoporosis", "CD no GFD NHL"] <- NHL_probability_noGFD
    for (i_age in 1:50){
      transition_matrices[, i_age, "CD no GFD osteoporosis", "Death"] <- death_probability_osteoporosis[i_age, 2]
    }
    transition_matrices[, , "CD no GFD osteoporosis", "CD no GFD osteoporosis"] <- 1 - 
      transition_matrices[, ,"CD no GFD osteoporosis", "CD no GFD NHL"] - transition_matrices[, , "CD no GFD osteoporosis", "Death"]
  
    
    transition_matrices[, , "CD no GFD NHL", "Death"] <- death_probability_NHL
    transition_matrices[, ,"CD no GFD NHL", "CD no GFD NHL"] <- 1 - death_probability_NHL



#Undiagnosed CD  
    transition_matrices[, ,"Undiagnosed CD no complications", "CD GFD no complications"] <- adherence * probability_late_diagnosis
    transition_matrices[, ,"Undiagnosed CD no complications", "CD no GFD no complications"] <- (1-adherence) * probability_late_diagnosis
   
    transition_matrices[, c(1:10), "Undiagnosed CD no complications", "Undiagnosed CD subfertility"] <- subfertility_probability_noGFD_all[, starting_age_column]
    transition_matrices[, c(11:20), "Undiagnosed CD no complications", "Undiagnosed CD subfertility"] <- subfertility_probability_noGFD_all[, starting_age_column + 1]
    transition_matrices[, c(21:30), "Undiagnosed CD no complications", "Undiagnosed CD subfertility"] <- subfertility_probability_noGFD_all[, starting_age_column + 2]
    transition_matrices[, c(31:40), "Undiagnosed CD no complications", "Undiagnosed CD subfertility"] <- subfertility_probability_noGFD_all[, starting_age_column + 3]
    transition_matrices[, c(41:50), "Undiagnosed CD no complications", "Undiagnosed CD subfertility"] <- subfertility_probability_noGFD_all[, starting_age_column + 4]
    
    transition_matrices[, c(1:10), "Undiagnosed CD no complications", "Undiagnosed CD osteoporosis"] <- osteoporosis_probability_noGFD_all[, starting_age_column]
    transition_matrices[, c(11:20), "Undiagnosed CD no complications", "Undiagnosed CD osteoporosis"] <- osteoporosis_probability_noGFD_all[, starting_age_column + 1]
    transition_matrices[, c(21:30), "Undiagnosed CD no complications", "Undiagnosed CD osteoporosis"] <- osteoporosis_probability_noGFD_all[, starting_age_column + 2]
    transition_matrices[, c(31:40), "Undiagnosed CD no complications", "Undiagnosed CD osteoporosis"] <- osteoporosis_probability_noGFD_all[, starting_age_column + 3]
    transition_matrices[, c(41:50), "Undiagnosed CD no complications", "Undiagnosed CD osteoporosis"] <- osteoporosis_probability_noGFD_all[, starting_age_column + 4]
    transition_matrices[, , "Undiagnosed CD no complications", "Undiagnosed CD NHL"] <- NHL_probability_noGFD
    for (i_age in 1:50){
      transition_matrices[, i_age, "Undiagnosed CD no complications", "Death"] <- death_probability_nocomplications[i_age, 2]
    }
    
       transition_matrices[, , "Undiagnosed CD no complications", "Undiagnosed CD no complications"] <- 1-  transition_matrices[, ,"Undiagnosed CD no complications", "CD GFD no complications"] - 
         transition_matrices[, ,"Undiagnosed CD no complications", "CD no GFD no complications"] - transition_matrices[, , "Undiagnosed CD no complications", "Undiagnosed CD subfertility"] -
         transition_matrices[, , "Undiagnosed CD no complications", "Undiagnosed CD osteoporosis"] -  transition_matrices[, , "Undiagnosed CD no complications", "Undiagnosed CD NHL"] -  transition_matrices[, , "Undiagnosed CD no complications", "Death"] 
    
       transition_matrices[, ,"Undiagnosed CD subfertility", "CD GFD subfertility"] <- adherence * probability_late_diagnosis
       transition_matrices[, ,"Undiagnosed CD subfertility", "CD no GFD subfertility"] <- (1-adherence) * probability_late_diagnosis
       transition_matrices[, c(1:10),"Undiagnosed CD subfertility", "Undiagnosed CD osteoporosis"] <- osteoporosis_probability_noGFD_all[, starting_age_column]
       transition_matrices[, c(11:20), "Undiagnosed CD subfertility", "Undiagnosed CD osteoporosis"] <- osteoporosis_probability_noGFD_all[, starting_age_column + 1]
       transition_matrices[, c(21:30), "Undiagnosed CD subfertility", "Undiagnosed CD osteoporosis"] <- osteoporosis_probability_noGFD_all[, starting_age_column + 2]
       transition_matrices[, c(31:40), "Undiagnosed CD subfertility", "Undiagnosed CD osteoporosis"] <- osteoporosis_probability_noGFD_all[, starting_age_column + 3]
       transition_matrices[, c(41:50), "Undiagnosed CD subfertility", "Undiagnosed CD osteoporosis"] <- osteoporosis_probability_noGFD_all[, starting_age_column + 4]
       transition_matrices[, ,"Undiagnosed CD subfertility", "Undiagnosed CD NHL"] <- NHL_probability_noGFD
       for (i_age in 1:50){
         transition_matrices[, i_age, "Undiagnosed CD subfertility", "Death"] <- death_probability_nocomplications[i_age, 2]
       }
      transition_matrices[, , "Undiagnosed CD subfertility", "Undiagnosed CD subfertility"] <- 1 - transition_matrices[, ,"Undiagnosed CD subfertility", "CD GFD subfertility"]  - 
        transition_matrices[, ,"Undiagnosed CD subfertility", "CD no GFD subfertility"]  - 
        transition_matrices[, ,"Undiagnosed CD subfertility", "Undiagnosed CD osteoporosis"] -  transition_matrices[, ,"Undiagnosed CD subfertility", "Undiagnosed CD NHL"]  - transition_matrices[, , "Undiagnosed CD subfertility", "Death"]
   
    transition_matrices[, ,"Undiagnosed CD osteoporosis", "CD GFD osteoporosis"] <- adherence * probability_late_diagnosis
    transition_matrices[, ,"Undiagnosed CD osteoporosis", "CD no GFD osteoporosis"] <- (1-adherence) * probability_late_diagnosis
    transition_matrices[, ,"Undiagnosed CD osteoporosis", "Undiagnosed CD NHL"] <- NHL_probability_noGFD
    for (i_age in 1:50){
      transition_matrices[, i_age, "Undiagnosed CD osteoporosis", "Death"] <- death_probability_osteoporosis[i_age, 2]
    }
    transition_matrices[, ,"Undiagnosed CD osteoporosis", "Undiagnosed CD osteoporosis"] <- 1 -  transition_matrices[, ,"Undiagnosed CD osteoporosis", "CD GFD osteoporosis"] - 
      transition_matrices[, ,"Undiagnosed CD osteoporosis", "CD no GFD osteoporosis"] - transition_matrices[, , "Undiagnosed CD osteoporosis", "Death"] - transition_matrices[, ,"Undiagnosed CD osteoporosis", "Undiagnosed CD NHL"]
    
    transition_matrices[, ,"Undiagnosed CD NHL", "CD GFD NHL"] <- adherence * probability_late_diagnosis
    transition_matrices[, ,"Undiagnosed CD NHL", "CD no GFD NHL"] <- (1-adherence) * probability_late_diagnosis
      transition_matrices[, , "Undiagnosed CD NHL", "Death"] <- death_probability_NHL
      transition_matrices[, ,"Undiagnosed CD NHL", "Undiagnosed CD NHL"] <- 1 -  transition_matrices[, ,"Undiagnosed CD NHL", "CD GFD NHL"] - 
        transition_matrices[, ,"Undiagnosed CD NHL", "CD no GFD NHL"] - transition_matrices[, , "Undiagnosed CD NHL", "Death"] 

    
    transition_matrices[, ,"Death", "Death"] <- 1
    
    transition_matrices[, , , ] [is.na(transition_matrices[, , , ] )] <- 0
    
    rowSums (transition_matrices[1, 2, , ], na.rm = FALSE, dims = 1)
    # State utilities
    # Anything between 0 and 1
    
    eq5dGFD <- 0.85
    eq5dGFdse <- (0.86-0.84/3.92)
    eq5dGFDalpha <- (eq5dGFD ^ 2 * (1 - eq5dGFD)/eq5dGFdse ^ 2) - eq5dGFD
    eq5dGFDbeta <- (eq5dGFDalpha / eq5dGFD) - eq5dGFDalpha
    state_qalys[, "CD GFD no complications"] <- rbeta(n = n_samples, shape1 = eq5dGFDalpha, shape2 = eq5dGFDbeta)
    state_qalys[, "CD GFD subfertility"] <- 0.85 - 0.07
    state_qalys[,"CD GFD osteoporosis"] <- 0.7
    state_qalys[,"CD GFD NHL"] <- 0.618
    state_qalys[,"CD no GFD no complications"] <- 0.71
    state_qalys[,"CD no GFD subfertility"] <- 0.71 - 0.07
    state_qalys[, "CD no GFD osteoporosis"] <- 0.5 
    state_qalys[,"CD no GFD NHL"] <- 0.618
    state_qalys[,"Undiagnosed CD no complications"] <- 0.65
    state_qalys[,"Undiagnosed CD subfertility"] <- 0.65 - 0.07
    state_qalys[,"Undiagnosed CD osteoporosis"] <- 0.5 
    state_qalys[,"Undiagnosed CD NHL"] <- 0.618
    state_qalys[,"Death"] <- 0
    
    
    # State costs
    # Assumed normal with sd small enough to avoid negative values
    state_costs[, "CD GFD no complications"] <- 50
    state_costs[, "CD GFD subfertility"] <- 60
    state_costs[,"CD GFD osteoporosis"] <- 10
      state_costs[,"CD GFD NHL"] <- 30
    state_costs[,"CD no GFD no complications"] <- 40
    state_costs[,"CD no GFD subfertility"] <- 70
    state_costs[, "CD no GFD osteoporosis"] <-80
      state_costs[,"CD no GFD NHL"] <- 90
    state_costs[,"Undiagnosed CD no complications"] <- 100
    state_costs[,"Undiagnosed CD subfertility"] <- 110
    state_costs[,"Undiagnosed CD osteoporosis"] <- 120
      state_costs[,"Undiagnosed CD NHL"] <- 30
    state_costs[,"Death"] <- 0
    
  #How to add costs and disutilities for anemia?
  
  
  # Define the treatment costs
  # One for each PSA sample and each treatment
  # Treatment costs are actually fixed but this allows flexibility if we
  # want to include uncertainty/randomness in the cost
  treatment_costs<-array(dim=c(n_treatments,n_samples),dimnames=list(t_names,NULL))
  
  # Cost of the smoking cessation website is a one-off subscription fee of ?50
  treatment_costs["Test", ] <-50
  # Zero cost for standard of care
  treatment_costs["Test + biopsy", ] <-1000
  treatment_costs["Double test", ] <-500
  
  #############################################################################
  ## Accuracy of tests ########################################################
  #############################################################################
  
  tp <- fn <- fp <- tn <- matrix(nrow=n_samples, ncol=n_treat)
  
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
  
  # Build an array to store the cohort vector at each cycle
  # Each cohort vector has n_states elements: probability of being in each state,
  # There is one cohort vector for each treatment, for each PSA sample, for each cycle.
  cohort_vectors<-array(dim=c(n_treatments,n_samples,n_cycles,n_states),  #do I need to add another n_states?
                        dimnames=list(t_names,NULL,NULL, state_names))
  
  # This will be related to decision tree accuracy
  cohort_vectors[1, , 1,"CD GFD no complications"] <- tp[, 1] * probability_nocomplications * adherence
  cohort_vectors[2, , 1,"CD GFD no complications"] <- tp[, 2] * probability_nocomplications * adherence
  cohort_vectors[3, , 1,"CD GFD no complications"] <- tp[, 3] * probability_nocomplications * adherence
  cohort_vectors[1, , 1,"CD GFD subfertility"] <- tp[, 1] * probability_subfertility * adherence
  cohort_vectors[2, , 1,"CD GFD subfertility"] <- tp[, 2] * probability_subfertility * adherence
  cohort_vectors[3, , 1,"CD GFD subfertility"] <- tp[, 3] * probability_subfertility * adherence
  cohort_vectors[1, , 1,"CD GFD osteoporosis"] <- tp[, 1] * probability_osteoporosis * adherence
  cohort_vectors[2, , 1,"CD GFD osteoporosis"] <- tp[, 2] * probability_osteoporosis * adherence
  cohort_vectors[3, , 1,"CD GFD osteoporosis"] <- tp[, 3] * probability_osteoporosis * adherence
  cohort_vectors[1, , 1,"CD GFD NHL"] <- tp[, 1] * probability_NHL * adherence
  cohort_vectors[2, , 1,"CD GFD NHL"] <- tp[, 2] * probability_NHL * adherence
  cohort_vectors[3, , 1,"CD GFD NHL"] <- tp[, 3] * probability_NHL * adherence
  
  cohort_vectors[1, , 1,"CD no GFD no complications"] <- tp[, 1] * probability_nocomplications * (1-adherence)
  cohort_vectors[2, , 1,"CD no GFD no complications"] <- tp[, 2] * probability_nocomplications * (1-adherence)
  cohort_vectors[3, , 1,"CD no GFD no complications"] <- tp[, 3] * probability_nocomplications * (1-adherence)
  cohort_vectors[1, , 1,"CD no GFD subfertility"] <- tp[, 1] * probability_subfertility * (1-adherence)
  cohort_vectors[2, , 1,"CD no GFD subfertility"] <- tp[, 2] * probability_subfertility * (1-adherence)
  cohort_vectors[3, , 1,"CD no GFD subfertility"] <- tp[, 3] * probability_subfertility * (1-adherence)
  cohort_vectors[1, , 1,"CD no GFD osteoporosis"] <- tp[, 1] * probability_osteoporosis * (1-adherence)
  cohort_vectors[2, , 1,"CD no GFD osteoporosis"] <- tp[, 2] * probability_osteoporosis * (1-adherence)
  cohort_vectors[3, , 1,"CD no GFD osteoporosis"] <- tp[, 3] * probability_osteoporosis * (1-adherence)
  cohort_vectors[1, , 1,"CD no GFD NHL"] <- tp[, 1] * probability_NHL * (1-adherence)
  cohort_vectors[2, , 1,"CD no GFD NHL"] <- tp[, 2] * probability_NHL * (1-adherence)
  cohort_vectors[3, , 1,"CD no GFD NHL"] <- tp[, 3] * probability_NHL * (1-adherence)
  
  cohort_vectors[1, , 1,"Undiagnosed CD no complications"] <- fn[, 1] * probability_nocomplications 
  cohort_vectors[2, , 1,"Undiagnosed CD no complications"] <- fn[, 2] * probability_nocomplications 
  cohort_vectors[3, , 1,"Undiagnosed CD no complications"] <- fn[, 3] * probability_nocomplications 
  cohort_vectors[1, , 1,"Undiagnosed CD subfertility"] <- fn[, 1] * probability_subfertility 
  cohort_vectors[2, , 1,"Undiagnosed CD subfertility"] <- fn[, 2] * probability_subfertility 
  cohort_vectors[3, , 1,"Undiagnosed CD subfertility"] <- fn[, 3] * probability_subfertility 
  cohort_vectors[1, , 1,"Undiagnosed CD osteoporosis"] <- fn[, 1] * probability_osteoporosis 
  cohort_vectors[2, , 1,"Undiagnosed CD osteoporosis"] <- fn[, 2] * probability_osteoporosis 
  cohort_vectors[3, , 1,"Undiagnosed CD osteoporosis"] <- fn[, 3] * probability_osteoporosis 
  cohort_vectors[1, , 1,"Undiagnosed CD NHL"] <- fn[, 1] * probability_NHL 
  cohort_vectors[2, , 1,"Undiagnosed CD NHL"] <- fn[, 2] * probability_NHL 
  cohort_vectors[3, , 1,"Undiagnosed CD NHL"] <- fn[, 3] * probability_NHL 
  
  cohort_vectors[, , , ] [is.na(cohort_vectors[, , , ] )] <- 0
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
  
  #i_treatment <- 1
  #i_sample <- 1
  #i_cycle <- 2
  
  disc_vec <- (1/1.035) ^ rep(c(0:(n_cycles/2-1)), each = 2)
  
  # The remainder of the cohort_vectors will be filled in by Markov updating below
  
  # Main model code
  # Loop over the treatment options
  

    for (i_treatment in 1:n_treatments)
    {
    # Loop over the PSA samples
    for(i_sample in 1:n_samples)
    {
      
      transition_matrices_sample <- transition_matrices[i_sample, , , ]
      # Loop over the cycles
      # Cycle 1 is already defined so only need to update cycles 2:n_cycles
      for(i_cycle in 2:n_cycles)
      {
        # Markov update
        # Multiply previous cycle's cohort vector by transition matrix
        # i.e. pi_j = pi_(j-1)*P
        cohort_vectors[i_treatment, i_sample, i_cycle, ] <-
          cohort_vectors[i_treatment, i_sample, i_cycle-1, ] %*%
          transition_matrices_sample[i_cycle, , ]
      }
      
      cohort_vectors_tr_sample <- cohort_vectors[i_treatment, i_sample, , ]
      
      # Now use the cohort vectors to calculate the 
      # total costs for each cycle
      cycle_costs[i_treatment, i_sample, ] <-
        cohort_vectors_tr_sample %*% state_costs[i_sample, ]
      # And total QALYs for each cycle
      cycle_qalys[i_treatment, i_sample, ] <-
        cohort_vectors_tr_sample %*% state_qalys[i_sample, ]
      
      # Combine the cycle_costs and treatment_costs to get total costs
      # Apply the discount factor 
      # (1 in first year, 1.035 in second, 1.035^2 in third, and so on)
      # Each year acounts for two cycles so need to repeat the discount values
      total_costs[i_treatment, i_sample] <- treatment_costs[i_treatment, i_sample] +
        cycle_costs[i_treatment, i_sample, ] %*%
        disc_vec
      
      # Combine the cycle_qalys to get total qalys
      # Apply the discount factor 
      # (1 in first year, 1.035 in second, 1.035^2 in third, and so on)
      # Each year acounts for two cycles so need to repeat the discount values
      total_qalys[i_treatment, i_sample] <- cycle_qalys[i_treatment, i_sample, ] %*%
        disc_vec
    }
  
    }
  #############################################################################
  ## Analysis of results ######################################################
  #############################################################################
  output <- list()
  # Average costs
  # These are ?50 on the website and 0 on standard of care as there are no
  # costs other than the website subscription cost
  output$average_costs <- rowMeans(total_costs)
  # Average effects (in QALY units)
  # These are slightly higher on the website as higher probability of 
  # quitting smoking
  output$average_effects <- rowMeans(total_qalys)
  
  # Incremental costs and effects relative to standard of care
  # No uncertainty in the costs as the website cost is fixed at ?50
  output$incremental_costs <- total_costs["Test + biopsy", ] - total_costs["Test", ]
  # In some samples the website leads to higher QALYs but in others it is negative
  # There is uncertainty as to whether the website is an improvement over SoC
  output$incremental_effects <- total_qalys["Test + biopsy", ] - total_qalys["Test", ]
  
  # The ICER comparing Standard of care with website to standard of care
  # This is much lower than the ?20,000 willingness-to-pay threshold indicating
  # good value for money
  output$ICER <- mean(output$incremental_costs)/mean(output$incremental_effects)
  
  # Incremental net benefit at the ?20,000 willingness-to-pay
  # Sometimes positive (website more cost-effective) and sometimes negative (SoC more cost-effective)
  # Need to look at averages and consider probabilities of cost-effectiveness
  output$incremental_net_benefit <- 20000*output$incremental_effects - output$incremental_costs
  
  # Average incremental net benefit
  # This is positive indicating cost-effectiveness at the ?20,000 threshold
  output$average_inb <- mean(output$incremental_net_benefit)
  
  # Probability cost-effective
  # This is the proportion of samples for which the incremental net benefit is positive
  # It is clost to 72%, representing good degree of certainty
  # in recommendation to adopt the smoking cessation website
  output$probability_cost_effective <- sum(output$incremental_net_benefit>0)/n_samples
  
  # Now use the BCEA package to analyse the results___
  output


