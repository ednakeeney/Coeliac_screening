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
  # HT: I think you need to divide by sqrt(sample size) to get SE, so fewer would be negative.
  # HT: However, I would switch to a lognormal distribution so they're always positive
  # HT: Why is the calculation using exponential commented out?
  duration_of_symptoms <- rnorm(n = n_samples, mean = 10.93, sd = 14.88)     #calculated from Violato, what about SD? Should some of these be negative?
  rate_of_symptoms <- 1 / duration_of_symptoms
  #probability_late_diagnosis <- 1 - exp(-rate_of_symptoms)
  probability_late_diagnosis <- 0.06

    prevalence <- read.csv("CPRD prevalence.csv")
  
  #Initial cohort at diagnosis - depends on age at diagnosis
    # HT: When accesseing named matrices, try to use the column names to avoid errors and improve readability (like you do for osteoporosis_probability)
 probability_subfertility <- rbeta(n=n_samples, shape1 = prevalence[prevalence$Age.categories == starting_age, 5], shape2 = prevalence[prevalence$Age.categories == starting_age, 2] - prevalence[prevalence$Age.categories == starting_age, 5])
 probability_osteoporosis <- rbeta(n=n_samples, shape1 = prevalence[prevalence$Age.categories == starting_age, 4], shape2 = prevalence[prevalence$Age.categories == starting_age, 2] - prevalence[prevalence$Age.categories == starting_age, 4])
 probability_NHL <- rbeta(n=n_samples, shape1 = prevalence[prevalence$Age.categories == starting_age, 3], shape2 = prevalence[prevalence$Age.categories == starting_age, 2] - prevalence[prevalence$Age.categories == starting_age, 3])
 probability_nocomplications <- rbeta(n=n_samples, shape1 = prevalence[prevalence$Age.categories == starting_age, 7], shape2 = prevalence[prevalence$Age.categories == starting_age, 2] - prevalence[prevalence$Age.categories == starting_age, 7])
  
  # Osteoporosis probabilities On GFD
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

 # Subfertility probabilities on GFD
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
  
  # NHL probabilities on GFD
  NHL_probability <- read.csv("NHL_rate_nice.csv")
  NHL_probability_GFD_18orless	<- rbeta(n=n_samples, shape1 = as.numeric(NHL_probability$NHL_GFD_alpha[1]), shape2 = as.numeric(NHL_probability$NHL_GFD_beta[1]))
  NHL_probability_GFD_18plus	<- rbeta(n=n_samples, shape1 = as.numeric(NHL_probability$NHL_GFD_alpha[2]), shape2 = as.numeric(NHL_probability$NHL_GFD_beta[2]))
  NHL_probability_GFD <- NHL_probability_GFD_18plus
  
  # Corresponding probabilities not on GFD
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
  # HT: I had an error trying to read the above file
  # Error in read.table(file = file, header = header, sep = sep, quote = quote,  : 
  # more columns than column names
  
  # Just to check the rest of the code I inserted placeholder values but please delete these two lines
  lifetables <- as.data.frame(matrix(c(1:100, rep(seq(0.0001, 0.38, (0.38-0.0001)/99), 2)), ncol = 3))
  colnames(lifetables) <- c("Age", "Males", "Females")
  percentage_male <- 0.5
  lifetables$Overall <- (percentage_male * lifetables$Males) + ((1-percentage_male) * lifetables$Females)
 # lifetables30plus <- subset(lifetables, lifetables$Age > 29)
  death_probability_nocomplications	<- data.frame(lifetables$Age, lifetables$Overall)
 
  # HT: This gave the error
  # Error in read.table(file = file, header = header, sep = sep, quote = quote,  : 
  # more columns than column names
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
  
   
   # HT: Everything up to this point should go in the generate_model_parameters() function which exports an input_parameters matrix
   # The matrix will be large but will make it easy later on to use SAVI or MLMC for EVPPI
   # HT: Also code below about costs and QALYs
   # HT: The code that follows should go in a a generate_transition_matrices() function which
   # takes input_parameters as input and is called within the (main model running) generate_net_benefit() function
   
   # There is one transition matrix for each cycle and each PSA sample
  # Store them in an array with (before filling in below) NA entries
  transition_matrices <- array(dim = c(n_samples, n_cycles, n_states, n_states),
                             dimnames = list(NULL, NULL, state_names, state_names))

  
  # HT: Are we going to use QALY norms so that the QALYs decrease with age? You'd multiple all
  # state QALYs by their norms to reduce with age.
  # HT: If yes, we'll need state_qalys to be cycle specific and include it in the protocol.
  # Now define the QALYS associated with the states per cycle
  # There is one for each PSA sample and each state
  # Store in an NA array and then fill in below
  state_qalys <- array(dim=c(n_samples, n_states), dimnames = list(NULL, state_names))
  
  # And finally define the state costs
  # There is one for each PSA sample and each state
  # Store in an NA array and then fill in below
  state_costs <- array(dim=c(n_samples, n_cycles, n_states), dimnames = list(NULL, NULL , state_names))
  
  
  # Define transition matrices, state utilities and costs 

  # HT: The lines below should instead go in this for loop. I've added it elsewhere and they can probably
  # all be combined into a single loop
  # HT: Use this idea whereever possible to reduce script length
  #CD GFD 
  for(i_age_category in c(0:4)) {
    transition_matrices[, (c(1:10) + i_age_category * 10), "CD GFD no complications", "CD GFD subfertility"] <- subfertility_probability_GFD_all[, starting_age_column + i_age_category]  
    transition_matrices[, (c(1:10) + i_age_category * 10), "CD GFD no complications", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_all[, starting_age_column + i_age_category]
  }
  
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

    # HT: Error in the below. I think i_age should be added to starting_age when accessing the lifetables.
    # HT: Same error elsewhere below
    # HT: Also, as for the loop over age categories, please put all of these into a single for loop over i_age.
    for (i_age in 1:50){
    transition_matrices[, i_age, "CD GFD no complications", "Death"] <- death_probability_nocomplications[i_age, 2]
    }
   
    # HT: This is setting the probability of staying in the same state to be 1 minus the probability
    # of going to any other state. It would be better practice if this were done at the end for every state.
    # At the end (when every probability has been defined) you also wouldn't need to write out the probabilities one by one. 
    transition_matrices[, , "CD GFD no complications", "CD GFD no complications"] <- 1 - transition_matrices[, , "CD GFD no complications", "CD GFD subfertility"] - 
      transition_matrices[, , "CD GFD no complications", "CD GFD osteoporosis"] - transition_matrices[, , "CD GFD no complications", "CD GFD NHL"] -  transition_matrices[ , , "CD GFD no complications", "Death"]


    # HT: Again can replace the lines below with this for loop
    for(i_age_category in c(0:4)) {
      transition_matrices[, (c(1:10) + i_age_category * 10), "CD GFD subfertility", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_all[, starting_age_column + i_age_category]
    }
    
    transition_matrices[, c(1:10),"CD GFD subfertility", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_all[, starting_age_column]
    transition_matrices[, c(11:20), "CD GFD subfertility", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_all[, starting_age_column + 1]
    transition_matrices[, c(21:30), "CD GFD subfertility", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_all[, starting_age_column + 2]
    transition_matrices[, c(31:40), "CD GFD subfertility", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_all[, starting_age_column + 3]
    transition_matrices[, c(41:50), "CD GFD subfertility", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_all[, starting_age_column + 4]
    
    transition_matrices[, ,"CD GFD subfertility", "CD GFD NHL"] <- NHL_probability_GFD
    # HT: Same error with indexing as above, and need to combine all these for loops
    for (i_age in 1:50){
      transition_matrices[, i_age, "CD GFD subfertility", "Death"] <- death_probability_nocomplications[i_age, 2]
    }
    
    # HT: Again put all these at the end and rewrite as sum of other transition probabilies
      transition_matrices[, , "CD GFD subfertility", "CD GFD subfertility"] <- 1 -  transition_matrices[, ,"CD GFD subfertility", "CD GFD osteoporosis"] - 
        transition_matrices[, ,"CD GFD subfertility", "CD GFD NHL"] - transition_matrices[, , "CD GFD subfertility", "Death"]
    
    
    transition_matrices[, ,"CD GFD osteoporosis", "CD GFD NHL"] <- NHL_probability_GFD
    for (i_age in 1:50){
    transition_matrices[, i_age, "CD GFD osteoporosis", "Death"] <- death_probability_osteoporosis[i_age, 2]
    }
    
    # HT: Move to end and rewrite
      transition_matrices[, , "CD GFD osteoporosis", "CD GFD osteoporosis"] <- 1 - 
        transition_matrices[, ,"CD GFD osteoporosis", "CD GFD NHL"] - transition_matrices[, , "CD GFD osteoporosis", "Death"]
   
    transition_matrices[, , "CD GFD NHL", "Death"] <- death_probability_NHL
    transition_matrices[, ,"CD GFD NHL", "CD GFD NHL"] <- 1 - death_probability_NHL
   

#CD no GFD    
    # HT: Again use this for loop instead of code below (and combine with other loops over age categories)
    for(i_age_category in c(0:4)) {
      transition_matrices[, (c(1:10) + i_age_category * 10), "CD no GFD no complications", "CD no GFD subfertility"] <- subfertility_probability_noGFD_all[, starting_age_column + i_age_category]
      transition_matrices[, (c(1:10) + i_age_category * 10), "CD no GFD no complications", "CD no GFD osteoporosis"] <- osteoporosis_probability_noGFD_all[, starting_age_column + i_age_category]
    }
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
    
    # HT: I had a an error that NHL_probability_noGFD was not defined
    transition_matrices[, , "CD no GFD no complications", "CD no GFD NHL"] <- NHL_probability_noGFD
    for (i_age in 1:50){
      transition_matrices[, i_age, "CD no GFD no complications", "Death"] <- death_probability_nocomplications[i_age, 2]
    }
    
    # HT: Again at end after all transitions defined
    transition_matrices[, , "CD no GFD no complications", "CD no GFD no complications"] <- 1 - transition_matrices[, , "CD no GFD no complications", "CD no GFD subfertility"] - 
      transition_matrices[, , "CD no GFD no complications", "CD no GFD osteoporosis"] - transition_matrices[, , "CD no GFD no complications", "CD no GFD NHL"] -  transition_matrices[ , , "CD no GFD no complications", "Death"]
    
    
    # HT: Again use this for loop instead of code below (and combine with other loops over age categories)
    for(i_age_category in c(0:4)) {
      transition_matrices[, (c(1:10) + i_age_category * 10), "CD no GFD subfertility", "CD no GFD osteoporosis"] <- osteoporosis_probability_noGFD_all[, starting_age_column + i_age_category]
    }
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
   
    # HT: Again use a loop and combine with above
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
    
    # HT: Correct index and combine with otehr age loops
    for (i_age in 1:50){
      transition_matrices[, i_age, "Undiagnosed CD no complications", "Death"] <- death_probability_nocomplications[i_age, 2]
    }
    
    # HT: Again move to end and sum over all probabilities to avoid errors
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
    
    # HT: This is a little dangerous as you may have missed something by accident.
    # Please look at one matrix (e.g. transition_matrices[1, 1, , ]) and check that each NA really should be NA. We could perhaps do this together on our next call.
    transition_matrices[, , , ] [is.na(transition_matrices[, , , ] )] <- 0
    
    # HT: Obviously remove this check, or include it in the generate_transition_matrices() function
    # and have it output a warning() if the sums are not 1.
    rowSums (transition_matrices[1, 2, , ], na.rm = FALSE, dims = 1)
    
    # HT: What follows will go in a generate_state_qalys() function
    # which takes input_parameters as input and is called from within generate_net_benefit()
    # State utilities
    # Anything between 0 and 1
    
    eq5dGFD <- 0.85
    eq5dGFdse <- ((0.86-0.84)/3.92)
    eq5dGFDalpha <- (eq5dGFD ^ 2 * (1 - eq5dGFD)/eq5dGFdse ^ 2) - eq5dGFD
    eq5dGFDbeta <- (eq5dGFDalpha / eq5dGFD) - eq5dGFDalpha
    state_qalys[, "CD GFD no complications"] <- rbeta(n = n_samples, shape1 = eq5dGFDalpha, shape2 = eq5dGFDbeta)
    
    disutility_subfertility <-  0.158 
    disutility_subfertility_se <- (0.173 - 0.143)/3.92
    disutility_subfertility_alpha <- (disutility_subfertility ^ 2 * (1 - disutility_subfertility)/disutility_subfertility_se ^ 2) - disutility_subfertility
    disutility_subfertility_beta <- (disutility_subfertility_alpha/disutility_subfertility) - disutility_subfertility_alpha
    disutility_subfertility <- rbeta(n = n_samples, shape1 = disutility_subfertility_alpha, shape2 = disutility_subfertility_beta)
    state_qalys[, "CD GFD subfertility"] <- state_qalys[, "CD GFD no complications"] - disutility_subfertility
   
    probability_hipfracture <- 0.00196
    probability_vertebralfracture <- 0.00071
    probability_wristfracture <- 0.00125
    disutility_hipfracture <- 0.817 - 0.59
    disutility_hipfractureSE <- (((0.817 - 0.65) - (0.817 - 0.54))/3.92)
    disutility_hipfracture_alpha <- (disutility_hipfracture ^ 2 * (1 - disutility_hipfracture)/disutility_hipfractureSE ^ 2) - disutility_hipfracture
    disutility_hipfracture_beta <- (disutility_hipfracture_alpha/disutility_hipfracture) - disutility_hipfracture_alpha
    disutility_hipfracture <- rbeta(n = n_samples, shape1 = disutility_hipfracture_alpha, shape2 = disutility_hipfracture_beta)
    disutility_wristfracture <- 0.817 - 0.55
    disutility_wristfractureSE <- (((0.817 - 0.60) - (0.817 - 0.50))/3.92)
    disutility_wristfracture_alpha <- (disutility_wristfracture ^ 2 * (1 - disutility_wristfracture)/disutility_wristfractureSE ^ 2) - disutility_wristfracture
    disutility_wristfracture_beta <- (disutility_wristfracture_alpha/disutility_wristfracture) - disutility_wristfracture_alpha
    disutility_wristfracture <- rbeta(n = n_samples, shape1 = disutility_wristfracture_alpha, shape2 = disutility_wristfracture_beta)
    disutility_vertebralfracture <- 0.817 - 0.78
    disutility_vertebralfractureSE <- (((0.817 - 0.84) - (0.817 - 0.72))/3.92)
    disutility_vertebralfracture_alpha <- (disutility_vertebralfracture ^ 2 * (1 - disutility_vertebralfracture)/disutility_vertebralfractureSE ^ 2) - disutility_vertebralfracture
    disutility_vertebralfracture_beta <- (disutility_vertebralfracture_alpha/disutility_vertebralfracture) - disutility_vertebralfracture_alpha
    disutility_vertebralfracture <- rbeta(n = n_samples, shape1 = disutility_vertebralfracture_alpha, shape2 = disutility_vertebralfracture_beta)
    disutility_osteoporosis <- (probability_hipfracture * disutility_hipfracture) + (probability_wristfracture * disutility_wristfracture) + (probability_vertebralfracture * disutility_vertebralfracture)
    state_qalys[,"CD GFD osteoporosis"] <- state_qalys[, "CD GFD no complications"] - disutility_osteoporosis
    
    disutility_NHL <- runif(n = n_samples, min = 0.036, max = 0.136)
    state_qalys[,"CD GFD NHL"] <- state_qalys[, "CD GFD no complications"] - disutility_NHL
    
    disutility_noGFD <-  0.14 
    disutility_noGFD_se <- (0.31 - (-0.03))/3.92
    disutility_noGFD_alpha <- (disutility_noGFD ^ 2 * (1 - disutility_noGFD)/disutility_noGFD_se ^ 2) - disutility_noGFD
    disutility_noGFD_beta <- (disutility_noGFD_alpha/disutility_noGFD) - disutility_noGFD_alpha
    disutility_noGFD <- rbeta(n = n_samples, shape1 = disutility_noGFD_alpha, shape2 = disutility_noGFD_beta)
    state_qalys[,"CD no GFD no complications"] <- state_qalys[, "CD GFD no complications"] - disutility_noGFD
    
    state_qalys[,"CD no GFD subfertility"] <-  state_qalys[,"CD no GFD no complications"] - disutility_subfertility
    
    state_qalys[, "CD no GFD osteoporosis"] <- state_qalys[,"CD no GFD no complications"] - disutility_osteoporosis
    
    state_qalys[,"CD no GFD NHL"] <- state_qalys[,"CD no GFD no complications"] - disutility_NHL
    
    disutility_undiagnosedCD <-  0.65 
    disutility_undiagnosedCD_se <- (0.67 - 0.63)/3.92
    disutility_undiagnosedCD_alpha <- (disutility_undiagnosedCD ^ 2 * (1 - disutility_undiagnosedCD)/disutility_undiagnosedCD_se ^ 2) - disutility_undiagnosedCD
    disutility_undiagnosedCD_beta <- (disutility_undiagnosedCD_alpha/disutility_undiagnosedCD) - disutility_undiagnosedCD_alpha
    disutility_undiagnosedCD <- rbeta(n = n_samples, shape1 = disutility_undiagnosedCD_alpha, shape2 = disutility_undiagnosedCD_beta)
    state_qalys[,"Undiagnosed CD no complications"] <- disutility_undiagnosedCD
    
    state_qalys[,"Undiagnosed CD subfertility"] <-  state_qalys[,"Undiagnosed CD no complications"] - disutility_subfertility
    
    state_qalys[,"Undiagnosed CD osteoporosis"] <- state_qalys[,"Undiagnosed CD no complications"] - disutility_osteoporosis
    
    state_qalys[,"Undiagnosed CD NHL"] <- state_qalys[,"Undiagnosed CD no complications"] - disutility_NHL
    
    state_qalys[,"Death"] <- 0
    
    
    # State costs
    # Assumed normal with sd small enough to avoid negative values
    cost_CDGFD <- 650
    cost_CDGFD_se <- 4.68
    cost_CDGFD_alpha <- (cost_CDGFD/cost_CDGFD_se)^2
    cost_CDGFD_beta <- (cost_CDGFD_se^2)/cost_CDGFD
    state_costs[, , "CD GFD no complications"] <- state_costs[, , "CD no GFD no complications"] <- rgamma(n = n_samples, shape = cost_CDGFD_alpha, scale = cost_CDGFD_beta)
   
    cost_subfertility <- 8079.75 
    cost_subfertility_se <- (8742.96 - 7432.11)/3.92
    cost_subfertility_alpha <- (cost_subfertility/cost_subfertility_se)^2
    cost_subfertility_beta <- (cost_subfertility_se^2)/cost_subfertility
    state_costs[, 1, "CD GFD subfertility"] <-  state_costs[, , "CD no GFD subfertility"] <- (rgamma(n = n_samples, shape = cost_subfertility_alpha, scale = cost_subfertility_beta)) +  rgamma(n = n_samples, shape = cost_CDGFD_alpha, scale = cost_CDGFD_beta)
    state_costs[, 2:n_cycles, "CD GFD subfertility"] <- state_costs[, 2:n_cycles, "CD GFD no complications"]
    
    cost_hipfracture <- 19073
    cost_hipfractureSE <- ((16515 * 1.17) - (16097 * 1.17)) / 3.92
    cost_hipfracture_alpha <- (cost_hipfracture / cost_hipfractureSE) ^ 2
    cost_hipfracture_beta <- (cost_hipfractureSE ^ 2) / cost_hipfracture
    cost_hipfracture <- rgamma(n = n_samples, shape = cost_hipfracture_alpha, scale = cost_hipfracture_beta)
    cost_wristfracture <- 842.40
    cost_wristfractureSE <- cost_wristfracture/10
    cost_wristfracture_alpha <- (cost_wristfracture / cost_wristfractureSE) ^ 2
    cost_wristfracture_beta <- (cost_wristfractureSE ^ 2) / cost_wristfracture
    cost_wristfracture <- rgamma(n = n_samples, shape = cost_wristfracture_alpha, scale = cost_wristfracture_beta)
    cost_vertebralfracture <- 862.20
    cost_vertebralfractureSE <- cost_vertebralfracture/10
    cost_vertebralfracture_alpha <- (cost_vertebralfracture / cost_vertebralfractureSE) ^ 2
    cost_vertebralfracture_beta <- (cost_vertebralfractureSE ^ 2) / cost_vertebralfracture
    cost_vertebralfracture <- rgamma(n = n_samples, shape = cost_vertebralfracture_alpha, scale = cost_vertebralfracture_beta)
    cost_osteoporosis <- (probability_hipfracture * cost_hipfracture) + (probability_wristfracture * cost_wristfracture) + (probability_vertebralfracture * cost_vertebralfracture)
    state_costs[, , "CD GFD osteoporosis"] <- state_costs[, , "CD no GFD osteoporosis"] <- cost_osteoporosis + rgamma(n = n_samples, shape = cost_CDGFD_alpha, scale = cost_CDGFD_beta)
    
    cost_NHL <- 18396
    cost_NHL_sd <- sqrt(271) * ((18415 - 18377)/3.92)
    cost_NHL_location <- log(cost_NHL^2 / sqrt(cost_NHL_sd^2 + cost_NHL^2))
    cost_NHL_shape <- sqrt(log(1 + (cost_NHL_sd^2 / cost_NHL^2)))
    state_costs[, 1, "CD GFD NHL"] <-  state_costs[, 1, "CD no GFD NHL"] <- rlnorm(n = n_samples, cost_NHL_location,  cost_NHL_shape) + rgamma(n = n_samples, shape = cost_CDGFD_alpha, scale = cost_CDGFD_beta)
    state_costs[, 2:n_cycles, "CD GFD NHL"] <-  state_costs[, 2:n_cycles, "CD no GFD NHL"] <- rlnorm(n = n_samples, cost_NHL_location,  cost_NHL_shape) + rgamma(n = n_samples, shape = cost_CDGFD_alpha, scale = cost_CDGFD_beta)
    
    
    cost_undiagnosedCD <- 340
    cost_undiagnosedCD_se <- 2.96
    cost_undiagnosedCD_alpha <- (cost_undiagnosedCD/cost_undiagnosedCD_se)^2
    cost_undiagnosedCD_beta <- (cost_undiagnosedCD_se^2)/cost_undiagnosedCD
    state_costs[, , "Undiagnosed CD no complications"] <- rgamma(n = n_samples, shape = cost_undiagnosedCD_alpha, scale = cost_undiagnosedCD_beta)
    
    state_costs[, , "Undiagnosed CD subfertility"] <- state_costs[, , "Undiagnosed CD no complications"] +  state_costs[, 1, "CD GFD subfertility"]
    state_costs[, , "Undiagnosed CD osteoporosis"] <- state_costs[, , "Undiagnosed CD no complications"] + state_costs[, , "CD GFD osteoporosis"]
    state_costs[, , "Undiagnosed CD NHL"] <- state_costs[, , "Undiagnosed CD no complications"] + state_costs[, 1, "CD GFD NHL"]
    state_costs[, , "Death"] <- 0
    
  #How to add costs and disutilities for anemia and biopsy?
  
  
  # Define the treatment costs
  # One for each PSA sample and each treatment
  # Treatment costs are actually fixed but this allows flexibility if we
  # want to include uncertainty/randomness in the cost
  treatment_costs<-array(dim=c(n_treatments,n_samples),dimnames=list(t_names,NULL))
  
  treatment_costs["Test", ] <- 10
  probability_biopsy <- runif(n = n_samples, min = 0.6, max = 0.8)
  treatment_costs["Test + biopsy", ] <- 10 + (probability_biopsy * 530)
  treatment_costs["Double test", ] <- 20
  
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
      
      for (i_cycle in 1:n_cycles)
      {
      # Now use the cohort vectors to calculate the 
      # total costs for each cycle
      cycle_costs[i_treatment, i_sample, ] <-
        cohort_vectors_tr_sample %*% state_costs[i_sample, i_cycle, ]
      }
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
  
 
  output$incremental_costs_testbiopsy_test <- total_costs["Test + biopsy", ] - total_costs["Test", ]
  output$incremental_effects_testbiopsy_test <- total_qalys["Test + biopsy", ] - total_qalys["Test", ]
  
  output$incremental_costs_doubletest_test <- total_costs["Double test", ] - total_costs["Test", ]
  output$incremental_effects_doubletest_test <- total_qalys["Double test", ] - total_qalys["Test", ]
  

  output$ICER_testbiopsy_test <- mean(output$incremental_costs_testbiopsy_test)/mean(output$incremental_effects_testbiopsy_test)
  output$ICER_doubletest_test <- mean(output$incremental_costs_doubletest_test)/mean(output$incremental_effects_doubletest_test)
  
  # Incremental net benefit at the ?20,000 willingness-to-pay
 
  output$incremental_net_benefit_testbiopsy_test <- 20000*output$incremental_effects_testbiopsy_test - output$incremental_costs_testbiopsy_test
  output$incremental_net_benefit_doubletest_test <- 20000*output$incremental_effects_doubletest_test - output$incremental_costs_doubletest_test
 
   # Average incremental net benefit
  output$average_inb_testbiopsy_test <- mean(output$incremental_net_benefit_testbiopsy_test)
  output$average_inb_doubletest_test <- mean(output$incremental_net_benefit_doubletest_test)

  # Probability cost-effective
  # This is the proportion of samples for which the incremental net benefit is positive
 
  output$probability_cost_effective_testbiopsy_test <- sum(output$incremental_net_benefit_testbiopsy_test > 0)/n_samples
  output$probability_cost_effective_doubletest_test <- sum(output$incremental_net_benefit_doubletest_test > 0)/n_samples
  
  # Now use the BCEA package to analyse the results___
  output


