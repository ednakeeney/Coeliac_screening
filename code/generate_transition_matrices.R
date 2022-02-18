generate_transition_matrices <- function(input_parameters, population = NULL) {
 
  n_samples <- dim(input_parameters)[1]
  if(population == "children") percentage_male <- 0.5
  if(population == "men") percentage_male <- 1.0
  if(population == "women") percentage_male <- 0.0
  
   starting_age <- ifelse(population == "men" | population == "women", 18, 10) #based on mean age in under and over 18s in CPRD cost data
  n_cycles <- 90 - starting_age
  
  transition_matrices <- array(dim = c(n_samples, n_cycles, n_states, n_states),
                               dimnames = list(NULL, NULL, state_names, state_names))
  
  # There is one transition matrix for each cycle and each PSA sample
  # Store them in an array with (before filling in below) NA entries
  # Define transition matrices
  
  starting_age_columnandrow <- read.csv("data/starting_age_column.csv")
  starting_age_columnandrow$starting_age_row <- c(1:11)
  #starting_age_column <- starting_age_columnandrow[starting_age_columnandrow$Starting.age == starting_age, 2]
  #Initial cohort at diagnosis - depends on age at diagnosis
  # e.g. 10 year old starts with 10-20 prevalence and 18 year old starts with 10-20 prevalence. 
  starting_age_column <- which(starting_age_columnandrow$Starting.age  > starting_age)[1] - 1
  
  
  
  # Osteoporosis probabilities 
  # Generate probabilities based on these rates and ratios
  osteoporosis_probability_GFD <- osteoporosis_probability_noGFD <- matrix(NA, nrow = n_samples, ncol = 10)
  
  for(i_age_category in 1:10) {
    osteoporosis_probability_GFD[, i_age_category] <- 1 - exp(-exp(input_parameters[, paste0("osteoporosis_lograte_", 10*(i_age_category - 1))] + input_parameters$log_or_osteoporosis_GFD))
    osteoporosis_probability_noGFD[, i_age_category] <- 1 - exp(-exp(input_parameters[, paste0("osteoporosis_lograte_", 10 * (i_age_category - 1))] + input_parameters$log_or_osteoporosis_noGFD))
  }
  osteoporosis_probability_GFD_all <- as.data.frame(osteoporosis_probability_GFD)
  colnames(osteoporosis_probability_GFD_all) <- c("osteoporosis_probability_GFD_0","osteoporosis_probability_GFD_10","osteoporosis_probability_GFD_20","osteoporosis_probability_GFD_30",
                                                  "osteoporosis_probability_GFD_40","osteoporosis_probability_GFD_50","osteoporosis_probability_GFD_60",
                                                  "osteoporosis_probability_GFD_70","osteoporosis_probability_GFD_80","osteoporosis_probability_GFD_90")
  osteoporosis_probability_noGFD_all <- as.data.frame(osteoporosis_probability_noGFD)
  colnames(osteoporosis_probability_noGFD_all) <- c("osteoporosis_probability_noGFD_0","osteoporosis_probability_noGFD_10","osteoporosis_probability_noGFD_20","osteoporosis_probability_noGFD_30",
                                                    "osteoporosis_probability_noGFD_40","osteoporosis_probability_noGFD_50","osteoporosis_probability_noGFD_60",
                                                    "osteoporosis_probability_noGFD_70","osteoporosis_probability_noGFD_80","osteoporosis_probability_noGFD_90")
  
  
  # Calculate NHL probabilities
  
  NHL_probability_GFD_18orless	<- 1 - exp(-exp(input_parameters$NHL_lograte_18orless + input_parameters$log_rr_NHL_GFD))
  NHL_probability_GFD_18plus	<-  1 - exp(-exp(input_parameters$NHL_lograte_18plus + input_parameters$log_rr_NHL_GFD))
  NHL_probability_noGFD_18orless	<- 1 - exp(-exp(input_parameters$NHL_lograte_18orless + input_parameters$log_rr_NHL_noGFD))
  NHL_probability_noGFD_18plus	<-  1 - exp(-exp(input_parameters$NHL_lograte_18plus + input_parameters$log_rr_NHL_noGFD))
  
  NHL_probability_GFD <- rep(0, times = n_samples)
  NHL_probability_noGFD <- rep(0, times = n_samples)
  for (i_sample in 1:n_samples) {
    NHL_probability_GFD[i_sample] <- ifelse(population == "men" | population == "women", NHL_probability_GFD_18plus[i_sample], NHL_probability_GFD_18orless[i_sample])
    NHL_probability_noGFD[i_sample] <- ifelse(population == "men" | population == "women", NHL_probability_noGFD_18plus[i_sample], NHL_probability_noGFD_18orless[i_sample])
  }
  
  lifetables <- read.csv("data/lifetables.csv")
  # Lifetables are rates so must be converted to probabilities
  lifetables$Males_lograte <- log(lifetables$Males)
  lifetables$Females_lograte <- log(lifetables$Females)
  lifetables$Males <- 1 - exp(-lifetables$Males)
  lifetables$Females <- 1 - exp(-lifetables$Females)
  
  lifetables$Overall <- (percentage_male * lifetables$Males) + ((1-percentage_male) * lifetables$Females)
  death_probability_nocomplications	<- data.frame(lifetables$Age, lifetables$Overall)
  
  
  # Combine on log scale
  # Rows are for samples, columns are for ages
  #1-exp(-exp(1.25 + lifetables$Males_lograte))
  death_probability_osteoporosis <- list()
  death_probability_osteoporosis$Males <- 1 - exp(-exp(matrix(rep(input_parameters$death_log_hr_osteoporosis_male, length(lifetables$Males_lograte)) +   rep(lifetables$Males_lograte, each = n_samples), nrow = n_samples)))
  death_probability_osteoporosis$Females <- 1 - exp(-exp(matrix(rep(input_parameters$death_log_hr_osteoporosis_female, length(lifetables$Females_lograte)) +   rep(lifetables$Females_lograte, each = n_samples), nrow = n_samples)))
  colnames(death_probability_osteoporosis$Males) <- colnames(death_probability_osteoporosis$Females) <- lifetables$Age
  
  death_probability_osteoporosis$Overall <- (percentage_male * death_probability_osteoporosis$Males) + ((1-percentage_male) * death_probability_osteoporosis$Females)
  

  # NHL probabilities will change depending on starting age
  # Only need to adjust if starting age is 18 or less
  if(starting_age <= 18) {
    transition_matrices[, 1:(18 - starting_age + 1), "CD GFD no complications", "CD GFD NHL"] <-  
      transition_matrices[, 1:(18 - starting_age + 1),"CD GFD osteoporosis", "CD GFD NHL"] <- NHL_probability_GFD_18orless
    transition_matrices[, (18 - starting_age + 2):n_cycles, "CD GFD no complications", "CD GFD NHL"] <-  
      transition_matrices[, (18 - starting_age + 2):n_cycles,"CD GFD osteoporosis", "CD GFD NHL"] <- NHL_probability_GFD_18plus
    
    transition_matrices[, 1:(18 - starting_age + 1), "Undiagnosed CD no complications", "Undiagnosed CD NHL"] <- 
      transition_matrices[, 1:(18 - starting_age + 1),"Undiagnosed CD osteoporosis", "Undiagnosed CD NHL"] <- (1 - input_parameters$probability_late_diagnosis) * 
      NHL_probability_noGFD_18orless
    transition_matrices[, (18 - starting_age + 2):n_cycles, "Undiagnosed CD no complications", "Undiagnosed CD NHL"] <- 
      transition_matrices[, (18 - starting_age + 2):n_cycles,"Undiagnosed CD osteoporosis", "Undiagnosed CD NHL"] <- (1 - input_parameters$probability_late_diagnosis) * 
      NHL_probability_noGFD_18plus
    
    transition_matrices[, 1:(18 - starting_age + 1), "Undiagnosed CD no complications", "CD GFD NHL"] <- 
      transition_matrices[, 1:(18 - starting_age + 1),"Undiagnosed CD osteoporosis", "CD GFD NHL"] <- input_parameters$probability_late_diagnosis * 
      NHL_probability_noGFD_18orless
    transition_matrices[, (18 - starting_age + 2):n_cycles, "Undiagnosed CD no complications", "CD GFD NHL"] <- 
      transition_matrices[, (18 - starting_age + 2):n_cycles,"Undiagnosed CD osteoporosis", "CD GFD NHL"] <- input_parameters$probability_late_diagnosis * 
      NHL_probability_noGFD_18plus
  } else {
    # Otherwise all cycles are for >18 year olds
    transition_matrices[, , "CD GFD no complications", "CD GFD NHL"] <-  
      transition_matrices[, ,"CD GFD osteoporosis", "CD GFD NHL"] <- NHL_probability_GFD_18plus
    transition_matrices[, , "Undiagnosed CD no complications", "Undiagnosed CD NHL"] <- 
      transition_matrices[, ,"Undiagnosed CD osteoporosis", "Undiagnosed CD NHL"] <- (1 - input_parameters$probability_late_diagnosis) * 
      NHL_probability_noGFD_18plus
    transition_matrices[, , "Undiagnosed CD no complications", "CD GFD NHL"] <- 
      transition_matrices[, ,"Undiagnosed CD osteoporosis", "CD GFD NHL"] <- input_parameters$probability_late_diagnosis * 
      NHL_probability_noGFD_18plus
  }
  transition_matrices[, ,"Undiagnosed CD NHL", "CD GFD NHL"] <- input_parameters$probability_late_diagnosis 
  
  n_agecategories <- (n_cycles/10) - 1
  for(i_age_category in c(0:n_agecategories)) {
    transition_matrices[, (c(1:10) + i_age_category * 10), "CD GFD no complications", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_all[, starting_age_column + i_age_category]
    transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed CD no complications", "Undiagnosed CD osteoporosis"] <- (1 - input_parameters$probability_late_diagnosis) * osteoporosis_probability_noGFD_all[, starting_age_column + i_age_category]
    transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed CD no complications", "CD GFD osteoporosis"] <-  input_parameters$probability_late_diagnosis * osteoporosis_probability_noGFD_all[, starting_age_column + i_age_category]
    transition_matrices[, (c(1:10) + i_age_category * 10),"Undiagnosed CD no complications", "CD GFD no complications"] <- input_parameters$probability_late_diagnosis - 
      (transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed CD no complications", "CD GFD osteoporosis"] +  (transition_matrices[1, (c(1:10) + i_age_category * 10), "Undiagnosed CD no complications", "CD GFD NHL"] ))
    transition_matrices[, ,"Undiagnosed CD osteoporosis", "CD GFD osteoporosis"] <- input_parameters$probability_late_diagnosis - ( transition_matrices[, ,"Undiagnosed CD osteoporosis", "CD GFD NHL"] )
  }
  
  n_ages <- 90 - starting_age
  
  # NHL and osteoporosis death probabilities very high so instead of assuming no competing risks we assume
  # it reduces all other transitions (e.g., the 33% who do not die of NHL are divided among remaining probabilities)
  # This avoids negative probabilities of remaining in the same state
  death_state_index <- which(state_names == "Death")
  for (i_age in 1:n_ages){
    transition_matrices[, i_age, "CD GFD no complications", "Death"] <- death_probability_nocomplications[starting_age + i_age, 2]
    transition_matrices[, i_age, "CD GFD osteoporosis", "Death"] <- death_probability_osteoporosis$Overall[, starting_age + i_age]
    transition_matrices[, i_age, "CD GFD osteoporosis", -death_state_index] <-  transition_matrices[, i_age, "CD GFD osteoporosis", -death_state_index] * (1 - death_probability_osteoporosis$Overall[, starting_age + i_age])
    transition_matrices[, i_age, "Undiagnosed CD no complications", "Death"] <- death_probability_nocomplications[starting_age + i_age, 2]
    transition_matrices[, i_age, "Undiagnosed CD osteoporosis", "Death"] <- death_probability_osteoporosis$Overall[, starting_age + i_age]
    transition_matrices[, i_age, "Undiagnosed CD osteoporosis", -death_state_index] <- transition_matrices[, i_age, "Undiagnosed CD osteoporosis", -death_state_index] * (1 - death_probability_osteoporosis$Overall[, starting_age + i_age])
  }

  # Probabilities of death form NHL
  transition_matrices[, , "CD GFD NHL", "Death"] <- 
    transition_matrices[, , "Undiagnosed CD NHL", "Death"] <- input_parameters$death_probability_NHL
  
  # Multiply other transitions to prevent probabilities exceeding 1  
  for(i_sample in 1:n_samples) {
    transition_matrices[i_sample, , "Undiagnosed CD NHL", -death_state_index] <- 
      transition_matrices[i_sample, , "Undiagnosed CD NHL", -death_state_index] * (1 - input_parameters$death_probability_NHL[i_sample])
    transition_matrices[i_sample, , "CD GFD NHL", -death_state_index] <- 
      transition_matrices[i_sample, , "CD GFD NHL", -death_state_index] * (1 - input_parameters$death_probability_NHL[i_sample])
  }
  
  
  for(i_state in 1:length(state_names)) {
    transition_matrices[, , i_state, i_state] <- 1 - 
      apply(transition_matrices[, , i_state, -i_state], c(1,2), sum, na.rm=TRUE)
  }
  

  transition_matrices[, , , ] [is.na(transition_matrices[, , , ] )] <- 0
  
  # Check that rows sum to 1
  # rowSums (transition_matrices[1, 4, , ], na.rm = FALSE, dims = 1)
  
  return(transition_matrices[, , , ])
}
  

