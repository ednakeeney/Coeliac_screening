generate_transition_matrices <- function(input_parameters) {
 
   starting_age <- ifelse(population == "adults", 50, 10) #based on mean age in under and over 18s in CPRD cost data
  n_cycles <- 90 - starting_age
  
  transition_matrices <- array(dim = c(n_samples, n_cycles, n_states, n_states),
                               dimnames = list(NULL, NULL, state_names, state_names))
  
  # There is one transition matrix for each cycle and each PSA sample
  # Store them in an array with (before filling in below) NA entries
  # Define transition matrices
  
  starting_age_columnandrow <- read.csv("starting_age_column.csv")
  starting_age_columnandrow$starting_age_row <- c(1:11)
  starting_age_column <- starting_age_columnandrow[starting_age_columnandrow$Starting.age == starting_age, 2]
  
  # Osteoporosis probabilities On GFD
  osteoporosis_probability_GFD_all <- data.frame(input_parameters$osteoporosis_probability_GFD_0, input_parameters$osteoporosis_probability_GFD_10, input_parameters$osteoporosis_probability_GFD_20, input_parameters$osteoporosis_probability_GFD_30,
                                                 input_parameters$osteoporosis_probability_GFD_40, input_parameters$osteoporosis_probability_GFD_50, input_parameters$osteoporosis_probability_GFD_60,
                                                 input_parameters$osteoporosis_probability_GFD_70, input_parameters$osteoporosis_probability_GFD_80, input_parameters$osteoporosis_probability_GFD_90)
  
  # Subfertility probabilities on GFD
  subfertility_probability_GFD_all <- data.frame(input_parameters$subfertility_probability_GFD_0, input_parameters$subfertility_probability_GFD_10, input_parameters$subfertility_probability_GFD_20, input_parameters$subfertility_probability_GFD_30,
                                                 input_parameters$subfertility_probability_GFD_40, input_parameters$subfertility_probability_GFD_50, input_parameters$subfertility_probability_GFD_60,
                                                 input_parameters$subfertility_probability_GFD_70, input_parameters$subfertility_probability_GFD_80, input_parameters$subfertility_probability_GFD_90)
  
  # Corresponding probabilities not on GFD
  osteoporosis_probability_noGFD_all <- data.frame(input_parameters$osteoporosis_probability_noGFD_0, input_parameters$osteoporosis_probability_noGFD_10, input_parameters$osteoporosis_probability_noGFD_20, input_parameters$osteoporosis_probability_noGFD_30,
                                                   input_parameters$osteoporosis_probability_noGFD_40, input_parameters$osteoporosis_probability_noGFD_50, input_parameters$osteoporosis_probability_noGFD_60,
                                                   input_parameters$osteoporosis_probability_noGFD_70, input_parameters$osteoporosis_probability_noGFD_80, input_parameters$osteoporosis_probability_noGFD_90)
  
  
  subfertility_probability_noGFD_all <- data.frame(input_parameters$subfertility_probability_noGFD_0, input_parameters$subfertility_probability_noGFD_10, input_parameters$subfertility_probability_noGFD_20, input_parameters$subfertility_probability_noGFD_30,
                                                   input_parameters$subfertility_probability_noGFD_40, input_parameters$subfertility_probability_noGFD_50, input_parameters$subfertility_probability_noGFD_60,
                                                   input_parameters$subfertility_probability_noGFD_70, input_parameters$subfertility_probability_noGFD_80, input_parameters$subfertility_probability_noGFD_90)
  lifetables <- read.csv("lifetables.csv")
  percentage_male <- 0.5
  lifetables$Overall <- (percentage_male * lifetables$Males) + ((1-percentage_male) * lifetables$Females)
  death_probability_nocomplications	<- data.frame(lifetables$Age, lifetables$Overall)
  
  death_probability_osteoporosis <-	lifetables
  death_probability_osteoporosis$Males <- death_probability_osteoporosis$Males*3.5   #Currently  not probabilistic as lifetables are not probabilistic
  death_probability_osteoporosis$Females <- death_probability_osteoporosis$Females*2.4
  death_probability_osteoporosis$Overall <- (percentage_male * death_probability_osteoporosis$Males) + ((1-percentage_male) * death_probability_osteoporosis$Females)
  death_probability_osteoporosis	<- data.frame(death_probability_osteoporosis$Age, death_probability_osteoporosis$Overall)
  
  
  
  transition_matrices[, , "CD GFD no complications", "CD GFD NHL"] <-  
    transition_matrices[, ,"CD GFD subfertility", "CD GFD NHL"] <- 
    transition_matrices[, ,"CD GFD osteoporosis", "CD GFD NHL"] <- input_parameters$NHL_probability_GFD
  transition_matrices[, , "Undiagnosed CD no complications", "Undiagnosed CD NHL"] <- 
    transition_matrices[, ,"Undiagnosed CD subfertility", "Undiagnosed CD NHL"] <- 
    transition_matrices[, ,"Undiagnosed CD osteoporosis", "Undiagnosed CD NHL"] <- (1 - input_parameters$probability_late_diagnosis) * input_parameters$NHL_probability_noGFD
  transition_matrices[, , "Undiagnosed CD no complications", "CD GFD NHL"] <- 
    transition_matrices[, ,"Undiagnosed CD subfertility", "CD GFD NHL"] <- 
    transition_matrices[, ,"Undiagnosed CD osteoporosis", "CD GFD NHL"] <- input_parameters$probability_late_diagnosis * input_parameters$NHL_probability_noGFD
  transition_matrices[, ,"Undiagnosed CD NHL", "CD GFD NHL"] <- input_parameters$probability_late_diagnosis 

  
  n_agecategories <- (n_cycles/10) - 1
  
  for(i_age_category in c(0:n_agecategories)) {
    transition_matrices[, (c(1:10) + i_age_category * 10), "CD GFD no complications", "CD GFD subfertility"] <- subfertility_probability_GFD_all[, starting_age_column + i_age_category]  
    transition_matrices[, (c(1:10) + i_age_category * 10), "CD GFD no complications", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_all[, starting_age_column + i_age_category]
    transition_matrices[, (c(1:10) + i_age_category * 10), "CD GFD subfertility", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_all[, starting_age_column + i_age_category]
    transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed CD no complications", "Undiagnosed CD subfertility"] <- (1 - input_parameters$probability_late_diagnosis) * subfertility_probability_noGFD_all[, starting_age_column + i_age_category]
    transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed CD no complications", "Undiagnosed CD osteoporosis"] <- (1 - input_parameters$probability_late_diagnosis) * osteoporosis_probability_noGFD_all[, starting_age_column + i_age_category]
    transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed CD subfertility", "Undiagnosed CD osteoporosis"] <- (1 - input_parameters$probability_late_diagnosis) * osteoporosis_probability_noGFD_all[, starting_age_column + i_age_category]
    transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed CD subfertility", "CD GFD osteoporosis"] <- input_parameters$probability_late_diagnosis * osteoporosis_probability_noGFD_all[, starting_age_column + i_age_category]
     transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed CD no complications", "CD GFD subfertility"] <- input_parameters$probability_late_diagnosis * subfertility_probability_noGFD_all[, starting_age_column + i_age_category]
    transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed CD no complications", "CD GFD osteoporosis"] <-  input_parameters$probability_late_diagnosis * osteoporosis_probability_noGFD_all[, starting_age_column + i_age_category]
    transition_matrices[, (c(1:10) + i_age_category * 10),"Undiagnosed CD no complications", "CD GFD no complications"] <- input_parameters$probability_late_diagnosis - (transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed CD no complications", "CD GFD subfertility"] +  transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed CD no complications", "CD GFD osteoporosis"] +  ((1 - input_parameters$probability_late_diagnosis) * input_parameters$NHL_probability_noGFD) )
    transition_matrices[, ,"Undiagnosed CD osteoporosis", "CD GFD osteoporosis"] <- input_parameters$probability_late_diagnosis - ( transition_matrices[, ,"Undiagnosed CD osteoporosis", "CD GFD NHL"] )
    transition_matrices[, (c(1:10) + i_age_category * 10),"Undiagnosed CD subfertility", "CD GFD subfertility"] <- input_parameters$probability_late_diagnosis - (transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed CD subfertility", "CD GFD osteoporosis"] + (input_parameters$probability_late_diagnosis * input_parameters$NHL_probability_noGFD) )
      
      
  }
  n_ages <- 90 - starting_age
  
  for (i_age in 1:n_ages){
    transition_matrices[, i_age, "CD GFD no complications", "Death"] <- death_probability_nocomplications[starting_age + i_age, 2]
    transition_matrices[, i_age, "CD GFD subfertility", "Death"] <- death_probability_nocomplications[starting_age + i_age, 2]
    transition_matrices[, i_age, "CD GFD osteoporosis", "Death"] <- death_probability_osteoporosis[starting_age + i_age, 2]
    transition_matrices[, i_age, "Undiagnosed CD no complications", "Death"] <- death_probability_nocomplications[starting_age + i_age, 2]
    transition_matrices[, i_age, "Undiagnosed CD subfertility", "Death"] <- death_probability_nocomplications[starting_age + i_age, 2]
    transition_matrices[, i_age, "Undiagnosed CD osteoporosis", "Death"] <- death_probability_osteoporosis[starting_age + i_age, 2]
  }
  

  
  transition_matrices[, , "CD GFD NHL", "Death"] <- 
    transition_matrices[, , "Undiagnosed CD NHL", "Death"] <- input_parameters$death_probability_NHL
  
  
  for(i_state in 1:length(state_names)) {
    transition_matrices[, , i_state, i_state] <- 1 - 
      apply(transition_matrices[, , i_state, -i_state], c(1,2), sum, na.rm=TRUE)
  }
  

  transition_matrices[, , , ] [is.na(transition_matrices[, , , ] )] <- 0
  
  rowSums (transition_matrices[1, 4, , ], na.rm = FALSE, dims = 1)
  return(transition_matrices[, , , ])
  
}
  

