generate_transition_matrices <- function(input_parameters) {
 
  transition_matrices <- array(dim = c(n_samples, n_cycles, n_states, n_states),
                               dimnames = list(NULL, NULL, state_names, state_names))
  
  # There is one transition matrix for each cycle and each PSA sample
  # Store them in an array with (before filling in below) NA entries
  # Define transition matrices
  
  transition_matrices[, , "CD GFD no complications", "CD GFD NHL"] <-  
    transition_matrices[, ,"CD GFD subfertility", "CD GFD NHL"] <- 
    transition_matrices[, ,"CD GFD osteoporosis", "CD GFD NHL"] <- NHL_probability_GFD
  #  transition_matrices[, , "CD no GFD no complications", "CD no GFD NHL"] <- 
  # transition_matrices[, ,"CD no GFD subfertility", "CD no GFD NHL"] <- 
  #transition_matrices[, ,"CD no GFD osteoporosis", "CD no GFD NHL"] <- 
  transition_matrices[, , "Undiagnosed CD no complications", "Undiagnosed CD NHL"] <- 
    transition_matrices[, ,"Undiagnosed CD subfertility", "Undiagnosed CD NHL"] <- 
    transition_matrices[, ,"Undiagnosed CD osteoporosis", "Undiagnosed CD NHL"] <- (1 - probability_late_diagnosis) * NHL_probability_noGFD
  transition_matrices[, , "Undiagnosed CD no complications", "CD GFD NHL"] <- 
    transition_matrices[, ,"Undiagnosed CD subfertility", "CD GFD NHL"] <- 
    transition_matrices[, ,"Undiagnosed CD osteoporosis", "CD GFD NHL"] <- probability_late_diagnosis * NHL_probability_noGFD
  transition_matrices[, ,"Undiagnosed CD NHL", "CD GFD NHL"] <- probability_late_diagnosis 
  #transition_matrices[, , "Undiagnosed CD no complications", "CD no GFD NHL"] <- (1 - adherence) * probability_late_diagnosis * NHL_probability_noGFD
  
  
  for(i_age_category in c(0:4)) {
    transition_matrices[, (c(1:10) + i_age_category * 10), "CD GFD no complications", "CD GFD subfertility"] <- subfertility_probability_GFD_all[, starting_age_column + i_age_category]  
    transition_matrices[, (c(1:10) + i_age_category * 10), "CD GFD no complications", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_all[, starting_age_column + i_age_category]
    transition_matrices[, (c(1:10) + i_age_category * 10), "CD GFD subfertility", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_all[, starting_age_column + i_age_category]
    #transition_matrices[, (c(1:10) + i_age_category * 10), "CD no GFD no complications", "CD no GFD subfertility"] <- subfertility_probability_noGFD_all[, starting_age_column + i_age_category]
    #transition_matrices[, (c(1:10) + i_age_category * 10), "CD no GFD no complications", "CD no GFD osteoporosis"] <- osteoporosis_probability_noGFD_all[, starting_age_column + i_age_category]
    #transition_matrices[, (c(1:10) + i_age_category * 10), "CD no GFD subfertility", "CD no GFD osteoporosis"] <- osteoporosis_probability_noGFD_all[, starting_age_column + i_age_category]
    transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed CD no complications", "Undiagnosed CD subfertility"] <- (1 - probability_late_diagnosis) * subfertility_probability_noGFD_all[, starting_age_column + i_age_category]
    transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed CD no complications", "Undiagnosed CD osteoporosis"] <- (1 - probability_late_diagnosis) * osteoporosis_probability_noGFD_all[, starting_age_column + i_age_category]
    transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed CD subfertility", "Undiagnosed CD osteoporosis"] <- (1 - probability_late_diagnosis) * osteoporosis_probability_noGFD_all[, starting_age_column + i_age_category]
    transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed CD subfertility", "CD GFD osteoporosis"] <- probability_late_diagnosis * osteoporosis_probability_noGFD_all[, starting_age_column + i_age_category]
     transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed CD no complications", "CD GFD subfertility"] <- probability_late_diagnosis * subfertility_probability_noGFD_all[, starting_age_column + i_age_category]
    transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed CD no complications", "CD GFD osteoporosis"] <-  probability_late_diagnosis * osteoporosis_probability_noGFD_all[, starting_age_column + i_age_category]
   # transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed CD no complications", "CD no GFD subfertility"] <- (1 - adherence) * probability_late_diagnosis * subfertility_probability_noGFD_all[, starting_age_column + i_age_category]
  #  transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed CD no complications", "CD no GFD osteoporosis"] <- (1 - adherence) * probability_late_diagnosis * osteoporosis_probability_noGFD_all[, starting_age_column + i_age_category]
    transition_matrices[, (c(1:10) + i_age_category * 10),"Undiagnosed CD no complications", "CD GFD no complications"] <- probability_late_diagnosis - (transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed CD no complications", "CD GFD subfertility"] +  transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed CD no complications", "CD GFD osteoporosis"] +  ((1 - probability_late_diagnosis) * NHL_probability_noGFD) )
    transition_matrices[, ,"Undiagnosed CD osteoporosis", "CD GFD osteoporosis"] <- probability_late_diagnosis - ( transition_matrices[, ,"Undiagnosed CD osteoporosis", "CD GFD NHL"] )
    transition_matrices[, (c(1:10) + i_age_category * 10),"Undiagnosed CD subfertility", "CD GFD subfertility"] <- probability_late_diagnosis - (transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed CD subfertility", "CD GFD osteoporosis"] + (probability_late_diagnosis * NHL_probability_noGFD) )
      
      
     }
  
  for (i_age in 1:50){
    transition_matrices[, i_age, "CD GFD no complications", "Death"] <- death_probability_nocomplications[starting_age + i_age, 2]
    transition_matrices[, i_age, "CD GFD subfertility", "Death"] <- death_probability_nocomplications[starting_age + i_age, 2]
    transition_matrices[, i_age, "CD GFD osteoporosis", "Death"] <- death_probability_osteoporosis[starting_age + i_age, 2]
   # transition_matrices[, i_age, "CD no GFD no complications", "Death"] <- death_probability_nocomplications[starting_age + i_age, 2]
  #  transition_matrices[, i_age, "CD no GFD subfertility", "Death"] <- death_probability_nocomplications[starting_age + i_age, 2]
   # transition_matrices[, i_age, "CD no GFD osteoporosis", "Death"] <- death_probability_osteoporosis[starting_age + i_age, 2]
    transition_matrices[, i_age, "Undiagnosed CD no complications", "Death"] <- death_probability_nocomplications[starting_age + i_age, 2]
    transition_matrices[, i_age, "Undiagnosed CD subfertility", "Death"] <- death_probability_nocomplications[starting_age + i_age, 2]
    transition_matrices[, i_age, "Undiagnosed CD osteoporosis", "Death"] <- death_probability_osteoporosis[starting_age + i_age, 2]
  }
  

  
  transition_matrices[, , "CD GFD NHL", "Death"] <- 
   # transition_matrices[, , "CD no GFD NHL", "Death"] <- 
    transition_matrices[, , "Undiagnosed CD NHL", "Death"] <- death_probability_NHL
  
  
 # transition_matrices[, ,"Undiagnosed CD no complications", "CD no GFD no complications"] <- 
  #  transition_matrices[, ,"Undiagnosed CD subfertility", "CD no GFD subfertility"] <- 
   # transition_matrices[, ,"Undiagnosed CD osteoporosis", "CD no GFD osteoporosis"] <- 
    #transition_matrices[, ,"Undiagnosed CD NHL", "CD no GFD NHL"] <- (1-adherence) * probability_late_diagnosis
  
  
  for(i_state in 1:length(state_names)) {
    transition_matrices[, , i_state, i_state] <- 1 - 
      apply(transition_matrices[, , i_state, -i_state], c(1,2), sum, na.rm=TRUE)
  }
  
  # HT: This is a little dangerous as you may have missed something by accident.
  # Please look at one matrix (e.g. transition_matrices[1, 1, , ]) and check that each NA really should be NA. We could perhaps do this together on our next call.
  # EK: Yes, let's talk through it
  transition_matrices[, , , ] [is.na(transition_matrices[, , , ] )] <- 0
  
  rowSums (transition_matrices[1, 4, , ], na.rm = FALSE, dims = 1)
  return(transition_matrices[, , , ])
  
}
  

View(transition_matrices[1, 20, , ])
