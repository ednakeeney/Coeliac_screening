generate_state_qalys <- function(input_parameters) {
  
  # State utilities
  # Anything between 0 and 1
  # Now define the QALYS associated with the states per cycle
  # There is one for each PSA sample and each state
  # Store in an NA array and then fill in below
  
  state_qalys <- array(dim=c(n_samples, n_cycles, n_states), dimnames = list(NULL, NULL, state_names))
  
  
  for(i_age_category in c(0:4)) {
    state_qalys[, (c(1:10) + i_age_category * 10), "CD GFD no complications"] <- eq5d_norms[starting_age_row + i_age_category, 2] * utility_GFD
    state_qalys[, (c(1:10) + i_age_category * 10), "Undiagnosed CD no complications"] <- eq5d_norms[starting_age_row + i_age_category, 2] * utility_undiagnosedCD
  }
  
  state_qalys[, , "CD GFD subfertility"] <- state_qalys[, , "CD GFD no complications"] - disutility_subfertility
  state_qalys[, , "CD GFD osteoporosis"] <- state_qalys[, , "CD GFD no complications"] - disutility_osteoporosis 
  state_qalys[, , "CD GFD NHL"] <- state_qalys[, , "CD GFD no complications"] - disutility_NHL
  state_qalys[, , "CD no GFD no complications"] <- state_qalys[, , "CD GFD no complications"] - disutility_noGFD
  state_qalys[, , "CD no GFD subfertility"] <-  state_qalys[, , "CD no GFD no complications"] - disutility_subfertility
  state_qalys[, , "CD no GFD osteoporosis"] <- state_qalys[, , "CD no GFD no complications"] - disutility_osteoporosis
  state_qalys[, , "CD no GFD NHL"] <- state_qalys[, , "CD no GFD no complications"] - disutility_NHL
  state_qalys[, , "Undiagnosed CD subfertility"] <-  state_qalys[, , "Undiagnosed CD no complications"] - disutility_subfertility
  state_qalys[, , "Undiagnosed CD osteoporosis"] <- state_qalys[, , "Undiagnosed CD no complications"] - disutility_osteoporosis
  state_qalys[, , "Undiagnosed CD NHL"] <- state_qalys[, , "Undiagnosed CD no complications"] - disutility_NHL
  state_qalys[, , "Death"] <- 0
  return(state_qalys[, , ])
}