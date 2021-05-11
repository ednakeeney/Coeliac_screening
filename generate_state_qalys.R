generate_state_qalys <- function(input_parameters) {
  starting_age <- ifelse(population == "adults", 50, 10) #based on mean age in under and over 18s in CPRD cost data
  n_cycles <- 90 - starting_age
  # State utilities
  # Anything between 0 and 1
  # Now define the QALYS associated with the states per cycle
  # There is one for each PSA sample and each state
  # Store in an NA array and then fill in below
  
  state_qalys <- array(dim=c(n_samples, n_cycles, n_states), dimnames = list(NULL, NULL, state_names))
  
  starting_age_columnandrow <- read.csv("starting_age_column.csv")
  starting_age_columnandrow$starting_age_row <- c(1:11)
  starting_age_row <- starting_age_columnandrow[starting_age_columnandrow$Starting.age == starting_age, 2]
  
  eq5d_norms <- read.csv("eq5d_norms.csv")
  eq5d_norms$age <- c(0, 10, 20, 30, 40, 50, 60, 70, 80)
  
  n_agecategories <- (n_cycles/10) - 1
  
  for(i_age_category in c(0:n_agecategories)) {
    state_qalys[, (c(1:10) + i_age_category * 10), "CD GFD no complications"] <- eq5d_norms[starting_age_row + i_age_category, 2] * input_parameters$utility_GFD
    state_qalys[, (c(1:10) + i_age_category * 10), "Undiagnosed CD no complications"] <- eq5d_norms[starting_age_row + i_age_category, 2] * input_parameters$utility_undiagnosedCD
  }
  
  state_qalys[, , "CD GFD subfertility"] <- state_qalys[, , "CD GFD no complications"] - input_parameters$disutility_subfertility
  state_qalys[, , "CD GFD osteoporosis"] <- state_qalys[, , "CD GFD no complications"] - input_parameters$disutility_osteoporosis 
  state_qalys[, , "CD GFD NHL"] <- state_qalys[, , "CD GFD no complications"] - input_parameters$disutility_NHL
  state_qalys[, , "Undiagnosed CD subfertility"] <-  state_qalys[, , "Undiagnosed CD no complications"] - input_parameters$disutility_subfertility
  state_qalys[, , "Undiagnosed CD osteoporosis"] <- state_qalys[, , "Undiagnosed CD no complications"] - input_parameters$disutility_osteoporosis
  state_qalys[, , "Undiagnosed CD NHL"] <- state_qalys[, , "Undiagnosed CD no complications"] - input_parameters$disutility_NHL
  state_qalys[, , "Death"] <- 0
  return(state_qalys[, , ])
}