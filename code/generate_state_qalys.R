generate_state_qalys <- function(input_parameters) {
  n_samples <- dim(input_parameters)[1]
  
  starting_age <- ifelse(population == "men" | population == "women", 18, 10) #based on mean age in under and over 18s in CPRD cost data
  n_cycles <- 90 - starting_age
  # State utilities
  # Anything between 0 and 1
  # Now define the QALYS associated with the states per cycle
  # There is one for each PSA sample and each state
  # Store in an NA array and then fill in below
  
  state_qalys <- array(dim=c(n_samples, n_cycles, n_states), dimnames = list(NULL, NULL, state_names))
  
  starting_age_columnandrow <- read.csv("data/starting_age_column.csv")
  starting_age_columnandrow$starting_age_row <- c(1:11)
  #starting_age_row <- starting_age_columnandrow[starting_age_columnandrow$Starting.age == starting_age, 2]
  #Initial cohort at diagnosis - depends on age at diagnosis
  # e.g. 10 year old starts with 10-20 prevalence and 18 year old starts with 10-20 prevalence. 
  starting_age_row <- which(starting_age_columnandrow$Starting.age  > starting_age)[1] - 1
  
  eq5d_norms <- read.csv("data/eq5d_norms.csv")
  eq5d_norms$age <- c(0, 10, 20, 30, 40, 50, 60, 70, 80)
  

  n_agecategories <- ceiling((n_cycles/10) - 1)
  
  for(i_age_category in c(0:n_agecategories)) {
    age_category_indices <- (c(1:10) + i_age_category * 10)
    age_category_indices <- age_category_indices[age_category_indices <= n_cycles]
    
    state_qalys[, age_category_indices, "CD GFD no complications"] <- eq5d_norms[starting_age_row + i_age_category, 2] * input_parameters$utility_GFD
    state_qalys[, age_category_indices, "Undiagnosed CD no complications"] <- eq5d_norms[starting_age_row + i_age_category, 2] * input_parameters$utility_undiagnosedCD - (input_parameters$probability_late_diagnosis * input_parameters$probability_biopsy * input_parameters$disutility_biopsy)
  }
  
  state_qalys[, , "CD GFD osteoporosis"] <- state_qalys[, , "CD GFD no complications"] - input_parameters$disutility_osteoporosis 
  state_qalys[, , "CD GFD NHL"] <- state_qalys[, , "CD GFD no complications"] - input_parameters$disutility_NHL
  state_qalys[, , "Undiagnosed CD osteoporosis"] <- state_qalys[, , "Undiagnosed CD no complications"] - input_parameters$disutility_osteoporosis
  state_qalys[, , "Undiagnosed CD NHL"] <- state_qalys[, , "Undiagnosed CD no complications"] - input_parameters$disutility_NHL
  state_qalys[, , "Death"] <- 0
  return(state_qalys[, , ])
}