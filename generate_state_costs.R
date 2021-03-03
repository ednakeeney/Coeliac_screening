generate_state_costs <- function(input_parameters) {
  
  state_costs <- array(dim=c(n_samples, n_cycles, n_states), dimnames = list(NULL, NULL, state_names))
  
  for(i_age_category in c(0:4)) {
    state_costs[, (c(1:10) + i_age_category * 10), "CD GFD no complications"] <- cost_CDGFD + (probability_IDA[, starting_age_column + i_age_category] * cost_IDA) + (transition_matrices[, (c(1:10) + i_age_category * 10), "CD GFD no complications", "CD GFD subfertility"] * cost_subfertility)
    state_costs[, (c(1:10) + i_age_category * 10), "CD no GFD no complications"] <- cost_CDGFD + (probability_IDA[, starting_age_column + i_age_category] * cost_IDA) +  (transition_matrices[, (c(1:10) + i_age_category * 10), "CD no GFD no complications", "CD no GFD subfertility"] * cost_subfertility)
    state_costs[, (c(1:10) + i_age_category * 10), "CD GFD osteoporosis"] <- state_costs[, (c(1:10) + i_age_category * 10), "CD no GFD osteoporosis"] <- cost_osteoporosis + cost_CDGFD + (probability_IDA[, starting_age_column + i_age_category] * cost_IDA) 
    state_costs[, (c(1:10) + i_age_category * 10), "Undiagnosed CD no complications"] <- cost_undiagnosedCD + (probability_IDA[, starting_age_column + i_age_category] * cost_IDA) + (probability_late_diagnosis * probability_biopsy * cost_biopsy) +
      (transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed CD no complications", "Undiagnosed CD subfertility"] * cost_subfertility) + (transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed CD no complications", "CD GFD subfertility"] * cost_subfertility) + 
      (transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed CD no complications", "CD no GFD subfertility"] * cost_subfertility)
    state_costs[, , "Undiagnosed CD osteoporosis"] <- cost_undiagnosedCD + cost_osteoporosis + (probability_IDA[, starting_age_column + i_age_category] * cost_IDA) 
  }
  
  
  state_costs[, , "Undiagnosed CD subfertility"] <-  state_costs[, , "Undiagnosed CD NHL"] <- cost_undiagnosedCD + (probability_IDA[, starting_age_column + i_age_category] * cost_IDA)
  state_costs[, , "CD GFD subfertility"] <- state_costs[, , "CD GFD NHL"] <-   state_costs[, , "CD no GFD subfertility"] <- state_costs[, , "CD no GFD NHL"] <- 
    cost_CDGFD + (probability_IDA[, starting_age_column + i_age_category] * cost_IDA) 
  state_costs[, , "Death"] <- 0
  return(state_costs[, , ])
}