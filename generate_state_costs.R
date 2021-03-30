generate_state_costs <- function(input_parameters) {
  
  state_costs <- array(dim=c(n_samples, n_cycles, n_states), dimnames = list(NULL, NULL, state_names))
  
  starting_age_columnandrow <- read.csv("starting_age_column.csv")
  starting_age_columnandrow$starting_age_row <- c(1:11)
  starting_age_column <- starting_age_columnandrow[starting_age_columnandrow$Starting.age == starting_age, 2]
  
  probability_IDA <- data.frame(input_parameters$probability_IDA_0, input_parameters$probability_IDA_10, input_parameters$probability_IDA_20, input_parameters$probability_IDA_30, input_parameters$probability_IDA_40, 
                                input_parameters$probability_IDA_50, input_parameters$probability_IDA_60, input_parameters$probability_IDA_70, input_parameters$probability_IDA_80, input_parameters$probability_IDA_90)
  
  
  for(i_age_category in c(0:4)) {
    state_costs[, (c(1:10) + i_age_category * 10), "CD GFD no complications"] <- input_parameters$cost_CDGFD + (probability_IDA[, starting_age_column + i_age_category] * input_parameters$cost_IDA) + (transition_matrices[, (c(1:10) + i_age_category * 10), "CD GFD no complications", "CD GFD subfertility"] * input_parameters$cost_subfertility) + (input_parameters$NHL_probability_GFD * input_parameters$cost_NHL)
    state_costs[, (c(1:10) + i_age_category * 10), "CD GFD osteoporosis"]  <- input_parameters$cost_osteoporosis + input_parameters$cost_CDGFD + (probability_IDA[, starting_age_column + i_age_category] * input_parameters$cost_IDA) + (input_parameters$NHL_probability_GFD + input_parameters$cost_NHL)
    state_costs[, (c(1:10) + i_age_category * 10), "Undiagnosed CD no complications"] <- input_parameters$cost_undiagnosedCD + (probability_IDA[, starting_age_column + i_age_category] * input_parameters$cost_IDA) + (input_parameters$probability_late_diagnosis * input_parameters$probability_biopsy * input_parameters$cost_biopsy) +
      (transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed CD no complications", "Undiagnosed CD subfertility"] * input_parameters$cost_subfertility) + (transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed CD no complications", "CD GFD subfertility"] * input_parameters$cost_subfertility) +  (input_parameters$NHL_probability_noGFD * input_parameters$cost_NHL)
    state_costs[, , "Undiagnosed CD osteoporosis"] <- input_parameters$cost_undiagnosedCD + input_parameters$cost_osteoporosis + (probability_IDA[, starting_age_column + i_age_category] * input_parameters$cost_IDA) + (input_parameters$NHL_probability_noGFD * input_parameters$cost_NHL)
 
     }
  
  
  state_costs[, , "Undiagnosed CD subfertility"] <-  input_parameters$cost_undiagnosedCD + (probability_IDA[, starting_age_column + i_age_category] * input_parameters$cost_IDA) + (input_parameters$NHL_probability_noGFD * input_parameters$cost_NHL)
  state_costs[, , "Undiagnosed CD NHL"] <- input_parameters$cost_undiagnosedCD + (probability_IDA[, starting_age_column + i_age_category] * input_parameters$cost_IDA)
  state_costs[, , "CD GFD subfertility"] <- state_costs[, , "CD GFD NHL"] <- 
    input_parameters$cost_CDGFD + (probability_IDA[, starting_age_column + i_age_category] * input_parameters$cost_IDA) 
  state_costs[, , "Death"] <- 0
  return(state_costs[, , ])
}