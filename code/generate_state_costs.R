generate_state_costs <- function(input_parameters) {
  n_samples <- dim(input_parameters)[1]
  
  starting_age <- ifelse(population == "men" | population == "women", 18, 10) #based on mean age in under and over 18s in CPRD cost data
  n_cycles <- 90 - starting_age
  state_costs <- array(dim=c(n_samples, n_cycles, n_states), dimnames = list(NULL, NULL, state_names))
  
  starting_age_columnandrow <- read.csv("data/starting_age_column.csv")
  starting_age_columnandrow$starting_age_row <- c(1:11)
  #starting_age_column <- starting_age_columnandrow[starting_age_columnandrow$Starting.age == starting_age, 2]
  #Initial cohort at diagnosis - depends on age at diagnosis
  # e.g. 10 year old starts with 10-20 prevalence and 18 year old starts with 10-20 prevalence. 
  starting_age_column <- which(starting_age_columnandrow$Starting.age  > starting_age)[1] - 1
  
  
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
  
  
  
  probability_IDA <- data.frame(input_parameters$probability_IDA_0, input_parameters$probability_IDA_10, input_parameters$probability_IDA_20, input_parameters$probability_IDA_30, input_parameters$probability_IDA_40, 
                                input_parameters$probability_IDA_50, input_parameters$probability_IDA_60, input_parameters$probability_IDA_70, input_parameters$probability_IDA_80, input_parameters$probability_IDA_90)
  
  n_agecategories <- ceiling((n_cycles/10) - 1)
  
  for(i_age_category in c(0:n_agecategories)) {
    age_category_indices <- (c(1:10) + i_age_category * 10)
    age_category_indices <- age_category_indices[age_category_indices <= n_cycles]
    
    state_costs[, age_category_indices, "CD GFD no complications"] <- input_parameters$cost_CDGFD + (probability_IDA[, starting_age_column + i_age_category] * input_parameters$cost_IDA) + (NHL_probability_GFD * input_parameters$cost_NHL)
    state_costs[, age_category_indices, "CD GFD osteoporosis"]  <- input_parameters$cost_osteoporosis + input_parameters$cost_CDGFD + (probability_IDA[, starting_age_column + i_age_category] * input_parameters$cost_IDA) + (NHL_probability_GFD * input_parameters$cost_NHL)
    state_costs[, age_category_indices, "Undiagnosed CD no complications"] <- input_parameters$cost_undiagnosedCD + (probability_IDA[, starting_age_column + i_age_category] * input_parameters$cost_IDA) + 
      (input_parameters$probability_late_diagnosis * input_parameters$probability_biopsy * input_parameters$cost_biopsy) +
     (NHL_probability_noGFD * input_parameters$cost_NHL)
    state_costs[, age_category_indices, "Undiagnosed CD osteoporosis"] <- input_parameters$cost_undiagnosedCD + input_parameters$cost_osteoporosis + 
      (probability_IDA[, starting_age_column + i_age_category] * input_parameters$cost_IDA) + (NHL_probability_noGFD * input_parameters$cost_NHL)
 
     }
  
  
  state_costs[, , "Undiagnosed CD NHL"] <- input_parameters$cost_undiagnosedCD + (probability_IDA[, starting_age_column + i_age_category] * input_parameters$cost_IDA)
  state_costs[, , "CD GFD NHL"] <- 
    input_parameters$cost_CDGFD + (probability_IDA[, starting_age_column + i_age_category] * input_parameters$cost_IDA) 
  state_costs[, , "Death"] <- 0
  return(state_costs[, , ])
}