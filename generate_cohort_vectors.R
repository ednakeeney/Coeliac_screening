generate_cohort_vectors <- function(input_parameters) {
  
  # Build an array to store the cohort vector at each cycle
  # Each cohort vector has n_states elements: probability of being in each state,
  # There is one cohort vector for each treatment, for each PSA sample, for each cycle.
  cohort_vectors<-array(dim=c(n_treatments,n_samples,n_cycles,n_states),  
                        dimnames=list(t_names,NULL,NULL, state_names))
  
  tp <- fn <- fp <- tn <- matrix(n=n_samples, ncol=n_treatments)
  
  # Probabilities for test 
  tp[,1] <- (n_samples * p_cd * input_parameters$sens_IgAEMA_adults)/n_samples
  fn[,1] <- 1 - tp[,1]  
  tn[,1] <- (n_samples * p_cd * input_parameters$spec_IgAEMA_adults)/n_samples
  fp[,1] <- 1 - tn[,1]
  
  # Probabilities for test + biopsy
  tp[,2] <- (n_samples * p_cd * input_parameters$sens_IgATTG)/n_samples
  fn[,2] <- 1 - tp[,2] 
  tn[,2] <- (n_samples * p_cd * input_parameters$spec_IgATTG)/n_samples
  fp[,2] <- 1 - tn[,2]
  
  # Probabilities for Double test
  tp[,3] <- (n_samples * p_cd * input_parameters$sens_doubletest)/n_samples
  fn[,3] <- 1 - tp[,3] 
  tn[,3] <- (n_samples * p_cd * input_parameters$spec_doubletest)/n_samples
  fp[,3] <- 1 - tn[,3]
  
  
  for (i_treatment in c(1:3)) { 
    cohort_vectors[i_treatment, , 1, "CD GFD no complications"] <- tp[, i_treatment] * input_parameters$probability_nocomplications 
    cohort_vectors[i_treatment, , 1, "CD GFD subfertility"] <- tp[, i_treatment] * input_parameters$probability_subfertility
    cohort_vectors[i_treatment, , 1, "CD GFD osteoporosis"] <- tp[, i_treatment] * input_parameters$probability_osteoporosis 
    cohort_vectors[i_treatment, , 1,"CD GFD NHL"] <- tp[, i_treatment] * input_parameters$probability_NHL
    
   # cohort_vectors[i_treatment, , 1,"CD no GFD no complications"] <- tp[, i_treatment] * probability_nocomplications * (1-adherence)
  #  cohort_vectors[i_treatment, , 1,"CD no GFD subfertility"] <- tp[, i_treatment] * probability_subfertility * (1-adherence)
  #  cohort_vectors[i_treatment, , 1,"CD no GFD osteoporosis"] <- tp[, i_treatment] * probability_osteoporosis * (1-adherence)
  #  cohort_vectors[i_treatment, , 1,"CD no GFD NHL"] <- tp[, i_treatment] * probability_NHL * (1-adherence)
    
    cohort_vectors[i_treatment, , 1,"Undiagnosed CD no complications"] <- fn[, i_treatment] * input_parameters$probability_nocomplications 
    cohort_vectors[i_treatment, , 1,"Undiagnosed CD subfertility"] <- fn[, i_treatment] * input_parameters$probability_subfertility 
    cohort_vectors[i_treatment, , 1,"Undiagnosed CD osteoporosis"] <- fn[, i_treatment] * input_parameters$probability_osteoporosis 
    cohort_vectors[i_treatment, , 1,"Undiagnosed CD NHL"] <- fn[, i_treatment] * input_parameters$probability_NHL 
  }
  
  cohort_vectors[, , , ] [is.na(cohort_vectors[, , , ] )] <- 0
  return(cohort_vectors[, , ,])
}