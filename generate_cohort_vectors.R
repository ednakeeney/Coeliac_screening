generate_cohort_vectors <- function(input_parameters) {
  
  # Build an array to store the cohort vector at each cycle
  # Each cohort vector has n_states elements: probability of being in each state,
  # There is one cohort vector for each treatment, for each PSA sample, for each cycle.
  cohort_vectors<-array(dim=c(n_treatments,n_samples,n_cycles,n_states),  
                        dimnames=list(t_names,NULL,NULL, state_names))
  
  
  for (i_treatment in c(1:3)) { 
    cohort_vectors[i_treatment, , 1, "CD GFD no complications"] <- tp[, i_treatment] * probability_nocomplications * adherence
    cohort_vectors[i_treatment, , 1, "CD GFD subfertility"] <- tp[, i_treatment] * probability_subfertility * adherence
    cohort_vectors[i_treatment, , 1, "CD GFD osteoporosis"] <- tp[, i_treatment] * probability_osteoporosis * adherence
    cohort_vectors[i_treatment, , 1,"CD GFD NHL"] <- tp[, i_treatment] * probability_NHL * adherence
    
    cohort_vectors[i_treatment, , 1,"CD no GFD no complications"] <- tp[, i_treatment] * probability_nocomplications * (1-adherence)
    cohort_vectors[i_treatment, , 1,"CD no GFD subfertility"] <- tp[, i_treatment] * probability_subfertility * (1-adherence)
    cohort_vectors[i_treatment, , 1,"CD no GFD osteoporosis"] <- tp[, i_treatment] * probability_osteoporosis * (1-adherence)
    cohort_vectors[i_treatment, , 1,"CD no GFD NHL"] <- tp[, i_treatment] * probability_NHL * (1-adherence)
    
    cohort_vectors[i_treatment, , 1,"Undiagnosed CD no complications"] <- fn[, i_treatment] * probability_nocomplications 
    cohort_vectors[i_treatment, , 1,"Undiagnosed CD subfertility"] <- fn[, i_treatment] * probability_subfertility 
    cohort_vectors[i_treatment, , 1,"Undiagnosed CD osteoporosis"] <- fn[, i_treatment] * probability_osteoporosis 
    cohort_vectors[i_treatment, , 1,"Undiagnosed CD NHL"] <- fn[, i_treatment] * probability_NHL 
  }
  
  cohort_vectors[, , , ] [is.na(cohort_vectors[, , , ] )] <- 0
  return(cohort_vectors[, , ,])
}