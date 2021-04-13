generate_net_benefit <- function(input_parameters) {
  
  transition_matrices <- generate_transition_matrices(input_parameters)
  
  state_qalys <- generate_state_qalys(input_parameters)
  
  state_costs <- generate_state_costs(input_parameters)
  
  treatment_costs <- array(dim=c(n_treatments,n_samples),dimnames=list(t_names,NULL))
  
  post_test_probability_IgAEMA_adults <- data.frame(input_parameters$X0.05, input_parameters$X0.15, input_parameters$X0.25, input_parameters$X0.35, input_parameters$X0.45)
  post_test_probability_IgATTGplusEMA <- data.frame(input_parameters$X0.05.1, input_parameters$X0.15.1, input_parameters$X0.25.1, input_parameters$X0.35.1, input_parameters$X0.45.1)
  
  for (i in 1:5) {
  treatment_costs[i, ] <-  ifelse(post_test_probability_IgAEMA_adults[i] >= 0.9, input_parameters$treatment_cost_IgAEMA, 
                                             input_parameters$treatment_cost_IgAEMA +  input_parameters$cost_biopsy)
  }

  for (i in 1:5) {
  treatment_costs[i+5, ] <-  ifelse(post_test_probability_IgATTGplusEMA[i] >= 0.9, input_parameters$treatment_cost_IgATTG, 
                                              input_parameters$treatment_cost_IgATTG +  input_parameters$cost_biopsy)
  } 

  treatment_costs[11:15, ] <- input_parameters$treatment_cost_doubletest
 
  
  tp <- fn <- fp <- tn <- matrix(n=n_samples, ncol=n_treatments)
  
  # Probabilities for IgAEMA plus biopsy 
  for (i in 1:5) {
  tp[,i] <- ifelse(post_test_probability_IgAEMA_adults[i] >= 0.9, (n_samples * pre_test_probability[i] * input_parameters$sens_IgAEMA_adults)/n_samples, 
                   (n_samples * pre_test_probability[i] * input_parameters$sens_biopsy)/n_samples)
  fn[,i] <- pre_test_probability[i] - tp[,i]  
  tn[,i] <- ifelse(post_test_probability_IgAEMA_adults[i] >= 0.9, (n_samples * (1 - pre_test_probability[i]) * input_parameters$spec_IgAEMA_adults)/n_samples,
                   (n_samples * (1 - pre_test_probability[i]) * input_parameters$spec_biopsy)/n_samples)
  fp[,i] <- (1 - pre_test_probability[i]) - tn[,i]
  }
  
  # Probabilities for IgATTG + IgAEMA

    for (i in 1:5) {
  tp[,i+5] <- ifelse(post_test_probability_IgATTGplusEMA[i] >= 0.9, (n_samples * pre_test_probability[i] * input_parameters$sens_IgATTGplusEMA)/n_samples,
                   (n_samples * pre_test_probability[i] * input_parameters$sens_biopsy)/n_samples)
  fn[,i+5] <- pre_test_probability[i] - tp[,i+5] 
  tn[,i+5] <- ifelse(post_test_probability_IgATTGplusEMA[i] >= 0.9, (n_samples * (1 - pre_test_probability[i]) * input_parameters$spec_IgATTGplusEMA)/n_samples,
                   (n_samples * (1 - pre_test_probability[i]) * input_parameters$spec_biopsy)/n_samples)
  fp[,i+5] <- (1 - pre_test_probability[i]) - tn[,i+5]
    }
  
  # Probabilities for Double test
 
    for (i in 1:5) {
      tp[,i+10] <- (n_samples * pre_test_probability[i] * input_parameters$sens_doubletest)/n_samples
      fn[,i+10] <- pre_test_probability[i] - tp[,i+10] 
      tn[,i+10] <- (n_samples * (1 - pre_test_probability[i]) * input_parameters$spec_doubletest)/n_samples
      fp[,i+10] <- (1 - pre_test_probability[i]) - tn[,i+10]
    }
  
  fp_costs <- array(dim = c(n_treatments, n_samples),
                    dimnames = list(t_names, NULL))
  
  fp_costs[,] <- (fp[,] * input_parameters$cost_gfp) 
  
  diagnosis_costs <- array(dim = c(n_treatments, n_samples),
                           dimnames = list(t_names, NULL))
  
  diagnosis_costs[,] <- (tp[,] + fp[,]) * input_parameters$cost_diagnosis
  
 
  
  # Build an array to store the cohort vector at each cycle
  # Each cohort vector has n_states elements: probability of being in each state,
  # There is one cohort vector for each treatment, for each PSA sample, for each cycle.
  cohort_vectors<-array(dim=c(n_treatments,n_samples,n_cycles,n_states),  
                        dimnames=list(t_names,NULL,NULL, state_names))
  
  #scaling up true positives and false negatives 
  tp <- tp/(tp+fn)
  fn <- 1 - tp
  
  for (i_treatment in c(1:n_treatments)) { 
    cohort_vectors[i_treatment, , 1, "CD GFD no complications"] <- tp[, i_treatment] * input_parameters$probability_nocomplications 
    cohort_vectors[i_treatment, , 1, "CD GFD subfertility"] <- tp[, i_treatment] * input_parameters$probability_subfertility
    cohort_vectors[i_treatment, , 1, "CD GFD osteoporosis"] <- tp[, i_treatment] * input_parameters$probability_osteoporosis 
    cohort_vectors[i_treatment, , 1, "CD GFD NHL"] <- tp[, i_treatment] * input_parameters$probability_NHL
    
    cohort_vectors[i_treatment, , 1, "Undiagnosed CD no complications"] <- fn[, i_treatment] * input_parameters$probability_nocomplications 
    cohort_vectors[i_treatment, , 1, "Undiagnosed CD subfertility"] <- fn[, i_treatment] * input_parameters$probability_subfertility 
    cohort_vectors[i_treatment, , 1, "Undiagnosed CD osteoporosis"] <- fn[, i_treatment] * input_parameters$probability_osteoporosis 
    cohort_vectors[i_treatment, , 1, "Undiagnosed CD NHL"] <- fn[, i_treatment] * input_parameters$probability_NHL 
  }
  
  cohort_vectors[, , , ] [is.na(cohort_vectors[, , , ] )] <- 0
  rowSums (cohort_vectors[, 2, , ], na.rm = FALSE, dims = 1)
  
  # Build an array to store the costs and QALYs accrued per cycle
  # One for each treatment, for each PSA sample, for each cycle
  # These will be filled in below in the main model code
  # Then discounted and summed to contribute to total costs and total QALYs
  cycle_costs <- array(dim = c(n_treatments, n_samples, n_cycles),
                       dimnames = list(t_names, NULL, NULL))
  cycle_qalys <- array(dim = c(n_treatments, n_samples, n_cycles),
                       dimnames = list(t_names, NULL, NULL))
  
  # Build arrays to store the total costs and total QALYs
  # There is one for each treatment and each PSA sample
  # These are filled in below using cycle_costs, 
  # treatment_costs, and cycle_qalys
  total_costs <- array(dim = c(n_treatments, n_samples),
                       dimnames = list(t_names, NULL))
  total_qalys <- array(dim = c(n_treatments, n_samples),
                       dimnames = list(t_names, NULL))
  
  
  disc_vec <- (1/1.035) ^ rep(c(0:(n_cycles/2-1)), each = 2)
  
  # Main model code
  # Loop over the treatment options
  
  for (i_treatment in 1:n_treatments)
  {
    # Loop over the PSA samples
    for(i_sample in 1:n_samples)
    {
      
      transition_matrices_sample <- transition_matrices[i_sample, , , ]
      # Loop over the cycles
      # Cycle 1 is already defined so only need to update cycles 2:n_cycles
      for(i_cycle in 2:n_cycles)
      {
        # Markov update
        # Multiply previous cycle's cohort vector by transition matrix
        # i.e. pi_j = pi_(j-1)*P
        cohort_vectors[i_treatment, i_sample, i_cycle, ] <-
          cohort_vectors[i_treatment, i_sample, i_cycle-1, ] %*%
          transition_matrices_sample[i_cycle, , ]
      }
      
      cohort_vectors_tr_sample <- cohort_vectors[i_treatment, i_sample, , ]
      
      for (i_cycle in 1:n_cycles)
      {
        # Now use the cohort vectors to calculate the 
        # total costs for each cycle
        cycle_costs[i_treatment, i_sample, ] <-
          cohort_vectors_tr_sample %*% state_costs[i_sample, i_cycle, ]
      }
      # And total QALYs for each cycle
      cycle_qalys[i_treatment, i_sample, ] <-
        cohort_vectors_tr_sample %*% state_qalys[i_sample, i_cycle, ]
      
      # Combine the cycle_costs and other costs to get total costs
      # Apply the discount factor 
      # (1 in first year, 1.035 in second, 1.035^2 in third, and so on)
      # Each year acounts for two cycles so need to repeat the discount values
      total_costs[i_treatment, i_sample] <- treatment_costs[i_treatment, i_sample] + fp_costs[i_treatment, i_sample] + diagnosis_costs[i_treatment, i_sample]
      + cycle_costs[i_treatment, i_sample, ] %*% disc_vec
      
      # Combine the cycle_qalys to get total qalys
      # Apply the discount factor 
      # (1 in first year, 1.035 in second, 1.035^2 in third, and so on)
      # Each year acounts for two cycles so need to repeat the discount values
      total_qalys[i_treatment, i_sample] <- cycle_qalys[i_treatment, i_sample, ] %*%
        disc_vec
    }
    
  }
  

  #############################################################################
  ## Analysis of results ######################################################
  #############################################################################
  output <- list()
  
  output$total_costs <- total_costs
  output$total_qalys <- total_qalys
  # Average costs
  output$average_costs <- rowMeans(total_costs)
  # Average effects (in QALY units)
  output$average_effects <- rowMeans(total_qalys)
  
  
  #output$incremental_costs_IgATTGplusIgAEMA_IgAEMA <- total_costs["IgATTGplusIgAEMA", ] - total_costs["IgAEMA", ]
  #output$incremental_effects_IgATTGplusIgAEMA_IgAEMA <- total_qalys["IgATTGplusIgAEMA", ] - total_qalys["IgAEMA", ]
  
  #output$incremental_costs_doubletest_IgAEMA <- total_costs["Double test", ] - total_costs["IgAEMA", ]
  #output$incremental_effects_doubletest_IgAEMA <- total_qalys["Double test", ] - total_qalys["IgAEMA", ]
  
  
  #output$ICER_IgATTGplusIgAEMA_IgAEMA <- mean(output$incremental_costs_IgATTGplusIgAEMA_IgAEMA)/mean(output$incremental_effects_IgATTGplusIgAEMA_IgAEMA)
  #output$ICER_doubletest_IgAEMA <- mean(output$incremental_costs_doubletest_IgAEMA)/mean(output$incremental_effects_doubletest_IgAEMA)
  
  # Incremental net benefit at the ?20,000 willingness-to-pay
  
  #output$incremental_net_benefit_IgATTGplusIgAEMA_IgAEMA <- 20000*output$incremental_effects_IgATTGplusIgAEMA_IgAEMA - output$incremental_costs_IgATTGplusIgAEMA_IgAEMA
  #output$incremental_net_benefit_doubletest_IgAEMA <- 20000*output$incremental_effects_doubletest_IgAEMA - output$incremental_costs_doubletest_IgAEMA
  
  # Average incremental net benefit
  #output$average_inb_IgATTGplusIgAEMA_IgAEMA <- mean(output$incremental_net_benefit_IgATTGplusIgAEMA_IgAEMA)
  #output$average_inb_doubletest_IgAEMA <- mean(output$incremental_net_benefit_doubletest_IgAEMA)
  
  # Probability cost-effective
  # This is the proportion of samples for which the incremental net benefit is positive
  
  #output$probability_cost_effective_IgATTGplusIgAEMA_IgAEMA <- sum(output$incremental_net_benefit_IgATTGplusIgAEMA_IgAEMA > 0)/n_samples
  #output$probability_cost_effective_doubletest_IgAEMA <- sum(output$incremental_net_benefit_doubletest_IgAEMA > 0)/n_samples
  return(output)
}
  
  