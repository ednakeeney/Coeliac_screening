generate_net_benefit <- function(transition_matrices, state_costs, state_qalys) {
  
  
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
      
      # Combine the cycle_costs and treatment_costs to get total costs
      # Apply the discount factor 
      # (1 in first year, 1.035 in second, 1.035^2 in third, and so on)
      # Each year acounts for two cycles so need to repeat the discount values
      total_costs[i_treatment, i_sample] <- treatment_costs[i_treatment, i_sample] +
        cycle_costs[i_treatment, i_sample, ] %*%
        disc_vec
      
      # Combine the cycle_qalys to get total qalys
      # Apply the discount factor 
      # (1 in first year, 1.035 in second, 1.035^2 in third, and so on)
      # Each year acounts for two cycles so need to repeat the discount values
      total_qalys[i_treatment, i_sample] <- cycle_qalys[i_treatment, i_sample, ] %*%
        disc_vec
    }
    
  }
  
  # HT: Weird error that cohort_vectors for final cycle don't sum to 1
  # e.g. sum(cohort_vectors[1,20,50,]) is less than one but sum(cohort_vectors[1,20,49,]) is one
  # It's caused by the transition_matrices for 50th cycle not summing to one
  # In particule the row for "CD no GFD no complications"
  #############################################################################
  ## Analysis of results ######################################################
  #############################################################################
  output <- list()
  # Average costs
  output$average_costs <- rowMeans(total_costs)
  # Average effects (in QALY units)
  output$average_effects <- rowMeans(total_qalys)
  
  
  output$incremental_costs_testbiopsy_test <- total_costs["Test + biopsy", ] - total_costs["Test", ]
  output$incremental_effects_testbiopsy_test <- total_qalys["Test + biopsy", ] - total_qalys["Test", ]
  
  output$incremental_costs_doubletest_test <- total_costs["Double test", ] - total_costs["Test", ]
  output$incremental_effects_doubletest_test <- total_qalys["Double test", ] - total_qalys["Test", ]
  
  
  output$ICER_testbiopsy_test <- mean(output$incremental_costs_testbiopsy_test)/mean(output$incremental_effects_testbiopsy_test)
  output$ICER_doubletest_test <- mean(output$incremental_costs_doubletest_test)/mean(output$incremental_effects_doubletest_test)
  
  # Incremental net benefit at the ?20,000 willingness-to-pay
  
  output$incremental_net_benefit_testbiopsy_test <- 20000*output$incremental_effects_testbiopsy_test - output$incremental_costs_testbiopsy_test
  output$incremental_net_benefit_doubletest_test <- 20000*output$incremental_effects_doubletest_test - output$incremental_costs_doubletest_test
  
  # Average incremental net benefit
  output$average_inb_testbiopsy_test <- mean(output$incremental_net_benefit_testbiopsy_test)
  output$average_inb_doubletest_test <- mean(output$incremental_net_benefit_doubletest_test)
  
  # Probability cost-effective
  # This is the proportion of samples for which the incremental net benefit is positive
  
  output$probability_cost_effective_testbiopsy_test <- sum(output$incremental_net_benefit_testbiopsy_test > 0)/n_samples
  output$probability_cost_effective_doubletest_test <- sum(output$incremental_net_benefit_doubletest_test > 0)/n_samples
  return(output)
}
  
  