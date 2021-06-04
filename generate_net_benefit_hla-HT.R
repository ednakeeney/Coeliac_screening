generate_net_benefit <- function(input_parameters) {
  
  starting_age <- ifelse(population == "adults", 50, 10) #based on mean age in under and over 18s in CPRD cost data
  n_cycles <- 90 - starting_age
  
  transition_matrices <- generate_transition_matrices(input_parameters)
  
  state_qalys <- generate_state_qalys(input_parameters)
  
  state_costs <- generate_state_costs(input_parameters)
  
 
  
  post_test_probability_IgAEMA <- input_parameters %>% select(X0.5.0.5 : X0.9999.0.9999)
  post_test_probability_IgATTGplusEMA <- input_parameters %>% select(X0.5.0.5.1 : X0.9999.0.9999.1)
  post_test_probability_IgATTG <- input_parameters %>% select(X0.5.0.5.2 : X0.9999.0.9999.2)
  pre_test_probability <- input_parameters %>% select(X0.5.0.5.3 : X0.9999.0.9999.3)
  fn_riskfactor <- input_parameters %>% select(X0.5.0.5.4 : X0.9999.0.9999.4)
  pre_test_probability_overall <- 0.01 #based on West 2014
  
 
  
  
   tp <- fn <- fp <- tn <- array(dim=c(n_samples, n_tests), dimnames = list(NULL, t_names))
  
  # Probabilities for IgAEMA
  for (i in 1:n_combinations) {
    tp[,1] <- 0
  tp[,i+1] <- ifelse(post_test_probability_IgAEMA[i] >= 0.9, (pre_test_probability[,i] * input_parameters$sens_IgAEMA), 
                   (pre_test_probability[,i] * input_parameters$sens_biopsy))
  fn[,1] <- pre_test_probability_overall  #in no screening all False Negatives
  fn[,i+1] <- pre_test_probability[,i] - tp[,i+1]  
  tn[,1] <- 1 - pre_test_probability_overall
  tn[,i+1] <- ifelse(post_test_probability_IgAEMA[i] >= 0.9, ((1 - pre_test_probability[,i]) * input_parameters$spec_IgAEMA),
                   ((1 - pre_test_probability[,i]) * input_parameters$spec_biopsy))
  fp[,1] <- 0
  fp[,i+1] <- (1 - pre_test_probability[,i]) - tn[,i+1]
  }
  
  # Probabilities for IgATTG + IgAEMA

    for (i in 1:n_combinations) {
  tp[,i+1+n_combinations] <- ifelse(post_test_probability_IgATTGplusEMA[i] >= 0.9, (pre_test_probability[,i] * input_parameters$sens_IgATTGplusEMA),
                   (pre_test_probability[,i] * input_parameters$sens_biopsy))
   fn[,i+1+n_combinations] <- pre_test_probability[,i] - tp[,i+1+n_combinations] 
   tn[,i+1+n_combinations] <- ifelse(post_test_probability_IgATTGplusEMA[i] >= 0.9, ((1 - pre_test_probability[,i]) * input_parameters$spec_IgATTGplusEMA),
                   ((1 - pre_test_probability[,i]) * input_parameters$spec_biopsy))
   fp[,i+1+n_combinations] <- (1 - pre_test_probability[,i]) - tn[,i+1+n_combinations]
    }
  
 # Probabilities for IgATTG

  for (i in 1:n_combinations) {
    tp[,i+1+n_combinations+n_combinations] <- ifelse(post_test_probability_IgATTG[i] >= 0.9, (pre_test_probability[,i] * input_parameters$sens_IgATTG),
                        (pre_test_probability[,i] * input_parameters$sens_biopsy))   
    fn[,i+1+n_combinations+n_combinations] <- pre_test_probability[,i] - tp[,i+1+n_combinations+n_combinations] 
    tn[,i+1+n_combinations+n_combinations] <- ifelse(post_test_probability_IgATTG[i] >= 0.9, ((1 - pre_test_probability[,i]) * input_parameters$spec_IgATTG),
                        ((1 - pre_test_probability[,i]) * input_parameters$spec_biopsy))
    fp[,i+1+n_combinations+n_combinations] <- (1 - pre_test_probability[,i]) - tn[,i+1+n_combinations+n_combinations]
  }
  
   
   pre_test_probability_HLA <- array(0, dim=c(n_samples, n_combinations*3), dimnames = list(NULL, t_names[2:109]))
   
   for (i in 1:n_combinations) {
     pre_test_probability_HLA[,i] <- tp[,i+1] /(tp[,i+1] + fp[,i+1])
     pre_test_probability_HLA[,i+n_combinations] <- tp[,i+1+n_combinations] /(tp[,i+1+n_combinations] + fp[,i+1+n_combinations])
     pre_test_probability_HLA[,i+n_combinations+n_combinations] <- tp[,i+1+n_combinations+n_combinations] /(tp[,i+1+n_combinations+n_combinations] + fp[,i+1+n_combinations+n_combinations])
   }
   
   pre_test_probability_HLA <- replace(pre_test_probability_HLA, pre_test_probability_HLA == 1, 0.99999)
   pre_test_odds_HLA <- array(0, dim=c(n_samples, n_combinations * 3), dimnames = list(NULL, t_names[2:109]))
   post_test_odds_HLA <- array(dim=c(n_samples, n_combinations * 3),dimnames=list(NULL, t_names[2:109]))
   post_test_probability_HLA <- array(dim=c(n_samples, n_combinations * 3),dimnames=list(NULL, t_names[2:109]))
  
   for (i in 1:108){
     pre_test_odds_HLA[,i] <- pre_test_probability_HLA[,i]/(1 - pre_test_probability_HLA[,i])
    post_test_odds_HLA[,i] <- pre_test_odds_HLA[,i] * input_parameters$LR_HLA
   post_test_probability_HLA[,i] <- post_test_odds_HLA[,i]/(1 + post_test_odds_HLA[,i])
   }
   
   # Probabilities for IgAEMA plus HLA
   for (i in 1:n_combinations) {
     tp[,i+1+(n_combinations*3)] <- ifelse(post_test_probability_HLA[,i] >= 0.9, (pre_test_probability_HLA[,i] * input_parameters$sens_HLA), 
                        (pre_test_probability_HLA[,i] * input_parameters$sens_biopsy))
     fn[,i+1+(n_combinations*3)] <- pre_test_probability_HLA[,i] - tp[,i+1+(n_combinations*3)]  
     tn[,i+1+(n_combinations*3)] <- ifelse(post_test_probability_HLA[,i] >= 0.9, ((1 - pre_test_probability_HLA[,i]) * input_parameters$spec_HLA),
                        ((1 - pre_test_probability_HLA[,i]) * input_parameters$spec_biopsy))
     fp[,i+1+(n_combinations*3)] <- (1 - pre_test_probability_HLA[,i]) - tn[,i+1+(n_combinations*3)]
   }
   
   # Probabilities for IgATTG plus IgaEMA plus HLA
   for (i in 1:n_combinations) {
     tp[,i+1+(n_combinations*4)] <- ifelse(post_test_probability_HLA[,i + n_combinations] >= 0.9, (pre_test_probability_HLA[,i + n_combinations] * input_parameters$sens_HLA), 
                                           (pre_test_probability_HLA[,i + n_combinations] * input_parameters$sens_biopsy))
     fn[,i+1+(n_combinations*4)] <- pre_test_probability_HLA[,i + n_combinations] - tp[,i+1+(n_combinations*4)]  
     tn[,i+1+(n_combinations*4)] <- ifelse(post_test_probability_HLA[,i + n_combinations] >= 0.9, ((1 - pre_test_probability_HLA[,i + n_combinations]) * input_parameters$spec_HLA),
                                           ((1 - pre_test_probability_HLA[,i + n_combinations]) * input_parameters$spec_biopsy))
     fp[,i+1+(n_combinations*4)] <- (1 - pre_test_probability_HLA[,i + n_combinations]) - tn[,i+1+(n_combinations*4)]
   }
   
   
   # Probabilities for IgATTG plus HLA
   for (i in 1:n_combinations) {
     tp[,i+1+(n_combinations*5)] <- ifelse(post_test_probability_HLA[,i + n_combinations + n_combinations] >= 0.9, (pre_test_probability_HLA[,i + n_combinations + n_combinations] * input_parameters$sens_HLA), 
                                           (pre_test_probability_HLA[,i + n_combinations + n_combinations] * input_parameters$sens_biopsy))
     fn[,i+1+(n_combinations*5)] <- pre_test_probability_HLA[,i + n_combinations + n_combinations] - tp[,i+1+(n_combinations*5)]  
     tn[,i+1+(n_combinations*5)] <- ifelse(post_test_probability_HLA[,i + n_combinations + n_combinations] >= 0.9, ((1 - pre_test_probability_HLA[,i + n_combinations + n_combinations]) * input_parameters$spec_HLA),
                                           ((1 - pre_test_probability_HLA[,i + n_combinations + n_combinations]) * input_parameters$spec_biopsy))
     fp[,i+1+(n_combinations*5)] <- (1 - pre_test_probability_HLA[,i + n_combinations + n_combinations]) - tn[,i+1+(n_combinations*5)]
   }
 
   write.csv(data.frame(colMeans(tp), colMeans(fn), colMeans(tn), colMeans(fp)), "test.csv")
   
     
  fp_costs <- array(dim = c(n_samples, n_tests),
                    dimnames = list(NULL, t_names))
  
  fp_costs[,] <- (fp[,] * input_parameters$cost_gfp) 
  
  diagnosis_costs <- array(dim = c(n_samples, n_tests),
                           dimnames = list(NULL, t_names))
  
  diagnosis_costs[,] <- (tp[,] + fp[,]) * input_parameters$cost_diagnosis
  
  
  test_costs <- array(dim=c(n_samples, n_tests), dimnames=list(NULL, t_names))
  for (i in 1:n_combinations) {
    test_costs[, 1] <- 0
    test_costs[, i+1] <-  ifelse(post_test_probability_IgAEMA[,i] >= 0.9, input_parameters$test_cost_IgAEMA, 
                                      input_parameters$test_cost_IgAEMA + input_parameters$cost_biopsy)
    test_costs[, n_combinations+i+1] <-  ifelse(post_test_probability_IgATTGplusEMA[,i] >= 0.9, (input_parameters$test_cost_IgATTG + input_parameters$test_cost_IgAEMA), 
                                                     (input_parameters$test_cost_IgATTG + input_parameters$test_cost_IgAEMA + input_parameters$cost_biopsy))
    test_costs[, n_combinations+n_combinations+i+1] <-  ifelse(post_test_probability_IgATTG[,i] >= 0.9, input_parameters$test_cost_IgATTG, 
                                                                    input_parameters$test_cost_IgATTG + input_parameters$cost_biopsy)
    test_costs[, n_combinations+n_combinations+ n_combinations+i+1] <- 
      ifelse(post_test_probability_IgAEMA[,i] >= 0.9,
        input_parameters$test_cost_IgAEMA,
      ifelse(post_test_probability_IgAEMA[,i] < 0.9 & post_test_probability_HLA[,i] >= 0.9,
          (input_parameters$test_cost_IgAEMA + input_parameters$test_cost_HLA),
        (input_parameters$test_cost_IgAEMA + input_parameters$test_cost_HLA + input_parameters$cost_biopsy)))
    
    
    test_costs[, n_combinations + n_combinations+ n_combinations + n_combinations + i + 1] <-  
      ifelse(post_test_probability_IgATTGplusEMA[,i] >= 0.9,
             input_parameters$test_cost_IgAEMA + input_parameters$test_cost_IgATTG,
      ifelse(post_test_probability_IgATTGplusEMA[,i] < 0.9 & post_test_probability_HLA[,i+n_combinations] >= 0.9, 
             input_parameters$test_cost_IgAEMA + input_parameters$test_cost_IgATTG + input_parameters$test_cost_HLA, 
      input_parameters$test_cost_IgAEMA + input_parameters$test_cost_IgATTG + input_parameters$test_cost_HLA + input_parameters$cost_biopsy))

    test_costs[, n_combinations + n_combinations+ n_combinations + n_combinations+ n_combinations+i+1] <-  
      ifelse(post_test_probability_IgATTG[,i] >= 0.9,
             input_parameters$test_cost_IgATTG,
      ifelse(post_test_probability_IgATTG[,i] < 0.9 & post_test_probability_HLA[,i+n_combinations+n_combinations] >= 0.9, 
             input_parameters$test_cost_IgATTG + input_parameters$test_cost_HLA, 
   input_parameters$test_cost_IgATTG + input_parameters$test_cost_HLA + input_parameters$cost_biopsy))
  
  }
  
  test_costs_applied <- array(dim = c(n_samples, n_tests),
                              dimnames = list(NULL, t_names))
  
  tp_riskfactor_table <- array(dim=c(n_samples, n_tests), dimnames = list(NULL, t_names))
  #ncol(tp_riskfactor_table)
  for (i in 1:n_combinations) {
    tp_riskfactor_table[,1] <- 1
    tp_riskfactor_table[,i+1] <- tp_riskfactor[,i]
    tp_riskfactor_table[,i+1+n_combinations] <- tp_riskfactor[,i]
    tp_riskfactor_table[,i+1+n_combinations+n_combinations] <- tp_riskfactor[,i]
    tp_riskfactor_table[,i+1+n_combinations*3] <- tp_riskfactor[,i]
    tp_riskfactor_table[,i+1+n_combinations*4] <- tp_riskfactor[,i]
    tp_riskfactor_table[,i+1+n_combinations*5] <- tp_riskfactor[,i]
  }
  
  fp_riskfactor_table <- array(dim=c(n_samples, n_tests), dimnames = list(NULL, t_names))
  #ncol(fp_riskfactor_table)
  for (i in 1:n_combinations) {
    fp_riskfactor_table[,1] <- 1
    fp_riskfactor_table[,i+1] <- fp_riskfactor[,i]
    fp_riskfactor_table[,i+1+n_combinations] <- fp_riskfactor[,i]
    fp_riskfactor_table[,i+1+n_combinations+n_combinations] <- fp_riskfactor[,i]
    fp_riskfactor_table[,i+1+n_combinations*3] <- fp_riskfactor[,i]
    fp_riskfactor_table[,i+1+n_combinations*4] <- fp_riskfactor[,i]
    fp_riskfactor_table[,i+1+n_combinations*5] <- fp_riskfactor[,i]
  }
  test_costs_applied <- test_costs * (tp[,] + fp[,] + fn[,] + tn[,] + tp_riskfactor_table[,] + fp_riskfactor_table[,])
  
  
 
fn_riskfactor_table <- array(dim=c(n_samples, n_tests), dimnames = list(NULL, t_names))
  #ncol(fn_riskfactor_table)
  for (i in 1:n_combinations) {
  fn_riskfactor_table[,1] <- 1
  fn_riskfactor_table[,i+1] <- fn_riskfactor[,i]
  fn_riskfactor_table[,i+1+n_combinations] <- fn_riskfactor[,i]
  fn_riskfactor_table[,i+1+n_combinations+n_combinations] <- fn_riskfactor[,i]
  fn_riskfactor_table[,i+1+n_combinations*3] <- fn_riskfactor[,i]
  fn_riskfactor_table[,i+1+n_combinations*4] <- fn_riskfactor[,i]
  fn_riskfactor_table[,i+1+n_combinations*5] <- fn_riskfactor[,i]
  }
fn_riskfactor_table <- fn_riskfactor_table * 1/(fp+tp)
   fn_all <- fn + fn_riskfactor_table
   

   #fptntocoeliac_ratio <- (tn+fp)/(tp+fn_all)
  
  #scaling up true positives and false negatives 
  tp <- tp/(tp + fn_all)
  fn <- 1 - tp
  
  # Build an array to store the cohort vector at each cycle
  # Each cohort vector has n_states elements: probability of being in each state,
  # There is one cohort vector for each test, for each PSA sample, for each cycle.
  cohort_vectors<-array(dim=c(n_tests,n_samples,n_cycles,n_states),  
                        dimnames=list(t_names,NULL,NULL, state_names))
  
  
  for (i_test in c(1:n_tests)) { 
    cohort_vectors[i_test, , 1, "CD GFD no complications"] <- tp[, i_test] * input_parameters$probability_nocomplications 
    cohort_vectors[i_test, , 1, "CD GFD subfertility"] <- tp[, i_test] * input_parameters$probability_subfertility
    cohort_vectors[i_test, , 1, "CD GFD osteoporosis"] <- tp[, i_test] * input_parameters$probability_osteoporosis 
    cohort_vectors[i_test, , 1, "CD GFD NHL"] <- tp[, i_test] * input_parameters$probability_NHL
    
    cohort_vectors[i_test, , 1, "Undiagnosed CD no complications"] <- fn[, i_test] * input_parameters$probability_nocomplications 
    cohort_vectors[i_test, , 1, "Undiagnosed CD subfertility"] <- fn[, i_test] * input_parameters$probability_subfertility 
    cohort_vectors[i_test, , 1, "Undiagnosed CD osteoporosis"] <- fn[, i_test] * input_parameters$probability_osteoporosis 
    cohort_vectors[i_test, , 1, "Undiagnosed CD NHL"] <- fn[, i_test] * input_parameters$probability_NHL 
    
    cohort_vectors[i_test, , 1, "Death"] <- 0
  }
  
  cohort_vectors[, , , ] [is.na(cohort_vectors[, , , ] )] <- 0
  rowSums (cohort_vectors[, 2, , ], na.rm = FALSE, dims = 1)
  
  # Build an array to store the costs and QALYs accrued per cycle
  # One for each test, for each PSA sample, for each cycle
  # These will be filled in below in the main model code
  # Then discounted and summed to contribute to total costs and total QALYs
  cycle_costs <- array(dim = c(n_tests, n_samples, n_cycles),
                       dimnames = list(t_names, NULL, NULL))
  cycle_qalys <- array(dim = c(n_tests, n_samples, n_cycles),
                       dimnames = list(t_names, NULL, NULL))
  
  # Build arrays to store the total costs and total QALYs
  # There is one for each test and each PSA sample
  # These are filled in below using cycle_costs, 
  # test_costs, and cycle_qalys
  total_costs <- array(dim = c(n_tests, n_samples),
                       dimnames = list(t_names, NULL))
  total_qalys <- array(dim = c(n_tests, n_samples),
                       dimnames = list(t_names, NULL))
  
  
  disc_vec <- (1/1.035) ^ rep(c(0:(n_cycles/2-1)), each = 2)
  
  # Main model code
  # Loop over the test options
  
  for (i_test in 1:n_tests)
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
        cohort_vectors[i_test, i_sample, i_cycle, ] <-
          cohort_vectors[i_test, i_sample, i_cycle-1, ] %*%
          transition_matrices_sample[i_cycle, , ]
      }
      
      cohort_vectors_tr_sample <- cohort_vectors[i_test, i_sample, , ]
      
      for (i_cycle in 1:n_cycles)
      {
        # This was inside the loop, but overwrote cycle_costs[i_test, i_sample,] for
        # each cycle so can't be right
        # Now use the cohort vectors to calculate the 
        # total costs for each cycle
        #cycle_costs[i_test, i_sample, ] <-
        #  cohort_vectors_tr_sample %*% state_costs[i_sample, i_cycle, ]
        
        # I think this needs to be cycle specific
        # The below is from the intro to Markov modelling code
        cycle_costs[i_test, i_sample, i_cycle] <-
          cohort_vectors[i_test, i_sample, i_cycle, ] %*% state_costs[i_sample, i_cycle, ]
        # And total QALYs for each cycle
        cycle_qalys[i_test, i_sample, i_cycle] <-
          cohort_vectors[i_test, i_sample, i_cycle, ]  %*% state_qalys[i_sample, i_cycle, ]
        
      }
      
      # Combine the cycle_costs and other costs to get total costs
      # Apply the discount factor 
      # (1 in first year, 1.035 in second, 1.035^2 in third, and so on)
      # Each year acounts for two cycles so need to repeat the discount values
      total_costs[i_test, i_sample] <- test_costs_applied[i_sample, i_test] + fp_costs[i_sample, i_test] + diagnosis_costs[i_sample, i_test] + cycle_costs[i_test, i_sample, ] %*% disc_vec
      
      # Combine the cycle_qalys to get total qalys
      # Apply the discount factor 
      # (1 in first year, 1.035 in second, 1.035^2 in third, and so on)
      # Each year acounts for two cycles so need to repeat the discount values
      total_qalys[i_test, i_sample] <- cycle_qalys[i_test, i_sample, ] %*%
        disc_vec
    }
    
  }
  

  #############################################################################
  ## Analysis of results ######################################################
  #############################################################################
  output <- list()
  sum_IgAEMA <- rep(0, n_combinations)
  for (i in 1:n_combinations){
    sum_IgAEMA[i] <- sum(post_test_probability_IgAEMA[,i] <= 0.9)/n_samples  
  }
  names(sum_IgAEMA) <- combinations_names
  
  sum_IgATTGplusEMA <- rep(0, n_combinations)
  for (i in 1:n_combinations){
    sum_IgATTGplusEMA[i] <- sum(post_test_probability_IgATTGplusEMA[,i] <= 0.9)/n_samples  
  }
  names(sum_IgATTGplusEMA) <- combinations_names
  
  sum_IgATTG <- rep(0, n_combinations)
  for (i in 1:n_combinations){
    sum_IgATTG[i] <- sum(post_test_probability_IgATTG[,i] <= 0.9)/n_samples  
  }
  names(sum_IgATTG) <- combinations_names
  
  sum_IgAEMAplusHLA <- rep(0, n_combinations)
  for (i in 1:n_combinations){
    sum_IgAEMAplusHLA[i] <- sum(post_test_probability_HLA[,i] <= 0.9)/n_samples  
  }
  names(sum_IgAEMAplusHLA) <- t_names[110:145]
  
  sum_IgATTGplusEMAplusHLA <- rep(0, n_combinations)
  for (i in 1:n_combinations){
    sum_IgATTGplusEMAplusHLA[i] <- sum(post_test_probability_HLA[,i+n_combinations] <= 0.9)/n_samples  
  }
  names(sum_IgATTGplusEMAplusHLA) <- t_names[146:181]
  
  sum_IgATTGplusHLA <- rep(0, n_combinations)
  for (i in 1:n_combinations){
    sum_IgATTGplusHLA[i] <- sum(post_test_probability_HLA[,i+n_combinations+n_combinations] <= 0.9)/n_samples  
  }
  names(sum_IgATTGplusHLA) <- t_names[182:217]
  
  output$percentage_biopsy_IgAEMA <- sum_IgAEMA
  output$percentage_biopsy_IgATTGplusEMA <- sum_IgATTGplusEMA
  output$percentage_biopsy_IgATTG <- sum_IgATTG
  output$percentage_biopsy_IgAEMAplusHLA <- sum_IgAEMAplusHLA
  output$percentage_biopsy_IgATTGplusEMAplusHLA <- sum_IgATTGplusEMAplusHLA
  output$percentage_biopsy_IgATTGplusHLA <- sum_IgATTGplusHLA
  
  output$total_costs <- total_costs
  output$total_qalys <- total_qalys
  # Average costs
  output$average_costs <- rowMeans(total_costs)
  # Average effects (in QALY units)
  output$average_effects <- rowMeans(total_qalys)
  
 
  output$incremental_costs <-  output$average_costs - output$average_costs["No screening"]
  output$incremental_effects <-  output$average_effects -  output$average_effects["No screening"]
  
  output$ICER <- output$incremental_costs/output$incremental_effects

  # Incremental net benefit at the ?20,000 willingness-to-pay
  
  output$incremental_net_benefit <- 20000*output$incremental_effects - output$incremental_costs
 write.csv(output$incremental_effects, "inb.csv")
  
  output$test_costs <- colMeans(t(test_costs)) #costs of test and biopsies
  output$fp_costs <- colMeans(fp_costs)
  output$diagnosis_costs <- colMeans(diagnosis_costs)
 
 # Average incremental net benefit
  #output$average_inb_IgATTGplusIgAEMA_IgAEMA <- mean(output$incremental_net_benefit_IgATTGplusIgAEMA_IgAEMA)
  #output$average_inb_doubletest_IgAEMA <- mean(output$incremental_net_benefit_doubletest_IgAEMA)
  
  # Probability cost-effective
  # This is the proportion of samples for which the incremental net benefit is positive
  
  #output$probability_cost_effective_IgATTGplusIgAEMA_IgAEMA <- sum(output$incremental_net_benefit_IgATTGplusIgAEMA_IgAEMA > 0)/n_samples
  #output$probability_cost_effective_doubletest_IgAEMA <- sum(output$incremental_net_benefit_doubletest_IgAEMA > 0)/n_samples
  return(output)
}
  
  