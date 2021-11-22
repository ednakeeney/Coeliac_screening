generate_net_benefit <- function(input_parameters, 
                                 strategies_of_interest,
                                 transition_matrices,
                                 combinations,
                                 population = NULL) {

  # Derived descriptions of inputs
  n_combinations <- dim(combinations)[1]
  n_samples <- dim(input_parameters)[1]
    
  starting_age <- ifelse(population == "men" | population == "women", 18, 10) #based on mean age in under and over 18s in CPRD cost data
  n_cycles <- 90 - starting_age
  
  state_qalys <- generate_state_qalys(input_parameters)
  
  state_costs <- generate_state_costs(input_parameters)
  
 

  # Define the strategies based on sampled test sensitivity and specificity
  pre_test_probability_overall <- input_parameters$pre_test_probability_overall #based on West 2014
  
  
  # Sensitivity and specificity parameters
  spec_IgAEMA_adults <- input_parameters$spec_IgAEMA_adults
  sens_IgAEMA_adults <- input_parameters$sens_IgAEMA_adults
  sens_IgATTGplusEMA_adults <- input_parameters$sens_IgATTGplusEMA_adults
  spec_IgATTGplusEMA_adults <- input_parameters$spec_IgATTGplusEMA_adults
  sens_IgATTG_adults <- input_parameters$sens_IgATTG_adults
  spec_IgATTG_adults <- input_parameters$spec_IgATTG_adults
  
  spec_IgAEMA_children <- input_parameters$spec_IgAEMA_children
  sens_IgAEMA_children <- input_parameters$sens_IgAEMA_children
  sens_IgATTGplusEMA_children <- input_parameters$sens_IgATTGplusEMA_children
  spec_IgATTGplusEMA_children <- input_parameters$spec_IgATTGplusEMA_children
  sens_IgATTG_children <- input_parameters$sens_IgATTG_children
  spec_IgATTG_children <- input_parameters$spec_IgATTG_children
  
  # HLA not split by age
  sens_HLA <- input_parameters$sens_HLA
  spec_HLA <- input_parameters$spec_HLA
  
  
  ## Define the true positives etc and post test probabilities of the serological testing strategies
  
  tp_riskfactor <- array(dim=c(n_samples, n_combinations), dimnames = list(NULL, paste(combinations_names, "tp_riskfactor")))
  fn_riskfactor <- array(dim=c(n_samples, n_combinations), dimnames = list(NULL, paste(combinations_names, "fn_riskfactor")))
  fp_riskfactor <- array(dim=c(n_samples, n_combinations), dimnames = list(NULL, paste(combinations_names, "fp_riskfactor")))
  tn_riskfactor <- array(dim=c(n_samples, n_combinations), dimnames = list(NULL, paste(combinations_names, "tn_riskfactor")))
  
  for (i in 1:n_combinations) {
    tp_riskfactor[,i] <- pre_test_probability_overall * combinations$sens_riskfactor[i]
    fn_riskfactor[,i] <- pre_test_probability_overall - tp_riskfactor[,i]  
    tn_riskfactor[,i] <- (1 - pre_test_probability_overall) * combinations$spec_riskfactor[i]
    fp_riskfactor[,i] <- (1 - pre_test_probability_overall) - tn_riskfactor[,i]
  }
  
  write.csv(data.frame(tp_riskfactor[1,], fn_riskfactor[1,], tn_riskfactor[1,], fp_riskfactor[1,]), "risk_factor.csv")
  pre_test_probability <- tp_riskfactor/(tp_riskfactor+fp_riskfactor)
  colnames(pre_test_probability) <- paste(combinations_names, "pre_test_probability")
  write.csv(colMeans(pre_test_probability), "pretestprob.csv")
  
  pre_test_odds <- array(0, dim=c(n_samples, n_combinations), dimnames = list(NULL, combinations_names))
  
  for (i in 1:n_combinations){
    pre_test_odds[,i] <- pre_test_probability[,i]/(1 - pre_test_probability[,i])
  }
  
  ##################################################################################################
  
  # IGAEmA
  LR_IgAEMA_adults <- sens_IgAEMA_adults/ (1 - spec_IgAEMA_adults)
  post_test_odds_IgAEMA_adults <- array(dim=c(n_samples, n_combinations),dimnames=list(NULL, combinations_names))
  post_test_probability_IgAEMA_adults <- array(dim=c(n_samples, n_combinations),dimnames=list(NULL, combinations_names))
  
  for (i in 1:n_combinations){
    post_test_odds_IgAEMA_adults[,i] <- pre_test_odds[,i] * LR_IgAEMA_adults
    post_test_probability_IgAEMA_adults[,i] <- post_test_odds_IgAEMA_adults[,i]/(1 + post_test_odds_IgAEMA_adults[,i])
  }
  
  LR_IgAEMA_children <- sens_IgAEMA_children/ spec_IgAEMA_children
  
  post_test_odds_IgAEMA_children <- array(dim=c(n_samples, n_combinations),dimnames=list(NULL, combinations_names))
  post_test_probability_IgAEMA_children <- array(dim=c(n_samples, n_combinations),dimnames=list(NULL, combinations_names))
  
  for (i in 1:n_combinations){
    post_test_odds_IgAEMA_children[,i] <- pre_test_odds[,i] * LR_IgAEMA_children
    post_test_probability_IgAEMA_children[,i] <- post_test_odds_IgAEMA_children[,i]/(1 + post_test_odds_IgAEMA_children[,i])
  }
  

  post_test_probability_IgAEMA <- array(dim=c(n_samples, n_combinations),dimnames=list(NULL, paste(combinations_names, "IgAEMA")))
  sens_IgAEMA <- rep(0, times = n_samples)
  spec_IgAEMA <- rep(0, times = n_samples)
  
  # Choose correct accuracies for this population
  for(i_sample in 1:n_samples){
    sens_IgAEMA[i_sample] <- ifelse(population == "men" | population == "women", sens_IgAEMA_adults[i_sample], sens_IgAEMA_children[i_sample])
    spec_IgAEMA[i_sample] <- ifelse(population == "men" | population == "women", spec_IgAEMA_adults[i_sample], spec_IgAEMA_children[i_sample])
  }
  
  for(i_sample in 1:n_samples){
    for (i in 1:n_combinations){
      post_test_probability_IgAEMA[i_sample,i] <- ifelse(population == "men" | population == "women", post_test_probability_IgAEMA_adults[i_sample,i], post_test_probability_IgAEMA_children[i_sample,i])
    }}
  
  
  
  ##################################################################################################
  # IgATTGplusEMA
  
  LR_IgATTGplusEMA_adults <- sens_IgATTGplusEMA_adults/ (1 - spec_IgATTGplusEMA_adults)
  
  post_test_odds_IgATTGplusEMA_adults <- array(dim=c(n_samples, n_combinations),dimnames=list(NULL, combinations_names))
  post_test_probability_IgATTGplusEMA_adults <- array(dim=c(n_samples, n_combinations),dimnames=list(NULL, paste(combinations_names, "IgATTGplusEMA")))
  
  for (i in 1:n_combinations){
    post_test_odds_IgATTGplusEMA_adults[,i] <- pre_test_odds[,i] * LR_IgATTGplusEMA_adults[1:n_samples]
    post_test_probability_IgATTGplusEMA_adults[,i] <- post_test_odds_IgATTGplusEMA_adults[,i]/(1 + post_test_odds_IgATTGplusEMA_adults[,i])
  }
  
  LR_IgATTGplusEMA_children <- sens_IgATTGplusEMA_children/ (1 - spec_IgATTGplusEMA_children)
  
  post_test_odds_IgATTGplusEMA_children <- array(dim=c(n_samples, n_combinations),dimnames=list(NULL, combinations_names))
  post_test_probability_IgATTGplusEMA_children <- array(dim=c(n_samples, n_combinations),dimnames=list(NULL, paste(combinations_names, "IgATTGplusEMA")))
  
  for (i in 1:n_combinations){
    post_test_odds_IgATTGplusEMA_children[,i] <- pre_test_odds[,i] * LR_IgATTGplusEMA_children[1:n_samples]
    post_test_probability_IgATTGplusEMA_children[,i] <- post_test_odds_IgATTGplusEMA_children[,i]/(1 + post_test_odds_IgATTGplusEMA_children[,i])
  }
  
  post_test_probability_IgATTGplusEMA <- array(dim=c(n_samples, n_combinations),dimnames=list(NULL, paste(combinations_names, "IgATTGplusEMA")))
  sens_IgATTGplusEMA <- rep(0, times = n_samples)
  spec_IgATTGplusEMA <- rep(0, times = n_samples)
  for(i_sample in 1:n_samples){
    for (i in 1:n_combinations){
      post_test_probability_IgATTGplusEMA[i_sample,i] <- ifelse(population == "men" | population == "women", post_test_probability_IgATTGplusEMA_adults[i_sample,i], post_test_probability_IgATTGplusEMA_children[i_sample,i])
    }}
  
  for(i_sample in 1:n_samples){
    sens_IgATTGplusEMA[i_sample] <- ifelse(population == "men" | population == "women", sens_IgATTGplusEMA_adults[i_sample], sens_IgATTGplusEMA_children[i_sample])
    spec_IgATTGplusEMA[i_sample] <- ifelse(population == "men" | population == "women", spec_IgATTGplusEMA_adults[i_sample], spec_IgATTGplusEMA_children[i_sample])
  }
  
  ##################################################################################################
  # IgATTG
  
  LR_IgATTG_adults <- sens_IgATTG_adults/ (1 - spec_IgATTG_adults)
  
  post_test_odds_IgATTG_adults <- array(dim=c(n_samples, n_combinations),dimnames=list(NULL, combinations_names))
  post_test_probability_IgATTG_adults <- array(dim=c(n_samples, n_combinations),dimnames=list(NULL, combinations_names))
  
  for (i in 1:n_combinations){
    post_test_odds_IgATTG_adults[,i] <- pre_test_odds[,i] * LR_IgATTG_adults
    post_test_probability_IgATTG_adults[,i] <- post_test_odds_IgATTG_adults[,i]/(1 + post_test_odds_IgATTG_adults[,i])
  }
  
  LR_IgATTG_children <- sens_IgATTG_children/ spec_IgATTG_children
  
  post_test_odds_IgATTG_children <- array(dim=c(n_samples, n_combinations),dimnames=list(NULL, combinations_names))
  post_test_probability_IgATTG_children <- array(dim=c(n_samples, n_combinations),dimnames=list(NULL, combinations_names))
  
  for (i in 1:n_combinations){
    post_test_odds_IgATTG_children[,i] <- pre_test_odds[,i] * LR_IgATTG_children
    post_test_probability_IgATTG_children[,i] <- post_test_odds_IgATTG_children[,i]/(1 + post_test_odds_IgATTG_children[,i])
  }
  
  
  post_test_probability_IgATTG <- array(dim=c(n_samples, n_combinations),dimnames=list(NULL, paste(combinations_names, "IgATTG")))
  sens_IgATTG <- rep(0, times = n_samples)
  spec_IgATTG <- rep(0, times = n_samples)
  for(i_sample in 1:n_samples){
    for (i in 1:n_combinations){
      post_test_probability_IgATTG[i_sample,i] <- ifelse(population == "men" | population == "women", post_test_probability_IgATTG_adults[i_sample,i], post_test_probability_IgATTG_children[i_sample,i])
    }}
  
  for(i_sample in 1:n_samples){
    sens_IgATTG[i_sample] <- ifelse(population == "men" | population == "women", sens_IgATTG_adults[i_sample], sens_IgATTG_children[i_sample])
    spec_IgATTG[i_sample] <- ifelse(population == "men" | population == "women", spec_IgATTG_adults[i_sample], spec_IgATTG_children[i_sample])
  }
  
  # Likelihood ratios
  LR_IgATTG <- sens_IgATTG/spec_IgATTG
  LR_IgAEMA <- sens_IgAEMA/spec_IgAEMA
  LR_IgATTGplusEMA <- sens_IgATTGplusEMA/spec_IgATTGplusEMA
  
  # HLA
  
  LR_HLA <-  sens_HLA / spec_HLA
  post_test_odds_HLAalone  <- post_test_probability_HLAalone <- array(dim=c(n_samples, n_combinations),dimnames=list(NULL, paste(combinations_names, "IgATTG")))
  for (i in 1:n_combinations){
    post_test_odds_HLAalone[,i] <- pre_test_odds[,i] * LR_HLA
    post_test_probability_HLAalone[,i] <- post_test_odds_HLAalone[,i]/(1 + post_test_odds_HLAalone[,i])
  }
  
   tp <- fn <- fp <- tn <- array(dim=c(n_samples, n_tests), dimnames = list(NULL, t_names))
  

  for (i in 1:n_combinations) {
   
     #Probabilities for no screening
    tp[,1] <- 0 
    fn[,1] <- pre_test_probability_overall  #in no screening all False Negatives
    tn[,1] <- 1 - pre_test_probability_overall
    fp[,1] <- 0
    
    # Probabilities for IgAEMA
  tp[,i+1] <- ifelse(post_test_probability_IgAEMA[,i] >= 0.9, (pre_test_probability[,i] * sens_IgAEMA), 
                   (pre_test_probability[,i] * input_parameters$sens_biopsy))
  fn[,i+1] <- pre_test_probability[,i] - tp[,i+1]  
  tn[,i+1] <- ifelse(post_test_probability_IgAEMA[,i] >= 0.9, ((1 - pre_test_probability[,i]) * spec_IgAEMA),
                   ((1 - pre_test_probability[,i]) * input_parameters$spec_biopsy))
  fp[,i+1] <- (1 - pre_test_probability[,i]) - tn[,i+1]
  
  # Probabilities for IgATTG + IgAEMA
    tp[,i+1+n_combinations] <- ifelse(post_test_probability_IgATTGplusEMA[,i] >= 0.9, (pre_test_probability[,i] * sens_IgATTGplusEMA),
                                      (pre_test_probability[,i] * input_parameters$sens_biopsy))
    fn[,i+1+n_combinations] <- pre_test_probability[,i] - tp[,i+1+n_combinations] 
    tn[,i+1+n_combinations] <- ifelse(post_test_probability_IgATTGplusEMA[,i] >= 0.9, ((1 - pre_test_probability[,i]) * spec_IgATTGplusEMA),
                                      ((1 - pre_test_probability[,i]) * input_parameters$spec_biopsy))
    fp[,i+1+n_combinations] <- (1 - pre_test_probability[,i]) - tn[,i+1+n_combinations]
    
    # Probabilities for IgATTG
      tp[,i+1+n_combinations+n_combinations] <- ifelse(post_test_probability_IgATTG[,i] >= 0.9, (pre_test_probability[,i] * sens_IgATTG),
                                                       (pre_test_probability[,i] * input_parameters$sens_biopsy))   
      fn[,i+1+n_combinations+n_combinations] <- pre_test_probability[,i] - tp[,i+1+n_combinations+n_combinations] 
      tn[,i+1+n_combinations+n_combinations] <- ifelse(post_test_probability_IgATTG[,i] >= 0.9, ((1 - pre_test_probability[,i]) * spec_IgATTG),
                                                       ((1 - pre_test_probability[,i]) * input_parameters$spec_biopsy))
      fp[,i+1+n_combinations+n_combinations] <- (1 - pre_test_probability[,i]) - tn[,i+1+n_combinations+n_combinations]
  
    
  }

   tp_nobiopsy <- fn_nobiopsy <- fp_nobiopsy <- tn_nobiopsy <- array(dim=c(n_samples, n_tests), dimnames = list(NULL, t_names))
   
   
   for (i in 1:n_combinations) {
     
     #Probabilities for no screening
     tp_nobiopsy[,1] <- 0
     fn_nobiopsy[,1] <- pre_test_probability_overall  #in no screening all False Negatives
     tn_nobiopsy[,1] <- 1 - pre_test_probability_overall
     fp_nobiopsy[,1] <- 0
     
     
     # Probabilities for IgAEMA without biopsy (to use for HLA strategies)
     tp_nobiopsy[,i+1] <- pre_test_probability[,i] * sens_IgAEMA
     fn_nobiopsy[,i+1] <- pre_test_probability[,i] - tp_nobiopsy[,i+1]  
     tn_nobiopsy[,i+1] <-  ((1 - pre_test_probability[,i]) * spec_IgAEMA)
     fp_nobiopsy[,i+1] <- (1 - pre_test_probability[,i]) - tn_nobiopsy[,i+1]
     
     # Probabilities for IgATTG + IgAEMA without biopsy
     tp_nobiopsy[,i+1+n_combinations] <- pre_test_probability[,i] * sens_IgATTGplusEMA
     fn_nobiopsy[,i+1+n_combinations] <- pre_test_probability[,i] - tp_nobiopsy[,i+1+n_combinations] 
     tn_nobiopsy[,i+1+n_combinations] <- (1 - pre_test_probability[,i]) * spec_IgATTGplusEMA
     fp_nobiopsy[,i+1+n_combinations] <- (1 - pre_test_probability[,i]) - tn_nobiopsy[,i+1+n_combinations]
     
     # Probabilities for IgATTG without biopsy
     tp_nobiopsy[,i+1+n_combinations+n_combinations] <- pre_test_probability[,i] * sens_IgATTG
     fn_nobiopsy[,i+1+n_combinations+n_combinations] <- pre_test_probability[,i] - tp_nobiopsy[,i+1+n_combinations+n_combinations] 
     tn_nobiopsy[,i+1+n_combinations+n_combinations] <- (1 - pre_test_probability[,i]) * spec_IgATTG
     fp_nobiopsy[,i+1+n_combinations+n_combinations] <- (1 - pre_test_probability[,i]) - tn_nobiopsy[,i+1+n_combinations+n_combinations]
     
     # Probabilities for HLA without biopsy
     # Whether paired with TTG, TTG+EMA or EMA they are the same
     for(j in 6:8) {
       tp_nobiopsy[,i+1+j*n_combinations] <- pre_test_probability[,i] * sens_HLA
       fn_nobiopsy[,i+1+j*n_combinations] <- pre_test_probability[,i] - tp_nobiopsy[,i+1]  
       tn_nobiopsy[,i+1+j*n_combinations] <-  ((1 - pre_test_probability[,i]) * spec_HLA)
       fp_nobiopsy[,i+1+j*n_combinations] <- (1 - pre_test_probability[,i]) - tn_nobiopsy[,i+1]
     }
     
   }
   
   # Pre and post test odds and probabilities for serological tests being first
   # So these are for HLA second
   # Changed every occurence of 109 (which is 36 times 3 plus 1) to general (1 + n_combinations * 3)
   pre_test_probability_HLA <- array(0, dim=c(n_samples, n_combinations*3), dimnames = list(NULL, t_names[2:(1 + n_combinations * 3)]))
   
   for (i in 1:n_combinations) {
     pre_test_probability_HLA[,i] <- tp_nobiopsy[,i+1] /(tp_nobiopsy[,i+1] + fp_nobiopsy[,i+1])
     pre_test_probability_HLA[,i+n_combinations] <- tp_nobiopsy[,i+1+n_combinations] /(tp_nobiopsy[,i+1+n_combinations] + fp_nobiopsy[,i+1+n_combinations])
     pre_test_probability_HLA[,i+n_combinations+n_combinations] <- tp_nobiopsy[,i+1+n_combinations+n_combinations] /(tp_nobiopsy[,i+1+n_combinations+n_combinations] + fp_nobiopsy[,i+1+n_combinations+n_combinations])
   }
   
   pre_test_probability_HLA <- replace(pre_test_probability_HLA, pre_test_probability_HLA == 1, 0.99999)
   pre_test_odds_HLA <- array(0, dim=c(n_samples, n_combinations * 3), dimnames = list(NULL, t_names[2:(1 + n_combinations * 3)]))
   post_test_odds_HLA <- array(dim=c(n_samples, n_combinations * 3),dimnames=list(NULL, t_names[2:(1 + n_combinations * 3)])) 
   post_test_probability_HLA <- array(dim=c(n_samples, n_combinations * 3),dimnames=list(NULL, t_names[2:(1 + n_combinations * 3)]))
   
   
   for (i in 1:(3*n_combinations)){
     pre_test_odds_HLA[,i] <- pre_test_probability_HLA[,i]/(1 - pre_test_probability_HLA[,i])
     post_test_odds_HLA[,i] <- pre_test_odds_HLA[,i] * LR_HLA
     post_test_probability_HLA[,i] <- post_test_odds_HLA[,i]/(1 + post_test_odds_HLA[,i])
   }
   
   # Corresponding pre and post test odds and probability for HLA being first
   pre_test_probability_HLAfirst <- array(0, dim=c(n_samples, n_combinations * 3), dimnames = list(NULL, t_names[2:(1 + n_combinations * 3)]))
   for(j in 1:3) {
   for (i in 1:n_combinations) {
     pre_test_probability_HLAfirst[, i + (j - 1) * n_combinations] <- tp_nobiopsy[,i + 1 + n_combinations * 6] /(tp_nobiopsy[,i + 1+ n_combinations * 6] + fp_nobiopsy[, i + 1 + n_combinations * 6])
   }
   }
   
   pre_test_probability_HLAfirst <- replace(pre_test_probability_HLAfirst, pre_test_probability_HLAfirst == 1, 0.99999)
   pre_test_odds_HLAfirst <- array(0, dim=c(n_samples, n_combinations * 3), dimnames = list(NULL, t_names[(2 + n_combinations * 6):(1 + n_combinations * 9)]))
   post_test_odds_HLAfirst <- array(dim=c(n_samples, n_combinations * 3),dimnames=list(NULL, t_names[(2 + n_combinations * 6):(1 + n_combinations * 9)])) 
   post_test_probability_HLAfirst <- array(dim=c(n_samples, n_combinations * 3),dimnames=list(NULL, t_names[(2 + n_combinations * 6):(1 + n_combinations * 9)]))
   
   for (i in 1:(3*n_combinations)){
     pre_test_odds_HLAfirst[, i] <- pre_test_probability_HLAfirst[, i]/(1 - pre_test_probability_HLAfirst[, i])
   }
   # Test dependent
   for (i in 1:n_combinations) {
     post_test_odds_HLAfirst[, i] <- pre_test_odds_HLAfirst[, i] * LR_IgAEMA
     post_test_odds_HLAfirst[, i + n_combinations] <- pre_test_odds_HLAfirst[, i] * LR_IgATTGplusEMA
     post_test_odds_HLAfirst[, i + n_combinations * 2] <- pre_test_odds_HLAfirst[, i] * LR_IgATTG
   }
   for (i in 1:(3*n_combinations)){
     post_test_probability_HLAfirst[, i] <- post_test_odds_HLAfirst[, i]/(1 + post_test_odds_HLAfirst[, i])
   }
   
   for (i in 1:n_combinations) {
     
     # Probabilities for IgAEMA plus HLA
     tp[,i+1+(n_combinations*3)] <- ifelse(post_test_probability_HLA[,i] >= 0.9, (pre_test_probability_HLA[,i] * input_parameters$sens_HLA), 
                                           (pre_test_probability_HLA[,i] * input_parameters$sens_biopsy))
     fn[,i+1+(n_combinations*3)] <- pre_test_probability_HLA[,i] - tp[,i+1+(n_combinations*3)]  
     tn[,i+1+(n_combinations*3)] <- ifelse(post_test_probability_HLA[,i] >= 0.9, ((1 - pre_test_probability_HLA[,i]) * input_parameters$spec_HLA),
                                           ((1 - pre_test_probability_HLA[,i]) * input_parameters$spec_biopsy))
     fp[,i+1+(n_combinations*3)] <- (1 - pre_test_probability_HLA[,i]) - tn[,i+1+(n_combinations*3)]
     
     # Probabilities for IgATTG plus IgaEMA plus HLA
     tp[,i+1+(n_combinations*4)] <- ifelse(post_test_probability_HLA[,i + n_combinations] >= 0.9, (pre_test_probability_HLA[,i + n_combinations] * input_parameters$sens_HLA), 
                                           (pre_test_probability_HLA[,i + n_combinations] * input_parameters$sens_biopsy))
     fn[,i+1+(n_combinations*4)] <- pre_test_probability_HLA[,i + n_combinations] - tp[,i+1+(n_combinations*4)]  
     tn[,i+1+(n_combinations*4)] <- ifelse(post_test_probability_HLA[,i + n_combinations] >= 0.9, ((1 - pre_test_probability_HLA[,i + n_combinations]) * input_parameters$spec_HLA),
                                           ((1 - pre_test_probability_HLA[,i + n_combinations]) * input_parameters$spec_biopsy))
     fp[,i+1+(n_combinations*4)] <- (1 - pre_test_probability_HLA[,i + n_combinations]) - tn[,i+1+(n_combinations*4)]
     
     # Probabilities for IgATTG plus HLA
     tp[,i+1+(n_combinations*5)] <- ifelse(post_test_probability_HLA[,i + n_combinations + n_combinations] >= 0.9, (pre_test_probability_HLA[,i + n_combinations + n_combinations] * input_parameters$sens_HLA), 
                                           (pre_test_probability_HLA[,i + n_combinations + n_combinations] * input_parameters$sens_biopsy))
     fn[,i+1+(n_combinations*5)] <- pre_test_probability_HLA[,i + n_combinations + n_combinations] - tp[,i+1+(n_combinations*5)]  
     tn[,i+1+(n_combinations*5)] <- ifelse(post_test_probability_HLA[,i + n_combinations + n_combinations] >= 0.9, ((1 - pre_test_probability_HLA[,i + n_combinations + n_combinations]) * input_parameters$spec_HLA),
                                           ((1 - pre_test_probability_HLA[,i + n_combinations + n_combinations]) * input_parameters$spec_biopsy))
     fp[,i+1+(n_combinations*5)] <- (1 - pre_test_probability_HLA[,i + n_combinations + n_combinations]) - tn[,i+1+(n_combinations*5)]
    
     
     # Probabilities for HLA plus IgAEMA
     tp[,i+1+(n_combinations * 6)] <- ifelse(post_test_probability_HLAfirst[,i] >= 0.9, (pre_test_probability_HLAfirst[,i] * sens_IgAEMA), 
                                             (pre_test_probability_HLAfirst[,i] * input_parameters$sens_biopsy))
     fn[,i+1+(n_combinations * 6)] <- pre_test_probability_HLAfirst[,i] - tp[,i+1+(n_combinations * 6)]  
     tn[,i+1+(n_combinations * 6)] <- ifelse(post_test_probability_HLAfirst[,i] >= 0.9, ((1 - pre_test_probability_HLAfirst[,i]) * spec_IgAEMA),
                                             ((1 - pre_test_probability_HLAfirst[,i]) * input_parameters$spec_biopsy))
     fp[,i+1+(n_combinations * 6)] <- (1 - pre_test_probability_HLAfirst[,i]) - tn[,i+1+(n_combinations * 6)]
     
     # Probabilities for HLA plus IgATTG plus IgaEMA
     tp[,i+1+(n_combinations * 7)] <- ifelse(post_test_probability_HLAfirst[,i + n_combinations] >= 0.9, (pre_test_probability_HLAfirst[,i + n_combinations] * sens_IgATTGplusEMA), 
                                             (pre_test_probability_HLAfirst[,i + n_combinations] * input_parameters$sens_biopsy))
     fn[,i+1+(n_combinations * 7)] <- pre_test_probability_HLAfirst[,i + n_combinations] - tp[,i+1+(n_combinations * 7)]  
     tn[,i+1+(n_combinations * 7)] <- ifelse(post_test_probability_HLAfirst[,i + n_combinations] >= 0.9, ((1 - pre_test_probability_HLAfirst[,i + n_combinations]) * spec_IgATTGplusEMA),
                                             ((1 - pre_test_probability_HLAfirst[,i + n_combinations]) * input_parameters$spec_biopsy))
     fp[,i+1+(n_combinations * 7)] <- (1 - pre_test_probability_HLAfirst[,i + n_combinations]) - tn[,i+1+(n_combinations * 7)]
     
     # Probabilities for HLA plus IgATTG
     tp[,i+1+(n_combinations * 8)] <- ifelse(post_test_probability_HLAfirst[,i + n_combinations + n_combinations] >= 0.9, (pre_test_probability_HLAfirst[,i + n_combinations + n_combinations] * sens_IgATTG), 
                                             (pre_test_probability_HLAfirst[,i + n_combinations + n_combinations] * input_parameters$sens_biopsy))
     fn[,i+1+(n_combinations * 8)] <- pre_test_probability_HLAfirst[,i + n_combinations + n_combinations] - tp[,i+1+(n_combinations * 8)]  
     tn[,i+1+(n_combinations * 8)] <- ifelse(post_test_probability_HLAfirst[,i + n_combinations + n_combinations] >= 0.9, ((1 - pre_test_probability_HLAfirst[,i + n_combinations + n_combinations]) * spec_IgATTG),
                                             ((1 - pre_test_probability_HLAfirst[,i + n_combinations + n_combinations]) * input_parameters$spec_biopsy))
     fp[,i+1+(n_combinations * 8)] <- (1 - pre_test_probability_HLAfirst[,i + n_combinations + n_combinations]) - tn[,i+1+(n_combinations * 8)]
     
   }
   
   # The following can be used for bug checking
   #write.csv(colMeans(post_test_probability_IgAEMA), "results/posttestIgAema.csv")
   #write.csv(colMeans(post_test_probability_IgATTGplusEMA), "results/posttestIgAttgplusema.csv")
   #write.csv(colMeans(post_test_probability_IgATTG), "results/posttestIgAttg.csv")
   #write.csv(colMeans(pre_test_probability_HLA), "results/hlapretest.csv")
   #write.csv(colMeans(post_test_probability_HLA), "results/hlaposttest.csv")
   #write.csv(data.frame(colMeans(tp), colMeans(fn), colMeans(tn), colMeans(fp)), "results/test.csv")
   
   #Adding costs   
   fp_costs <- array(dim = c(n_samples, n_tests),
                     dimnames = list(NULL, t_names))
   
   fp_costs[,] <- (fp[,] * input_parameters$cost_gfp) 
   
   diagnosis_costs <- array(dim = c(n_samples, n_tests),
                            dimnames = list(NULL, t_names))
   
   diagnosis_costs[,] <- (tp[,] + fp[,]) * input_parameters$cost_diagnosis  #cost of diagnosis is cost of nurse/GP/specialist consultations to get to a diagnosis
   
   
   test_costs <- array(dim=c(n_samples, n_tests), dimnames=list(NULL, t_names))
   for (i in 1:n_combinations) {
     test_costs[, 1] <- 0
     
     # EMA
     test_costs[, i+1] <-  ifelse(post_test_probability_IgAEMA[,i] >= 0.9, input_parameters$test_cost_IgAEMA, 
                                  input_parameters$test_cost_IgAEMA + input_parameters$cost_biopsy)
     # TTG plus EMA
     test_costs[, n_combinations+i+1] <-  ifelse(post_test_probability_IgATTGplusEMA[,i] >= 0.9, (input_parameters$test_cost_IgATTG + input_parameters$test_cost_IgAEMA), 
                                                 (input_parameters$test_cost_IgATTG + input_parameters$test_cost_IgAEMA + input_parameters$cost_biopsy))
     # TTG
     test_costs[, n_combinations+n_combinations+i+1] <-  ifelse(post_test_probability_IgATTG[,i] >= 0.9, input_parameters$test_cost_IgATTG, 
                                                                input_parameters$test_cost_IgATTG + input_parameters$cost_biopsy)
     # EMA plus HLA
     test_costs[, n_combinations+n_combinations+ n_combinations+i+1] <- 
       ifelse(post_test_probability_IgAEMA[,i] >= 0.9,
              input_parameters$test_cost_IgAEMA,
              ifelse(post_test_probability_IgAEMA[,i] < 0.9 & post_test_probability_HLA[,i] >= 0.9,
                     (input_parameters$test_cost_IgAEMA + input_parameters$test_cost_HLA),
                     (input_parameters$test_cost_IgAEMA + input_parameters$test_cost_HLA + input_parameters$cost_biopsy)))
     
     # TTG plus EMA plus HLA
     test_costs[, n_combinations + n_combinations+ n_combinations + n_combinations + i + 1] <-  
       ifelse(post_test_probability_IgATTGplusEMA[,i] >= 0.9,
              input_parameters$test_cost_IgAEMA + input_parameters$test_cost_IgATTG ,
              ifelse(post_test_probability_IgATTGplusEMA[,i] < 0.9 & post_test_probability_HLA[,i+n_combinations] >= 0.9, 
                     input_parameters$test_cost_IgAEMA + input_parameters$test_cost_IgATTG  + input_parameters$test_cost_HLA, 
                     input_parameters$test_cost_IgAEMA + input_parameters$test_cost_IgATTG + input_parameters$test_cost_HLA + input_parameters$cost_biopsy))

     # TTG plus HLA     
     test_costs[, n_combinations + n_combinations+ n_combinations + n_combinations+ n_combinations+i+1] <-  
       ifelse(post_test_probability_IgATTG[,i] >= 0.9,
              input_parameters$test_cost_IgATTG,
              ifelse(post_test_probability_IgATTG[,i] < 0.9 & post_test_probability_HLA[,i+n_combinations+n_combinations] >= 0.9, 
                     input_parameters$test_cost_IgATTG + input_parameters$test_cost_HLA, 
                     input_parameters$test_cost_IgATTG + input_parameters$test_cost_HLA + input_parameters$cost_biopsy))
    
     # HLA plus EMA
     test_costs[, (n_combinations * 6) + i + 1] <- 
       ifelse(post_test_probability_HLAalone[,i] >= 0.9, 
              input_parameters$test_cost_HLA,
              ifelse(post_test_probability_HLAalone[,i] < 0.9 & post_test_probability_HLAfirst[,i] >= 0.9,
                     (input_parameters$test_cost_HLA + input_parameters$test_cost_IgAEMA),
                     (input_parameters$test_cost_HLA + input_parameters$test_cost_IgAEMA + input_parameters$cost_biopsy)))
     
     # HLA plus TTG plus EMA
     test_costs[, (n_combinations * 7) + i + 1] <-  
       ifelse(post_test_probability_HLAalone[,i] >= 0.9, # This should be corrected to post test probability using HLA alone
              input_parameters$test_cost_HLA,
              ifelse(post_test_probability_HLAalone[,i] < 0.9 & post_test_probability_HLAfirst[,i+n_combinations] >= 0.9, 
                     input_parameters$test_cost_IgAEMA + input_parameters$test_cost_IgATTG  + input_parameters$test_cost_HLA, 
                     input_parameters$test_cost_IgAEMA + input_parameters$test_cost_IgATTG + input_parameters$test_cost_HLA + input_parameters$cost_biopsy))
     
     # HLA plus TTG
     test_costs[, (n_combinations * 8) + i + 1] <-  
       ifelse(post_test_probability_HLAalone[,i] >= 0.9, # This should be corrected to post test probability using HLA alone
              input_parameters$test_cost_HLA,
              ifelse(post_test_probability_HLAalone[,i] < 0.9 & post_test_probability_HLAfirst[,i+n_combinations+n_combinations] >= 0.9, 
                     input_parameters$test_cost_IgATTG + input_parameters$test_cost_HLA, 
                     input_parameters$test_cost_IgATTG + input_parameters$test_cost_HLA + input_parameters$cost_biopsy))
     
   }
   
   biopsy_costs <- array(dim=c(n_samples, n_tests), dimnames=list(NULL, t_names))
   for (i in 1:n_combinations) {
     biopsy_costs[, 1] <- 0
     # EMA
     biopsy_costs[, i+1] <-  ifelse(post_test_probability_IgAEMA[,i] >= 0.9, 0, input_parameters$cost_biopsy)
     # TTG plus EMA
     biopsy_costs[, n_combinations+i+1] <-  ifelse(post_test_probability_IgATTGplusEMA[,i] >= 0.9, 0, input_parameters$cost_biopsy)
     # TTG
     biopsy_costs[, n_combinations+n_combinations+i+1] <-  ifelse(post_test_probability_IgATTG[,i] >= 0.9, 0, input_parameters$cost_biopsy)
     # EMA plus HLA
     biopsy_costs[, n_combinations+n_combinations+ n_combinations+i+1] <- 
       ifelse(post_test_probability_IgAEMA[,i] >= 0.9,
              0,
              ifelse(post_test_probability_IgAEMA[,i] < 0.9 & post_test_probability_HLA[,i] >= 0.9,
                     0,
                     (input_parameters$cost_biopsy)))
     # TTG plus EMA plus HLA
     biopsy_costs[, n_combinations + n_combinations+ n_combinations + n_combinations + i + 1] <-  
       ifelse(post_test_probability_IgATTGplusEMA[,i] >= 0.9,
              0,
              ifelse(post_test_probability_IgATTGplusEMA[,i] < 0.9 & post_test_probability_HLA[,i+n_combinations] >= 0.9, 
                     0, 
                     input_parameters$cost_biopsy))
     
     # TTG plus HLA
     biopsy_costs[, n_combinations + n_combinations+ n_combinations + n_combinations+ n_combinations+i+1] <-  
       ifelse(post_test_probability_IgATTG[,i] >= 0.9,
              0,
              ifelse(post_test_probability_IgATTG[,i] < 0.9 & post_test_probability_HLA[,i+n_combinations+n_combinations] >= 0.9, 
                     0, 
                     input_parameters$cost_biopsy))

     # HLA plus EMA
     biopsy_costs[, (n_combinations * 6) + i + 1] <- 
       ifelse(post_test_probability_HLAalone[,i] >= 0.9, # This should be corrected to post test probability using HLA alone
              0,
              ifelse(post_test_probability_HLAalone[,i] < 0.9 & post_test_probability_HLAfirst[,i] >= 0.9,
                     0,
                     (input_parameters$cost_biopsy)))
     # HLA plus TTG plus EMA
     biopsy_costs[, (n_combinations * 7) + i + 1] <-  
       ifelse(post_test_probability_HLAalone[,i] >= 0.9,
              0,
              ifelse(post_test_probability_HLAalone[,i] < 0.9 & post_test_probability_HLAfirst[,i+n_combinations] >= 0.9, 
                     0, 
                     input_parameters$cost_biopsy))
     
     # HLA plus TTG
     biopsy_costs[, (n_combinations * 8) + i + 1] <-  
       ifelse(post_test_probability_HLAalone[,i] >= 0.9,
              0,
              ifelse(post_test_probability_HLAalone[,i] < 0.9 & post_test_probability_HLAfirst[,i+n_combinations+n_combinations] >= 0.9, 
                     0, 
                     input_parameters$cost_biopsy))
     
     
   }
   
   #Adding disutilities
   
   disutility_biopsy_screen <- array(dim=c(n_samples,n_tests),dimnames=list(NULL, t_names))
   for (i in 1:n_combinations) {
     disutility_biopsy_screen[, 1] <- 0
     # Sero tests alone
     disutility_biopsy_screen[, i+1] <-  ifelse(post_test_probability_IgAEMA[,i] >= 0.9, 0, input_parameters$disutility_biopsy)
     disutility_biopsy_screen[, n_combinations+i+1] <-  ifelse(post_test_probability_IgATTGplusEMA[,i] >= 0.9, 0, input_parameters$disutility_biopsy)
     disutility_biopsy_screen[, n_combinations+n_combinations+i+1] <-  ifelse(post_test_probability_IgATTG[,i] >= 0.9, 0, input_parameters$disutility_biopsy)
     # Sero plus HLA strategies
     disutility_biopsy_screen[, n_combinations+n_combinations+ n_combinations+i+1] <-  ifelse(post_test_probability_HLA[,i] >= 0.9, 0, input_parameters$disutility_biopsy)
     disutility_biopsy_screen[, n_combinations + n_combinations+ n_combinations + n_combinations +i+1] <-  ifelse(post_test_probability_HLA[,i+n_combinations] >= 0.9, 0, input_parameters$disutility_biopsy)
     disutility_biopsy_screen[, n_combinations + n_combinations+ n_combinations + n_combinations+ n_combinations+i+1] <-  ifelse(post_test_probability_HLA[,i+n_combinations+n_combinations] >= 0.9, 0, input_parameters$disutility_biopsy)
     # HLA plus EMA
     disutility_biopsy_screen[, (n_combinations * 6) + i + 1] <-  ifelse(post_test_probability_HLAfirst[,i] >= 0.9, 0, input_parameters$disutility_biopsy)
     # HLA plus TTG plus EMA
     disutility_biopsy_screen[, (n_combinations * 7) + i + 1] <-  ifelse(post_test_probability_HLAfirst[,i+n_combinations] >= 0.9, 0, input_parameters$disutility_biopsy)
     # HLS plus TTG
     disutility_biopsy_screen[, (n_combinations * 8) + i + 1] <-  ifelse(post_test_probability_HLAfirst[,i+n_combinations+n_combinations] >= 0.9, 0, input_parameters$disutility_biopsy)
     
   }
   
   
   disutility_biopsy_screen_wait <- array(dim=c(n_samples, n_tests),dimnames=list(NULL, t_names))
   for (i in 1:n_combinations) {
     disutility_biopsy_screen_wait[,1] <- 0
     # Sero tests
     disutility_biopsy_screen_wait[, i+1] <-  ifelse(post_test_probability_IgAEMA[,i] >= 0.9, 0, input_parameters$disutility_biopsy_wait)
     disutility_biopsy_screen_wait[, n_combinations+i+1] <-  ifelse(post_test_probability_IgATTGplusEMA[,i] >= 0.9, 0, input_parameters$disutility_biopsy_wait)
     disutility_biopsy_screen_wait[, n_combinations+n_combinations+i+1] <-  ifelse(post_test_probability_IgATTG[,i] >= 0.9, 0, input_parameters$disutility_biopsy_wait)
     # Sero plus HLA strategies
     disutility_biopsy_screen_wait[, n_combinations+n_combinations+ n_combinations+i+1] <-  ifelse(post_test_probability_HLA[,i] >= 0.9, 0, input_parameters$disutility_biopsy_wait)
     disutility_biopsy_screen_wait[, n_combinations + n_combinations+ n_combinations + n_combinations +i+1] <-  ifelse(post_test_probability_HLA[,i+n_combinations] >= 0.9, 0, input_parameters$disutility_biopsy_wait)
     disutility_biopsy_screen_wait[, n_combinations + n_combinations+ n_combinations + n_combinations+ n_combinations+i+1] <-  ifelse(post_test_probability_HLA[,i+n_combinations+n_combinations] >= 0.9, 0, input_parameters$disutility_biopsy_wait)
     # HLA plus EMA
     disutility_biopsy_screen_wait[, (n_combinations * 6) + i + 1] <-  ifelse(post_test_probability_HLAfirst[,i] >= 0.9, 0, input_parameters$disutility_biopsy_wait)
     # HLA plus TTG plus EMA
     disutility_biopsy_screen_wait[, (n_combinations * 7) + i + 1] <-  ifelse(post_test_probability_HLAfirst[,i+n_combinations] >= 0.9, 0, input_parameters$disutility_biopsy_wait)
     # HLA plus TTG
     disutility_biopsy_screen_wait[, (n_combinations * 8) + i + 1] <-  ifelse(post_test_probability_HLAfirst[,i+n_combinations+n_combinations] >= 0.9, 0, input_parameters$disutility_biopsy_wait)
   }
   
   disutility_fp <- array(dim = c(n_samples, n_tests),
                          dimnames = list(NULL, t_names))
   
   disutility_fp[,] <- (fp[,] * input_parameters$disutility_fp) 
   
   # True positives etc due to risk prediction are the same across testing strategies
   tp_riskfactor_table <- array(dim=c(n_samples, n_tests), dimnames = list(NULL, t_names))
   fp_riskfactor_table <- array(dim=c(n_samples, n_tests), dimnames = list(NULL, t_names))
   fn_riskfactor_table <- array(dim=c(n_samples, n_tests), dimnames = list(NULL, t_names))
   for (i in 1:n_combinations) {
     tp_riskfactor_table[,1] <- 1
     fp_riskfactor_table[,1] <- 1
     fn_riskfactor_table[,1] <- 1
     for(j in 1:length(tests)) {
       tp_riskfactor_table[,i + 1 + (j - 1) * n_combinations] <- tp_riskfactor[, i]
       fp_riskfactor_table[,i + 1 + (j - 1) * n_combinations] <- fp_riskfactor[, i]
       fn_riskfactor_table[,i + 1 + (j - 1) * n_combinations] <- fn_riskfactor[, i]
     }
   }
   
   test_costs_applied <- array(dim = c(n_samples, n_tests),
                               dimnames = list(NULL, t_names))
   
   test_costs_applied <- test_costs * (tp[,] + fp[,] + fn[,] + tn[,] + tp_riskfactor_table[,] + fp_riskfactor_table[,])
   false_positive_costs_applied <- test_costs * (fn[,] + tn[,]) #cost of testing in those without CD after risk factor test
   biopsy_disutility_applied <- disutility_biopsy_screen * (tp[,] + fp[,] + fn[,] + tn[,] + tp_riskfactor_table[,] + fp_riskfactor_table[,])
   biopsy_Wait_disutility_applied <- disutility_biopsy_screen_wait * (tp[,] + fn[,])
   
   
   fn_riskfactor_table <- fn_riskfactor_table * 1/(fp+tp)
   fn_all <- fn + fn_riskfactor_table
   
   
   #scaling up true positives and false negatives 
   tp <- tp/(tp + fn_all)
   fn <- 1 - tp
   # tp <- array(0, dim=c(n_samples, n_tests), dimnames = list(NULL, t_names))  
   #fn <- array(1, dim=c(n_samples, n_tests), dimnames = list(NULL, t_names))
   
   # Build an array to store the cohort vector at each cycle
   # Each cohort vector has n_states elements: probability of being in each state,
   # There is one cohort vector for each test, for each PSA sample, for each cycle.
   cohort_vectors<-array(dim=c(n_tests,n_samples,n_cycles,n_states),  
                         dimnames=list(t_names,NULL,NULL, state_names))
   
   
   for (i_test in c(1:n_tests)) { 
     cohort_vectors[i_test, , 1, "CD GFD no complications"] <- tp[, i_test] * input_parameters$probability_nocomplications 
     cohort_vectors[i_test, , 1, "CD GFD osteoporosis"] <- tp[, i_test] * input_parameters$probability_osteoporosis 
     cohort_vectors[i_test, , 1, "CD GFD NHL"] <- tp[, i_test] * input_parameters$probability_NHL
     
     cohort_vectors[i_test, , 1, "Undiagnosed CD no complications"] <- fn[, i_test] * input_parameters$probability_nocomplications 
     cohort_vectors[i_test, , 1, "Undiagnosed CD osteoporosis"] <- fn[, i_test] * input_parameters$probability_osteoporosis 
     cohort_vectors[i_test, , 1, "Undiagnosed CD NHL"] <- fn[, i_test] * input_parameters$probability_NHL 
     
     cohort_vectors[i_test, , 1, "Death"] <- 0
   }
   
   cohort_vectors[, , , ] [is.na(cohort_vectors[, , , ] )] <- 0
   # Checks to ensure making sense
   #colSums (cohort_vectors[, 2, , ], na.rm = FALSE, dims = 1)
   #rowSums (cohort_vectors[, 2, , ], na.rm = FALSE, dims = 1)
   
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
   
   
   disc_vec <- (1/1.035) ^ rep(c(0:(n_cycles-1)))
   
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
         cycle_costs[i_test, i_sample, i_cycle] <-
           cohort_vectors[i_test, i_sample, i_cycle, ] %*% state_costs[i_sample, i_cycle, ]
         # And total QALYs for each cycle
         cycle_qalys[i_test, i_sample, i_cycle] <-
           cohort_vectors[i_test, i_sample, i_cycle, ]  %*% state_qalys[i_sample, i_cycle, ]
         
       }
       
       # Combine the cycle_costs and other costs to get total costs
       # Apply the discount factor 
       # (1 in first year, 1.035 in second, 1.035^2 in third, and so on)
       total_costs[i_test, i_sample] <- test_costs_applied[i_sample, i_test] + fp_costs[i_sample, i_test] + diagnosis_costs[i_sample, i_test] + 
         cycle_costs[i_test, i_sample, ] %*% disc_vec
       
       # total_costs[i_test, i_sample] <-  cycle_costs[i_test, i_sample, ] %*% disc_vec #for checking all FPs/TNs
       
       # Combine the cycle_qalys to get total qalys
       # Apply the discount factor 
       # (1 in first year, 1.035 in second, 1.035^2 in third, and so on)
       total_qalys[i_test, i_sample] <- (cycle_qalys[i_test, i_sample, ]  - disutility_fp[i_sample, i_test] - biopsy_disutility_applied[i_sample, i_test] - biopsy_Wait_disutility_applied[i_sample, i_test])  %*% disc_vec
       
       # total_qalys[i_test, i_sample] <- (cycle_qalys[i_test, i_sample, ] 
       #%*% disc_vec)  #for checking all FPs/TNs
     }
     
   }
   
   
   
   
   
   #############################################################################
   ## Analysis of results ######################################################
   #############################################################################
   output <- list()
   
   
   
   #Percentage having biopsy
   sum_IgAEMA <- rep(0, n_combinations)
   for (i in 1:n_combinations){
     sum_IgAEMA[i] <- sum(post_test_probability_IgAEMA[,i] <= 0.9)/n_samples  
   }
   names(sum_IgAEMA) <- t_names[1 + 1:n_combinations]
   
   sum_IgATTGplusEMA <- rep(0, n_combinations)
   for (i in 1:n_combinations){
     sum_IgATTGplusEMA[i] <- sum(post_test_probability_IgATTGplusEMA[,i] <= 0.9)/n_samples  
   }
   names(sum_IgATTGplusEMA) <- t_names[1 + (1 * n_combinations) + (1:n_combinations)]
   
   sum_IgATTG <- rep(0, n_combinations)
   for (i in 1:n_combinations){
     sum_IgATTG[i] <- sum(post_test_probability_IgATTG[,i] <= 0.9)/n_samples  
   }
   names(sum_IgATTG) <- t_names[1 + (2 * n_combinations) + (1:n_combinations)]
   
   # Proportion getting biopsy or HLA on sero plus HLA strategies
   sum_IgAEMAplusHLA <- rep(0, n_combinations)
   sum_IgAEMAplusHLA_test <- rep(0, n_combinations)
   for (i in 1:n_combinations){
     sum_IgAEMAplusHLA[i] <- sum(post_test_probability_HLA[,i] <= 0.9)/n_samples  
     sum_IgAEMAplusHLA_test[i] <- sum(post_test_probability_IgAEMA[,i] <= 0.9)/n_samples
   }
   names(sum_IgAEMAplusHLA) <- names(sum_IgAEMAplusHLA_test) <- t_names[1 + (3 * n_combinations) + (1:n_combinations)]
   
   sum_IgATTGplusEMAplusHLA <- rep(0, n_combinations)
   sum_IgATTGplusEMAplusHLA_test <- rep(0, n_combinations)
   for (i in 1:n_combinations){
     sum_IgATTGplusEMAplusHLA[i] <- sum(post_test_probability_HLA[,i+n_combinations] <= 0.9)/n_samples 
     sum_IgATTGplusEMAplusHLA_test[i] <- sum(post_test_probability_IgATTGplusEMA[,i] <= 0.9)/n_samples
   }
   names(sum_IgATTGplusEMAplusHLA) <- names(sum_IgATTGplusEMAplusHLA_test) <- t_names[1 + (4 * n_combinations) + (1:n_combinations)]
   
   sum_IgATTGplusHLA <- rep(0, n_combinations)
   sum_IgATTGplusHLA_test <- rep(0, n_combinations)
   for (i in 1:n_combinations){
     sum_IgATTGplusHLA[i] <- sum(post_test_probability_HLA[,i+n_combinations+n_combinations] <= 0.9)/n_samples  
     sum_IgATTGplusHLA_test[i] <- sum(post_test_probability_IgATTG[,i] <= 0.9)/n_samples  
     
   }
   names(sum_IgATTGplusHLA) <-  names(sum_IgATTGplusHLA_test) <- t_names[1 + (5 * n_combinations) + (1:n_combinations)]
   
   # Proportion getting biopsy on HLA plus Sero strategies
   # HLA plus EMA
   sum_HLAplusIgAEMA <- rep(0, n_combinations)
   for (i in 1:n_combinations){
     sum_HLAplusIgAEMA[i] <- sum(post_test_probability_HLAfirst[,i] <= 0.9)/n_samples  
   }
   names(sum_HLAplusIgAEMA ) <- t_names[1 + (6 * n_combinations) + (1:n_combinations)]
   
   # HLA plus TTG plus EMA
   sum_HLAplusIgATTGplusEMA <- rep(0, n_combinations)
   for (i in 1:n_combinations){
     sum_HLAplusIgATTGplusEMA[i] <- sum(post_test_probability_HLAfirst[,i+n_combinations] <= 0.9)/n_samples 
   }
   names(sum_HLAplusIgATTGplusEMA) <- t_names[1 + (7 * n_combinations) + (1:n_combinations)]
   
   # HLA plus TTG
   sum_HLAplusIgATTG <- rep(0, n_combinations)
   for (i in 1:n_combinations){
     sum_HLAplusIgATTG[i] <- sum(post_test_probability_HLAfirst[,i+n_combinations+n_combinations] <= 0.9)/n_samples  
   }
   names(sum_HLAplusIgATTG) <- t_names[1 + (8 * n_combinations) + (1:n_combinations)]
   
   
   output$percentage_biopsy_IgAEMA <- sum_IgAEMA
   output$percentage_biopsy_IgATTGplusEMA <- sum_IgATTGplusEMA
   output$percentage_biopsy_IgATTG <- sum_IgATTG
   output$percentage_biopsy_IgAEMAplusHLA <- sum_IgAEMAplusHLA
   output$percentage_hla_IgAEMAplusHLA <- sum_IgAEMAplusHLA_test
   output$percentage_biopsy_IgATTGplusEMAplusHLA <- sum_IgATTGplusEMAplusHLA
   output$percentage_hla_IgATTGplusEMAplusHLA <- sum_IgATTGplusEMAplusHLA_test
   output$percentage_biopsy_IgATTGplusHLA <- sum_IgATTGplusHLA
   output$percentage_hla_IgATTGplusHLA <- sum_IgATTGplusHLA_test
   
   # For HLA plus sero strategies
   output$percentage_biopsy_HLAplusIgAEMA<- sum_HLAplusIgAEMA
   output$percentage_biopsy_HLAplusIgATTGplusEMA <- sum_HLAplusIgATTGplusEMA
   output$percentage_biopsy_HLAplusIgATTG <- sum_HLAplusIgATTG
   
   write.csv(sum_IgATTGplusEMAplusHLA_test, "x.csv")
   write.csv(sum_IgAEMAplusHLA, "x.csv")
   write.csv(sum_IgATTGplusEMAplusHLA, "x.csv")
   
   #Costs and QALYs
   output$total_costs <- total_costs
   output$total_qalys <- total_qalys
   # Average costs
   output$average_costs <- rowMeans(total_costs)
   # Average effects (in QALY units)
   output$average_effects <- rowMeans(total_qalys)
   
   #for all TPs/FNs comparison
   quantile(output$total_costs, 0.025)
   quantile(output$total_costs, 0.975)
   quantile(output$total_qalys, 0.025)
   quantile(output$total_qalys, 0.975)
   
   output$incremental_costs <-  output$average_costs - output$average_costs["No screening"]
   output$incremental_effects <-  output$average_effects -  output$average_effects["No screening"]
   
   output$all_incremental_costs <-   output$all_incremental_effects <- array(dim=c(n_tests,n_samples),  
                                                                             dimnames=list(t_names,NULL))
   
   for (i_sample in 1:n_samples) {
     for (i_test in 1:n_tests) {
       output$all_incremental_costs[i_test, ] <-  output$total_costs[i_test,] - output$total_costs["No screening",]
       output$all_incremental_effects[i_test, ] <-  output$total_qalys[i_test,] -  output$total_qalys["No screening",]
     }}
   
   
   output$ICER <- output$incremental_costs/output$incremental_effects
   
   # Incremental net benefit at the ?20,000 willingness-to-pay
   
   output$incremental_net_benefit <- 20000*output$incremental_effects - output$incremental_costs
   
   output$net_benefit<- 20000*output$average_effects - output$average_costs
   output$all_net_benefit <- 20000*total_qalys - total_costs
   
   quantile(output$all_net_benefit, 0.025)
   quantile(output$all_net_benefit, 0.975)
   
   output$ceac_calculation <- array(dim=c(n_tests,n_samples),  
                                    dimnames=list(t_names,NULL))
   
   for (i_sample in 1:n_samples) {
     for (i_test in 1:n_tests) {
       output$ceac_calculation[i_test,i_sample] <- ifelse(output$all_net_benefit[i_test,i_sample] == max(output$all_net_benefit[,i_sample]), 1, 0)
     }}
   
   output$probability_best <- rowMeans(output$ceac_calculation)
   
   #time spent in states
   dim(apply(cohort_vectors, c(1, 2, 4), sum))
   x<-apply(cohort_vectors, c(1, 2, 4), sum)
   x[1,1,]/n_cycles #total time spent in each state over all cycles
   time_in_states <- apply(x, c(1, 3), mean)/n_cycles
   apply(x, c(1, 3), quantile, probs = 0.025)
   
   if(!is.null(strategies_of_interest)) {
     output$time_in_states <- time_in_states[strategies_of_interest, ]
   } else {
     output$time_in_states <- NULL
   }
   
   
   
   #cost breakdown
   output$test_costs_applied <- test_costs_applied
   output$test_costs <- colMeans(test_costs_applied) #costs of test and biopsies
   output$false_positive_costs_applied <- false_positive_costs_applied
   output$fp_costs <- fp_costs
   output$diagnosis_costs <- diagnosis_costs
   output$cycle_costs <- rowMeans(cycle_costs)
   
   
   #utility breakdown
   output$cycle_qalys <- rowMeans(cycle_qalys)
   output$disutility_biopsy <- biopsy_disutility_applied
   output$disutility_biopsy_wait <- biopsy_Wait_disutility_applied
   output$disutility_fp <- disutility_fp
   
   
   
   # Average incremental net benefit
   #output$average_inb_IgATTGplusIgAEMA_IgAEMA <- mean(output$incremental_net_benefit_IgATTGplusIgAEMA_IgAEMA)
   #output$average_inb_doubletest_IgAEMA <- mean(output$incremental_net_benefit_doubletest_IgAEMA)
   
   # Probability cost-effective
   # This is the proportion of samples for which the incremental net benefit is positive
   output$all_incremental_net_benefit <- 20000*output$all_incremental_effects - output$all_incremental_costs
   output$incremental_net_benefit <- rowMeans(output$all_incremental_net_benefit)
   
   
   output$inb_lci <-   output$inb_uci <- array(dim=c(n_tests,n_samples),  
                                               dimnames=list(t_names,NULL))
   
   for (i_sample in 1:n_samples) {
     for (i_test in 1:n_tests) {
       output$inb_lci[i_test, i_sample] <- quantile(output$all_incremental_net_benefit[i_test,], 0.025)
       output$inb_uci[i_test, i_sample] <- quantile(output$all_incremental_net_benefit[i_test,], 0.975)
       
     }
   }
   
   output$inb_lci <- rowMeans(output$inb_lci)
   
   output$inb_uci <- rowMeans(output$inb_uci)
   
   output$pce_calculation <- array(dim=c(n_tests,n_samples),  
                                   dimnames=list(t_names,NULL))
   
   for (i_sample in 1:n_samples) {
     for (i_test in 1:n_tests) {
       output$pce_calculation[i_test,i_sample] <- sum(output$all_incremental_net_benefit[i_test,] > 0)/n_samples
     }}
   
   output$probability_cost_effective <- rowMeans(output$pce_calculation)
   #output$probability_cost_effective_doubletest_IgAEMA <- sum(output$incremental_net_benefit_doubletest_IgAEMA > 0)/n_samples
   output$pce_lci <-   output$pce_uci <- array(dim=c(n_tests,n_samples),  
                                               dimnames=list(t_names,NULL))
   
   for (i_sample in 1:n_samples) {
     for (i_test in 1:n_tests) {
       output$pce_lci[i_test, i_sample] <- quantile(output$pce_calculation[i_test,], 0.025)
       output$pce_uci[i_test, i_sample] <- quantile(output$pce_calculation[i_test,], 0.975)
       
     }
   }
   
   return(output)
}

