# Coeliac disease Markov model
# Edna Keeney

# Notes to Edna:
# Left disutility of NHL as an annual decrement as it represents maintenance therapy
# TODO: Add INB to time_in_states (line ~524)
# TODO: info.rank needs to plot only top parameters (line ~584) as most have zero impact
# and only include sens/spec for feasible strategies

# Used C++, lapply and vectorisation to optimise generate_net_benefit. Decreased runtime by ~50%

library(tictoc)
library(BCEA)
library(SimDesign)
library(BCEA)
library(dplyr)
library(EnvStats)
library(Rmisc)
library(ggplot2)
library(readxl)
library(car)
library(triangle)
library(mvtnorm)

#setwd("C:/Users/ek14588/Downloads/Coeliac_screening")


rm(list=ls())

# Define simulation parameters
# This is the number of PSA samples to use
n_samples <- 1000

perspective <- "NHS" #Options are "NHS" or "NHS+OOP" if out-of-pocket costs for iron supplements and gluten free products are to be included

# Set to FALSE if don't want to run MLMC
run_mlmc_evppi <- FALSE
#########################################################################################################################
source("code/generate_state_costs.R")
source("code/generate_state_qalys.R")
source("code/generate_model_parameters.R")
source("code/generate_transition_matrices.R")
source("code/convert_transition_matrices_to_df.R")
source("code/convert_cohort_vectors_to_df.R")
source("code/generate_net_benefit_hla-HT.R")
source("code/generate_net_benefit_cpp.R")

# MLMC scripts
source("code/mlmc.R")
source("code/mlmc.test.R")




#population <- "children" #Options are "men" or "women" or "children"
populations <- c("children", "men" , "women")


tests <- c("IgAEMA", "IgATTGplusEMA", "IgATTG", 
           "IgAEMA plus HLA", "IgATTGplusEMA plus HLA", "IgATTG plus HLA",
           "HLA plus IgAEMA", "HLA plus IgATTGplusEMA", "HLA plus IgATTG")
n_sero_tests <- length(tests)

for(population in populations) {
  # Set the seed for each population run to be the same
  set.seed(14143)
  
  print(population)
  #pre-test probabilities of coeliac disease 
  sens_riskfactor <- c(0.5, 0.6, 0.7, 0.8, 0.9, 0.9999)
  spec_riskfactor <- c(0.5, 0.6, 0.7, 0.8, 0.9, 0.9999)
  combinations <- expand.grid(sens_riskfactor = sens_riskfactor, spec_riskfactor = spec_riskfactor)
  # Add specific combinations from the risk prediction tool
  included_strategies <- as.data.frame(read.csv("data/included_strategies.csv"))
  new_combinations <- as.data.frame(matrix(nrow = 6, ncol = 2))
  colnames(new_combinations) <- c("sens_riskfactor", "spec_riskfactor")
  new_combinations$sens_riskfactor <- included_strategies[included_strategies$Population == population , "Sens"]
  new_combinations$spec_riskfactor <- included_strategies[included_strategies$Population == population, "Spec"]
  new_combinations$sens_riskfactor[new_combinations$sens_riskfactor == 1] <- 0.9999
  new_combinations$spec_riskfactor[new_combinations$spec_riskfactor == 1] <- 0.9999
  new_combinations <- unique(new_combinations)
  combinations <- rbind(combinations, new_combinations)
  
  
  # Used for naming the risk prediction strategies
  risk_percentages <- c("1%", "1.5%", "2%", "5%", "10%", "20%")
  
  # Specify the strategies of interest based on prediction modelling and exloratory economic modelling
  # Only use "IgAEMA" and "IgAEMA plus HLA" and "HLA plus IgATTG" 
  tests_of_interest <- c("IgATTG", "IgAEMA plus HLA", "HLA plus IgATTG")
  
  strategies_of_interest <- c("No screening",paste0(rep(paste0(new_combinations$sens_riskfactor, " ", new_combinations$spec_riskfactor), 2), 
                                                    " ", rep(tests_of_interest, each = dim(new_combinations)[1])))
  strategies_of_interest <- unique(strategies_of_interest)
  names(strategies_of_interest) <- c("No screening", paste(rep(risk_percentages, 2), rep(tests_of_interest, each = length(risk_percentages))))
  
  
  # Table explaining the strategies
  strategy_table <- matrix(NA, nrow = length(strategies_of_interest), ncol = 5)
  rownames(strategy_table) <- strategies_of_interest
  colnames(strategy_table) <- c("Strategy name", "Pre-test probability for blood test",
                                "Test", "Sensitivity", "Specificity")
  strategy_table[, 1] <- names(strategies_of_interest)
  strategy_table[, 2] <- c("", rep(rep(risk_percentages), length(tests_of_interest)))
  strategy_table[, 3] <- c("", rep(tests_of_interest, each = length(risk_percentages)))
  strategy_table[, 4] <- c("", rep(new_combinations$sens_riskfactor, length(tests_of_interest)))
  strategy_table[, 5] <- c("", rep(new_combinations$spec_riskfactor, length(tests_of_interest)))
  
  # Codes and numbers of all risk prediction strategies 
  combinations$x <- paste(combinations$sens_riskfactor, combinations$spec_riskfactor)
  combinations_names <- combinations$x
  n_combinations <- length(combinations$x)
  
  # Define the number and names of tests
  
  n_tests <- (n_sero_tests * n_combinations) + 1
  
  t_names <-  c("No screening", outer(combinations_names, tests, FUN = "paste")[1:n_tests-1])
  
  
  # Define the number and names of states of the model
  n_states <- 7
  state_names <- c("CD GFD no complications",
                   "CD GFD osteoporosis",
                   "CD GFD NHL",
                   "Undiagnosed CD no complications",
                   "Undiagnosed CD osteoporosis",
                   "Undiagnosed CD NHL",
                   "Death")
  
  
  
  #sensitivity analyses
  disutility_fp_diagnosis <- "Yes" #Options are "Yes" or "No". Relates to long term annual disutility of -0.009 for people falsely diagnosed from NICE guideline
  
  

  
  #generate input parameters
  input_parameters <- generate_model_parameters(starting_age, population = population,
                                                combinations = combinations, n_samples = n_samples)
  
  transition_matrices <- generate_transition_matrices(input_parameters, population = population)
  
  
  #generate results
  #system.time({
  #output <- generate_net_benefit(input_parameters, strategies_of_interest, 
  #                               transition_matrices,
  #                               combinations = combinations,
  #                               population = population)
  #})
  #system.time({
  output <- generate_net_benefit_cpp(input_parameters, strategies_of_interest, 
                                 transition_matrices,
                                 combinations = combinations,
                                 population = population)
  #})
  
  # When updating, did a check that the optimised and unoptimised results agreed
  #rowMeans(output$total_costs)[1:10]
  #rowMeans(output$total_qalys)[1:10]
  #rowMeans(output_cpp$total_costs)[1:10]
  #rowMeans(output_cpp$total_qalys)[1:10]
  
  strategies_excluded <- names(subset(output$incremental_net_benefit,output$incremental_net_benefit < 0)) #strategies with ENB less than no screening
  strategies_included <- names(subset(output$incremental_net_benefit,output$incremental_net_benefit > 0)) #strategies with ENB greater than no screening
  
  
  # Simple graphical comparison of biospy proportions
  jpeg(paste0("results/", population, "/", population, " biopsy plot.jpg"))
  par(mfrow=c(1,1))
  plot(output$percentage_biopsy_IgAEMA, col = 0, xlab = "Risk factor strategy", ylab = "Proportion Biospy")
  lines(output$percentage_biopsy_IgAEMA , col = 1)
  lines(output$percentage_biopsy_IgATTGplusEMA, col = 2)
  lines(output$percentage_biopsy_IgATTG, col = 3)
  lines(output$percentage_biopsy_IgAEMAplusHLA, col = 4)
  lines(output$percentage_biopsy_IgATTGplusEMAplusHLA, col = 5)
  lines(output$percentage_biopsy_IgATTGplusHLA, col = 6)
  lines(output$percentage_biopsy_HLAplusIgAEMA, col = 7)
  lines(output$percentage_biopsy_HLAplusIgATTGplusEMA, col = 8)
  lines(output$percentage_biopsy_HLAplusIgATTG, col = 9)
  legend("bottomleft", col = c(1:9), cex = 0.6, lty = 1,
         legend = c("IgAEMA", "IgATTGplusEMA",
                    "IgATTG", "IgAEMAplusHLA",
                    "IgATTGplusEMAplusHLA",
                    "IgATTGplusHLA",
                    "HLAplusIgAEMA", "HLAplusIgATTGplusEMA",
                    "HLAplusIgATTG"))
  
  dev.off()
  
  
  #plots of incremental net benefit
  # Need a simple function to pick out sensitivity, specificity and serological strategy
  inb_strategies <- function(sensitivity_vector, specificity_vector, sero_test) {
    # Identify results for this serological test
    temp <- grepl(sero_test, names(output$incremental_net_benefit))
    # And remove serological tests that might have strings in common
    for(wrong_test in tests[tests != sero_test]) {
      # Only remove if wrong_test is not a substring of sero_test
      if(!grepl(wrong_test, sero_test)) {
        temp <- temp & !grepl(wrong_test, names(output$incremental_net_benefit))  
      }
    }
    # And allow user to pass vector of sensitivity and specificity
    sens_spec_temp <- rep(FALSE, length(output$incremental_net_benefit))
    for(sensitivity in sensitivity_vector) {
      for(specificity in specificity_vector) {
        sens_spec_temp <- sens_spec_temp |
          grepl(paste0(sensitivity, " ", specificity, " "), names(output$incremental_net_benefit))
      }
    }
    # Right serological test and sens/spec in user provided vectors
    temp <- temp & sens_spec_temp
    return(temp)
  }
  
  pdf(paste0("results/", population, "/", population, " inbplot.pdf"))
  par(mfrow=c(3,3))
  par(mar = c(2, 1, 1, 1))
  
  # All tests
  formatted_test_names <- c("IgA EMA", "IgA TTG plus EMA", "IgA TTG", 
                            "IgA EMA plus HLA", "IgA TTG plus EMA plus HLA", "IgA TTG plus HLA",
                            "HLA plus IgA EMA", "HLA plus IgA TTG plus EMA", "HLA plus IgA TTG")
  names(formatted_test_names) <- tests
  # Each of the lines includes six sensitivities
  # One line for each specificity
  for(sero_test in tests) {
    plot(output$incremental_net_benefit[inb_strategies(sensitivity_vector = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.9999), 
                                                       specificity_vector = 0.5, sero_test = sero_test)],
         pch=19, ylim=c(-15000,50000), ylab = "Incremental net benefit", xlab = "sensitivity", main = formatted_test_names[sero_test], xaxt="n" , col = 0, font.main = 1)
    abline(h=0)
    axis(1,                         # Define x-axis manually
         at = 1:6,
         labels = c(0.5, 0.6, 0.7, 0.8, 0.9, 1), font.axis = 1)
    specificity_vector <- c(0.5, 0.6, 0.7, 0.8, 0.9)
    # All specificities
    for(i_specificity in 1:length(specificity_vector)) {
      lines(output$incremental_net_benefit[inb_strategies(sensitivity_vector = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.9999), 
                                                          specificity_vector = specificity_vector[i_specificity], sero_test = sero_test)],
            col = i_specificity, lwd=2, lty= i_specificity)  
      # Lower credible limit
      lines(output$inb_lci[inb_strategies(sensitivity_vector = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.9999), 
                                          specificity_vector = specificity_vector[i_specificity], sero_test = sero_test)],
            col = i_specificity, lwd=1, lty= i_specificity) 
      # Upper credible limit
      lines(output$inb_uci[inb_strategies(sensitivity_vector = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.9999), 
                                          specificity_vector = specificity_vector[i_specificity], sero_test = sero_test)],
            col = i_specificity , lwd=1, lty= i_specificity) 
    } # End loop over specificity
  }
  dev.off()
  
  pdf(paste0("results/", population, "/", population, " legend_inbplot.pdf"))
  par(mfrow=c(1,1))
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend("topleft", legend =c('Specificity 0.5', '0.6', '0.7',
                              '0.8', '0.9'),
         col = c(1:5), lty = (1:5))
  dev.off()
  
  #plots of probability CE
  pdf(paste0("results/", population, "/", population, " probce.pdf"))
  par(mfrow=c(3,3))
  par(mar = c(2, 1, 1, 1))
  
  for(sero_test in tests) {
    plot(output$probability_cost_effective[inb_strategies(sensitivity_vector = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.9999),
                                                          specificity_vector = 0.5,
                                                          sero_test = sero_test)],pch=19, ylim=c(0,1), ylab = "Incremental net benefit", xlab = "sensitivity", main = formatted_test_names[sero_test], xaxt="n" )
    #lines(output$probability_cost_effective[2:7], lwd=2)
    abline(h=0)
    axis(1,                         # Define x-axis manually
         at = 1:6,
         labels = c(0.5, 0.6, 0.7, 0.8, 0.9, 1))
    
    specificity_vector <- c(0.5, 0.6, 0.7, 0.8, 0.9)
    # All specificities
    for(i_specificity in 1:length(specificity_vector)) {
      lines(output$probability_cost_effective[inb_strategies(sensitivity_vector = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.9999),
                                                             specificity_vector = specificity_vector[i_specificity],
                                                             sero_test = sero_test)],
            col = i_specificity, lwd = 2, lty= i_specificity)
    }
  }
  
  dev.off()
  
  write.csv(t(output$total_qalys), paste0("results/", population, "/", population, " total qalys.csv"))
  write.csv(t(output$total_costs), paste0("results/", population, "/", population, " total costs.csv"))
  
  
  
  #table of costs and QALYs
  format_results<-function(x, n_digits = 2)
  {
    paste(round(mean(x, na.rm = TRUE),digits = n_digits), " (",
          round(quantile(x, probs = 0.025, na.rm = TRUE), digits = n_digits), ", ", 
          round(quantile(x, probs = 0.975, na.rm = TRUE), digits = n_digits),")",sep="")
  }
  
  # This where specific strategies are selected
  
  strategies_of_interest_qalys <- output$total_qalys[strategies_of_interest, ]
  strategies_of_interest_costs <- output$total_costs[strategies_of_interest, ]
  strategies_of_interest_inetbenefit <- output$all_incremental_net_benefit[strategies_of_interest, ]
  strategies_of_interest_qalys_table <- array(dim=c(1, length(strategies_of_interest)), dimnames = list(NULL, NULL))
  strategies_of_interest_costs_table <- array(dim=c(1, length(strategies_of_interest)), dimnames = list(NULL, NULL))
  # strategies_of_interest_netbenefit_table <- array(dim=c(1, length(strategies_of_interest)), dimnames = list(NULL, NULL))
  strategies_of_interest_inetbenefit_table <- array(dim=c(1, length(strategies_of_interest)), dimnames = list(NULL, NULL))
  
  for (i in 1:length(strategies_of_interest)) { 
    strategies_of_interest_qalys_table[, i] <- format_results(strategies_of_interest_qalys[i, ], n_digits = 2)
    strategies_of_interest_costs_table[, i] <- format_results(strategies_of_interest_costs[i, ], n_digits = 0)
    # strategies_of_interest_netbenefit_table[, i] <- format_results(strategies_of_interest_netbenefit[i, ], n_digits = 0)
    strategies_of_interest_inetbenefit_table[, i] <- format_results(strategies_of_interest_inetbenefit[i, ], n_digits = 0)
  }
  
  # Create a table comparing costs, QALYs, and incrmental net benefit for strategies of interest
  results_table <- data.frame(c("", rep(risk_percentages, length(tests_of_interest))), 
                      c("", rep(tests_of_interest, each = length(risk_percentages))),
                      c("", rep(new_combinations$sens_riskfactor, length(tests_of_interest))),
                      c("", rep(new_combinations$spec_riskfactor, length(tests_of_interest))),
                      t(strategies_of_interest_costs_table), t(strategies_of_interest_qalys_table), t(strategies_of_interest_inetbenefit_table))
  colnames(results_table) <- c("Pre-test probability for blood test", "Test", 
                       "Sensitivity", "Specificity",
                       "Costs", "QALYs", "Incremental net benefit v no screening")
  # Change the sensitivity and specificity of 0.9999 back to 1 for ease of presentation
  results_table$Sensitivity[results_table$Sensitivity == 0.9999] <- 1
  results_table$Specificity[results_table$Specificity == 0.9999] <- 1
  rownames(results_table) <- strategies_of_interest
  #table <-  table[order( table$`Incremental net benefit v no screening`), ]
  write.csv(results_table, paste0("results/", population, "/", population, " table of results.csv"))
  
  #for CEAC, if doing in excel_ can also do using BCEA
  write.csv(t(output$total_qalys[strategies_of_interest, ]), paste0("results/", population, "/", population, " qalys_feasible.csv"))
  write.csv(t(output$total_costs[strategies_of_interest, ]), paste0("results/", population, "/", population, " costs_feasible.csv"))
  #cost breakdown
  strategies_of_interest_testcosts <- output$test_costs_applied[, strategies_of_interest]
  strategies_of_interest_fpcosts <- output$fp_costs[, strategies_of_interest]
  strategies_of_interest_diagnosiscosts <- output$diagnosis_costs[, strategies_of_interest]
  
  strategies_of_interest_cyclecosts <- output$cycle_costs[strategies_of_interest]
  strategies_of_interest_testcosts_table <- array(dim=c(1, length(strategies_of_interest)), dimnames = list(NULL, NULL))
  strategies_of_interest_fpcosts_table <- array(dim=c(1, length(strategies_of_interest)), dimnames = list(NULL, NULL))
  strategies_of_interest_diagnosiscosts_table <- array(dim=c(1, length(strategies_of_interest)), dimnames = list(NULL, NULL))
  strategies_of_interest_cyclecosts_table <- array(dim=c(1, length(strategies_of_interest)), dimnames = list(NULL, NULL))
  for (i in 1:length(strategies_of_interest)) { 
    strategies_of_interest_testcosts_table[, i] <- format_results(strategies_of_interest_testcosts[, i])
    strategies_of_interest_fpcosts_table[, i] <- format_results(strategies_of_interest_fpcosts[, i])
    strategies_of_interest_diagnosiscosts_table[, i] <- format_results(strategies_of_interest_diagnosiscosts[, i])
    strategies_of_interest_cyclecosts_table[i] <- format_results(strategies_of_interest_cyclecosts[i])
  }
  
  cost_breakdown_table <- data.frame(t(strategies_of_interest_testcosts_table), t(strategies_of_interest_fpcosts_table), t(strategies_of_interest_diagnosiscosts_table), t(strategies_of_interest_cyclecosts_table))
  colnames(cost_breakdown_table) <- c("Test costs", "False positive costs", "Diagnosis costs", "Cycle costs")
  rownames(cost_breakdown_table) <- strategies_of_interest
  
  write.csv(cost_breakdown_table, paste0("results/", population, "/", population, " cost breakdown table.csv"))
  
  
  #utility breakdown
  strategies_of_interest_disutilitybiopsy <- output$disutility_biopsy[, strategies_of_interest]
  strategies_of_interest_disutilitybiopsywait <- output$disutility_biopsy_wait[, strategies_of_interest]
  strategies_of_interest_disutilityfp <- output$disutility_fp[, strategies_of_interest]
  
  strategies_of_interest_cycleqalys <- output$cycle_qalys[strategies_of_interest]
  strategies_of_interest_disutilitybiopsy_table <- array(dim=c(1, length(strategies_of_interest)), dimnames = list(NULL, NULL))
  strategies_of_interest_disutilitybiopsywait_table <- array(dim=c(1, length(strategies_of_interest)), dimnames = list(NULL, NULL))
  strategies_of_interest_disutilityfp_table <- array(dim=c(1, length(strategies_of_interest)), dimnames = list(NULL, NULL))
  strategies_of_interest_cycleqalys_table <- array(dim=c(1, length(strategies_of_interest)), dimnames = list(NULL, NULL))
  for (i in 1:length(strategies_of_interest)) { 
    strategies_of_interest_disutilitybiopsy_table[, i] <- format_results(strategies_of_interest_disutilitybiopsy[, i], n_digits = 4)
    strategies_of_interest_disutilitybiopsywait_table[, i] <- format_results(strategies_of_interest_disutilitybiopsywait[, i], n_digits = 4)
    strategies_of_interest_disutilityfp_table[, i] <- format_results(strategies_of_interest_disutilityfp[, i], n_digits = 4)
    strategies_of_interest_cycleqalys_table[i] <- format_results(strategies_of_interest_cycleqalys[i], n_digits = 4)
  }
  
  qaly_breakdown_table <- data.frame(t(strategies_of_interest_cycleqalys_table), t(strategies_of_interest_disutilityfp_table), t(strategies_of_interest_disutilitybiopsywait_table), t(strategies_of_interest_disutilitybiopsy_table))
  colnames(qaly_breakdown_table) <- c("Cycle QALYs", "Disutility FP", "Disutility waiting for biopsy", "Disutility biopsy")
  rownames(qaly_breakdown_table) <- strategies_of_interest
  
  write.csv(qaly_breakdown_table, paste0("results/", population, "/", population, " qaly breakdown table.csv"))
  
  #time in states
  
  time_in_states <- output$time_in_states
  time_in_states <- cbind(strategy_table, time_in_states)
  
  write.csv(time_in_states, paste0("results/", population, "/", population, " time in states.csv"))
  
  # Now use the BCEA package to analyse the results
  # pkgs <- c("MASS","Rtools","devtools")
  #repos <- c("https://cran.rstudio.com", "https://www.math.ntnu.no/inla/R/stable") 
  #install_packages(pkgs,repos=repos,dependencies = "Depends")
  #devtools::install_github("giabaio/BCEA")
  
  
  #m <- bcea(e = t(output$total_qalys), c = t(output$total_costs), ref = 1, interventions = t_names)
  #summary(m)
  
  m_of_interest <- bcea(e = t(output$total_qalys[strategies_of_interest, ]), 
                        c = t(output$total_costs[strategies_of_interest, ]), ref = 1, 
                        interventions = strategy_table[, "Strategy name"])
  summary(m_of_interest) 
  
  

  mce <- multi.ce(m_of_interest)
 
  pdf(file = paste0("results/", population, "/", population, "_ceac_of_interest.pdf"))
  ceac.plot(mce, graph = c("base"), pos = "topright",
           cex.main = 0.2,
           line = list(colors = c(1:length(strategies_of_interest)))) #CEAC 
  dev.off()
  
  
  
  # Calculate population EVPI

  # Prevalence of CD (different for adults and children)
  cd_prevalence <- 0.01
  # Total population susceptible to CD (i.e. how many adults or children in UK)
  # 0.21 are below the age of 18 https://www.ethnicity-facts-figures.service.gov.uk/uk-population-by-ethnicity/demographics/age-groups/latest#:~:text=by%20ethnicity%20Summary-,The%20data%20shows%20that%3A,aged%2060%20years%20and%20over
  # Total population 67081000 # https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates
  # CD prevalence is 0.01 Bingley 2004
  if(population != "children") total_population <- 66700000 * cd_prevalence * (1 - 0.21) * 0.5
  if(population == "children") total_population <- 66700000 * cd_prevalence * 0.21
  
  
  # Assume test technology remains relevant for at least 10 years
  technology_horizon <- 10
  discounted_population_size <- sum((1/1.035)^(0:(technology_horizon - 1))) * total_population 
  # Use EVPI in only strategies of interest because in practice the others (similar tests or theoretical risk
  # strategies are not real options)
  population_evpi <- m_of_interest$evi[201] * discounted_population_size 
  
  evpi_table <- matrix(NA, nrow = 8, ncol = 2)
  colnames(evpi_table) <- c("Per person", "Population")
  rownames(evpi_table) <- c("total population", "discounted population size", 
                            "evpi",
                            "Utilities and disutilities", "Rates of osteoporosis and NHL",
                            "GFD effect", "Test accuracies", "Probability of late diagnosis")
  # Again the EVPI considering only the strategies of interest (this is less than total)
  evpi_table[1:3, "Per person"] <- c(total_population, discounted_population_size,
                                  m_of_interest$evi[201])
  evpi_table[3, "Population"] <- population_evpi
  write.csv(evpi_table, paste0("results/", population, "/", population, "_evpi_summary.csv"))
  
  
  #info-rank
  # Subset of input parameters actually used for this population subset
  # e.g. Don't need sensitivity and specificity in adults if focusing on children
  if(population == "children")input_parameters_info_rank <- input_parameters[, -which(grepl("_adult", colnames(input_parameters)))]
  if(population != "children")input_parameters_info_rank <- input_parameters[, -which(grepl("_children", colnames(input_parameters)))]
  
  
  # If adults don't need log rates in 0-10 age categories
  
  input_parameters_info_rank <- input_parameters_info_rank[, -which(grepl( "probability_IDA_0", colnames(input_parameters_info_rank)) | grepl("osteoporosis_lograte_0", colnames(input_parameters_info_rank)))]
  # Remove parameters that are fixed at 1 (sens/spec biopsy) or 0 (some costs)
  input_parameters_info_rank <- input_parameters_info_rank[, -which(colSums(input_parameters_info_rank == 1) == n_samples)]
  input_parameters_info_rank <- input_parameters_info_rank[, -which(colSums(input_parameters_info_rank == 0) == n_samples)]
  # Remove probability of no complications as it is functional
  input_parameters_info_rank <- input_parameters_info_rank[, -which(colnames(input_parameters_info_rank) == "probability_nocomplications")]
  print("Conducting info rank analysis")
  
  # Need to plot most influential parameters and say the rest have no impact.
  coeliac_inp <- createInputs(input = input_parameters_info_rank)
  pdf(paste0("results/", population, "/", population, "_info_rank.pdf"))
  info.rank(inp = coeliac_inp, 
            m_of_interest, xlim = c(0,1), wtp = 20000, howManyPars = 10)
  dev.off()
  
  # Run the MLMC for EVPPI analyses
  # This is time consuming (takes about 5 hours per population)
  if(run_mlmc_evppi) {
    source("code/evppi_mlmc_3.R")  
    source("code/evppi_gam_gp_1.R")
  }
} # End loop over populations



