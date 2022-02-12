# Coeliac disease Markov model
# Edna Keeney

# Notes to Edna:
# Left disutility of NHL as an annual decrement as it represents maintenance therapy
# TODO: Update subfertility rates to be female only.
# TODO: Add subdirectories to results for each subpopulation and the sensitivity
# TODO: Add INB to time_in_states (line ~524)
# TODO: info.rank needs to plot only top parameters (line ~584) as most have zero impact
# and only include sens/spec for feasible strategies
# TODO: Check what to do about percentage_male as subfertility should not affect men

library(tictoc)
library(BCEA)
library(SimDesign)
library(BCEA)
library(dplyr)
library(EnvStats)
library(Rmisc)
library(ggplot2)

#setwd("C:/Users/ek14588/Downloads/Coeliac_screening")

tic()
rm(list=ls())
set.seed(14143)
  
# Define simulation parameters
# This is the number of PSA samples to use
n_samples <- 1000

perspective <- "NHS" #Options are "NHS" or "NHS+OOP" if out-of-pocket costs for iron supplements and gluten free products are to be included

population <- "adults" #Options are "adults" or "children"

 
  
  tests <- c("IgAEMA", "IgATTGplusEMA", "IgATTG", "IgAEMA plus HLA", "IgATTGplusEMA plus HLA", "IgATTG plus HLA")
  n_sero_tests <- length(tests)
  
  #pre-test probabilities of coeliac disease 
  sens_riskfactor <- c(0.5, 0.6, 0.7, 0.8, 0.9, 0.9999)
  spec_riskfactor <- c(0.5, 0.6, 0.7, 0.8, 0.9, 0.9999)
  combinations <- expand.grid(sens_riskfactor = sens_riskfactor, spec_riskfactor = spec_riskfactor)
  # Add specific combinations from the risk prediction tool
  included_strategies <- as.data.frame(read.csv("data/included_strategies.csv"))
  if(population == "adults") {
    new_combinations <- as.data.frame(matrix(nrow = 10, ncol = 2))
    colnames(new_combinations) <- c("sens_riskfactor", "spec_riskfactor")
    new_combinations$sens_riskfactor <- included_strategies[included_strategies$Population == "Women" | included_strategies$Population == "Men" , "Sens"]
    new_combinations$spec_riskfactor <- included_strategies[included_strategies$Population == "Women" | included_strategies$Population == "Men" , "Spec"]
    new_combinations$sens_riskfactor[new_combinations$sens_riskfactor == 1] <- 0.9999
    new_combinations$spec_riskfactor[new_combinations$spec_riskfactor == 1] <- 0.9999
    new_combinations <- unique(new_combinations)
    combinations <- rbind(combinations, new_combinations)
    

  } else {
    new_combinations <- as.data.frame(matrix(nrow = 5, ncol = 2))
    colnames(new_combinations) <- c("sens_riskfactor", "spec_riskfactor")
    new_combinations$sens_riskfactor <- included_strategies[included_strategies$Population == "Children" , "Sens"]
    new_combinations$spec_riskfactor <- included_strategies[included_strategies$Population == "Children", "Spec"]
    new_combinations <- unique(new_combinations)
    combinations <- rbind(combinations, new_combinations)
  }
  # Specify the strategies of interest based on prediction modelling and exloratory economic modelling
  # Only use "IgAEMA" and "IgAEMA plus HLA"
  strategies_of_interest <- c("No screening", paste0(rep(new_combinations$sens_riskfactor, each = 2), " ", 
                                                     rep(new_combinations$spec_riskfactor, each = 2), " ",
                                                     rep(c("IgAEMA", "IgAEMA plus HLA"), dim(new_combinations)[1])))
  strategies_of_interest <- unique(strategies_of_interest)
  
  combinations$x <- paste(combinations$sens_riskfactor, combinations$spec_riskfactor)
  combinations_names <- combinations$x
  n_combinations <- length(combinations$x)
  
  # Define the number and names of tests
  
  n_tests <- (n_sero_tests * n_combinations) + 1
  
  t_names <-  c("No screening", outer(combinations_names, tests, FUN = "paste")[1:n_tests-1])

 
  # Define the number and names of states of the model
  n_states <- 9
  state_names <- c("CD GFD no complications",
                  "CD GFD subfertility",
                  "CD GFD osteoporosis",
                  "CD GFD NHL",
                  "Undiagnosed CD no complications",
                  "Undiagnosed CD subfertility",
                  "Undiagnosed CD osteoporosis",
                  "Undiagnosed CD NHL",
                  "Death")
  
  
  
  #sensitivity analyses
disutility_fp_diagnosis <- "Yes" #Options are "Yes" or "No". Relates to long term annual disutility of -0.009 for people falsely diagnosed from NICE guideline

  
  #########################################################################################################################
  source("code/generate_state_costs.R")
  source("code/generate_state_qalys.R")
  source("code/generate_model_parameters.R")
  source("code/generate_transition_matrices.R")
  source("code/generate_net_benefit_hla-HT.R")
  
  
  #generate input parameters
  input_parameters <- generate_model_parameters(starting_age)
  
 transition_matrices <- generate_transition_matrices(input_parameters)
  
  #generate results
  output <- generate_net_benefit(input_parameters, strategies_of_interest)

  
  
  strategies_excluded <- names(subset(output$incremental_net_benefit,output$incremental_net_benefit < 0)) #strategies with ENB less than no screening
  strategies_included <- names(subset(output$incremental_net_benefit,output$incremental_net_benefit > 0)) #strategies with ENB greater than no screening

  
  
  # Simple graphical comparison of biospy proportions
  jpeg("results/ biopsy plot.jpeg")
  par(mfrow=c(1,1))
  plot(output$percentage_biopsy_IgAEMA, col = 0, xlab = "Risk factor strategy", ylab = "Proportion Biospy")
  lines(output$percentage_biopsy_IgAEMA , col = 1)
  lines(output$percentage_biopsy_IgATTGplusEMA, col = 2)
  lines(output$percentage_biopsy_IgATTG, col = 3)
  lines(output$percentage_biopsy_IgAEMAplusHLA, col = 4)
  lines(output$percentage_biopsy_IgATTGplusEMAplusHLA, col = 5)
  lines(output$percentage_biopsy_IgATTGplusHLA, col = 6)
  legend("bottomleft", col = c(1:6), cex = 0.6, lty = 1,
         legend = c("IgAEMA", "IgATTGplusEMA",
                    "IgATTG", "IgAEMAplusHLA",
                    "IgATTGplusEMAplusHLA",
                    "IgATTGplusHLA"))
  
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
  
 jpeg("results/inbplot.jpeg")
  par(mfrow=c(3,2))
  par(mar = c(2, 1, 1, 1))
  
  # All tests
  # Each of the lines includes six sensitivities
  # One line for each specificity
  for(sero_test in tests) {
    plot(output$incremental_net_benefit[inb_strategies(sensitivity_vector = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.9999), 
                                                       specificity_vector = 0.5, sero_test = sero_test)],
         pch=19, ylim=c(-15000,50000), ylab = "Incremental net benefit", xlab = "sensitivity", main = sero_test, xaxt="n" , col = 0)
    abline(h=0)
    axis(1,                         # Define x-axis manually
         at = 1:6,
         labels = c(0.5, 0.6, 0.7, 0.8, 0.9, 1))
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

jpeg("results/legend_inbplot.jpeg")
par(mfrow=c(1,1))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c('Specificity 0.5', '0.6', '0.7',
                            '0.8', '0.9'),
       col = c(1:5), lty = (1:5))
dev.off()

#plots of probability CE
jpeg("results/probce.jpeg")
par(mfrow=c(3,2))
par(mar = c(2, 1, 1, 1))

for(sero_test in tests) {
  plot(output$probability_cost_effective[inb_strategies(sensitivity_vector = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.9999),
                                                        specificity_vector = 0.5,
                                                        sero_test = sero_test)],pch=19, ylim=c(0,1), ylab = "Incremental net benefit", xlab = "sensitivity", main = "IGA EMA", xaxt="n" )
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

write.csv(t(output$total_qalys), "results/total qalys.csv")
write.csv(t(output$total_costs), "results/total costs.csv")



#table of costs and QALYs
format_results<-function(x, n_digits = 2)
{
  paste(round(mean(x),digits = n_digits), " (",
        round(quantile(x, probs = 0.025), digits = n_digits), ", ", 
        round(quantile(x, probs = 0.975), digits = n_digits),")",sep="")
}

# This where specific strategies are selected
# Need to update to use results of risk prediction model
# Define a feasible_strategies object used to select consistently

strategies_of_interest_qalys <- output$total_qalys[strategies_of_interest, ]
strategies_of_interest_costs <- output$total_costs[strategies_of_interest, ]
strategies_of_interest_inetbenefit <- output$all_incremental_net_benefit[strategies_of_interest, ]
strategies_of_interest_qalys_table <- array(dim=c(1, length(strategies_of_interest)), dimnames = list(NULL, NULL))
strategies_of_interest_costs_table <- array(dim=c(1, length(strategies_of_interest)), dimnames = list(NULL, NULL))
# strategies_of_interest_netbenefit_table <- array(dim=c(1, length(strategies_of_interest)), dimnames = list(NULL, NULL))
strategies_of_interest_inetbenefit_table <- array(dim=c(1, length(strategies_of_interest)), dimnames = list(NULL, NULL))

for (i in 1:length(strategies_of_interest)) { 
  strategies_of_interest_qalys_table[, i] <- format_results(strategies_of_interest_qalys[i, ], n_digits = 4)
  strategies_of_interest_costs_table[, i] <- format_results(strategies_of_interest_costs[i, ], n_digits = 0)
  # strategies_of_interest_netbenefit_table[, i] <- format_results(strategies_of_interest_netbenefit[i, ], n_digits = 0)
  strategies_of_interest_inetbenefit_table[, i] <- format_results(strategies_of_interest_inetbenefit[i, ], n_digits = 0)
}

table <- data.frame(t(strategies_of_interest_costs_table), t(strategies_of_interest_qalys_table), t(strategies_of_interest_inetbenefit_table))
colnames(table) <- c("Costs", "QALYs", "Incremental net benefit v no screening")
rownames(table) <- strategies_of_interest
table <-  table[order( table$`Incremental net benefit v no screening`), ]
write.csv(table, "results/table of results.csv")

#for CEAC, if doing in excel_ can also do using BCEA
write.csv(t(output$total_qalys[strategies_of_interest, ]), "qalys_feasible.csv")
write.csv(t(output$total_costs[strategies_of_interest, ]), "costs_feasible.csv")
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

write.csv(cost_breakdown_table, "results/cost breakdown table.csv")


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

write.csv(qaly_breakdown_table, "results/qaly breakdown table.csv")

#time in states

time_in_states <- output$time_in_states

write.csv(output$time_in_states, "results/time in states.csv")

# Now use the BCEA package to analyse the results
# pkgs <- c("MASS","Rtools","devtools")
#repos <- c("https://cran.rstudio.com", "https://www.math.ntnu.no/inla/R/stable") 
#install_packages(pkgs,repos=repos,dependencies = "Depends")
#devtools::install_github("giabaio/BCEA")


m <- bcea(e = t(output$total_qalys), c = t(output$total_costs), ref = 1, interventions = t_names)
summary(m)

m_of_interest <- bcea(e = t(output$total_qalys[strategies_of_interest, ]), 
                      c = t(output$total_costs[strategies_of_interest, ]), ref = 1, 
                      interventions = strategies_of_interest)
summary(m_of_interest) 


eib.plot(m_of_interest, comparison = NULL, pos =
           c(1, 0), size = NULL, plot.cri = NULL, graph
         = c("ggplot2"))
evi.plot(m_of_interest, graph = c("base", "ggplot2",
                       "plotly"))
ceac.plot(m_of_interest, comparison = NULL,
          pos = FALSE, graph = c("ggplot2"))
par(mfrow = c(1,1))
mce <- multi.ce(m_of_interest)
ceaf.plot(mce, graph = c("ggplot2"))
mce.plot(mce, color = c(1:13)) #CEAC 

ceplane.plot(m_of_interest, comparison =
               NULL, pos = c(1, 0), graph = c("ggplot2"), point_colors = c(1:13))
sim.table(m_of_interest)
toc()

# Calculate population EVPI
m <- bcea(e = t(output$total_qalys), c = t(output$total_costs), ref = 1, interventions = t_names, wtp = 20000)
m$evi
# Prevalence of CD (different for adults and children)
cd_prevalence <- 333767
# Total population susceptible to CD (i.e. how many adults or children in UK)
total_population <- 1200000
# Assume test technology remains relevant for at least 10 years
technology_horizon <- 10
discounted_population_size <- sum((1/1.035)^(0:(technology_horizon - 1))) * total_population * cd_prevalence
population_evpi <- m$evi * discounted_population_size 

#info-rank
# Need to plot most influential parameters and say the rest have no impact.
# input_parameters needs to include sens and spec only relevant to feasible_strategies
utility_parameters <- data.frame(input_parameters$utility_GFD, input_parameters$utility_undiagnosedCD, input_parameters$disutility_subfertility, 
                                 input_parameters$disutility_osteoporosis, input_parameters$disutility_NHL, input_parameters$disutility_biopsy, 
                                 input_parameters$disutility_biopsy_wait, input_parameters$disutility_fp)
cost_parameters <- data.frame( input_parameters$cost_CDGFD, input_parameters$cost_osteoporosis, input_parameters$cost_undiagnosedCD, 
                               input_parameters$cost_IDA, input_parameters$cost_biopsy, 
                               input_parameters$cost_subfertility, input_parameters$cost_NHL,input_parameters$cost_diagnosis,
                               input_parameters$test_cost_IgAEMA, input_parameters$test_cost_IgATTG, input_parameters$test_cost_HLA)
info.rank(parameter = colnames(input_parameters), input = input_parameters, m_of_interest, xlim = c(0,0.005), wtp=30000)
info.rank(parameter = colnames(utility_parameters), input = utility_parameters, m_of_interest, xlim = c(0,0.0005), wtp=30000)
info.rank(parameter = colnames(cost_parameters), input = cost_parameters, m_of_interest, xlim = c(0,0.00005), wtp=30000)




