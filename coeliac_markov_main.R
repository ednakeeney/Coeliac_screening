# Coeliac disease Markov model
# Edna Keeney
library(tictoc)
library(BCEA)
library(SimDesign)
library(BCEA)
library(dplyr)

setwd("C:/Users/ek14588/Downloads/Coeliac_screening")

tic()
rm(list=ls())
set.seed(14143)
  
 
  
  tests <- c("IgAEMA", "IgATTGplusEMA", "IgATTG", "IgAEMA plus HLA", "IgATTGplusEMA plus HLA", "IgATTG plus HLA")
  n_tests <- length(tests)
  
  #pre-test probabilities of coeliac disease 
  sens_riskfactor <- c(0.5, 0.6, 0.7, 0.8, 0.9, 0.9999)
  spec_riskfactor <- c(0.5, 0.6, 0.7, 0.8, 0.9, 0.9999)
  combinations <- expand.grid(sens_riskfactor = sens_riskfactor, spec_riskfactor = spec_riskfactor)
  combinations$x <- paste(combinations$sens_riskfactor, combinations$spec_riskfactor)
  combinations_names <- combinations$x
  n_combinations <- length(combinations$x)
  
  # Define the number and names of tests
  
  n_tests <- (n_tests * n_combinations) + 1
  
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
  
  
  # Define simulation parameters
  # This is the number of PSA samples to use
  n_samples <- 1000
 
  perspective <- "NHS" #Options are "NHS" or "NHS+OOP" if out-of-pocket costs for iron supplements and gluten free products are to be included
  
  population <- "children" #Options are "adults" or "children"


  
  #########################################################################################################################
  source("generate_state_costs.R")
  source("generate_state_qalys.R")
  source("generate_model_parameters.R")
  source("generate_transition_matrices.R")
  source("generate_net_benefit_hla-HT.R")
  
  
  #generate input parameters
  input_parameters <- generate_model_parameters(starting_age)
  
  transition_matrices <- generate_transition_matrices(input_parameters)
  
  #generate results
  output <- generate_net_benefit(input_parameters)
  output
  
  strategies_excluded <- names(subset(output$incremental_net_benefit,output$incremental_net_benefit < 0)) #strategies with ENB less than no screening
  strategies_included <- names(subset(output$incremental_net_benefit,output$incremental_net_benefit > 0)) #strategies with ENB less than no screening

  
  output$percentage_biopsy_IgAEMA 
  output$percentage_biopsy_IgATTGplusEMA 
  output$percentage_biopsy_IgATTG 
  output$percentage_biopsy_IgAEMAplusHLA 
  output$percentage_biopsy_IgATTGplusEMAplusHLA 
  output$percentage_biopsy_IgATTGplusHLA 
  
  write.csv(data.frame(output$test_costs, output$diagnosis_costs, output$fp_costs, output$cycle_costs, output$average_costs), "costs.csv")
  
  #ICER_table <- as.data.frame(output$ICER[2:7])
  plot(output$ICER[2:7],pch=19, ylim=c(0,80000), ylab = "ICER", main = "IGA EMA", )
  text(output$ICER[2:7], labels=c("0.5 0.5", "0.6 0.5", "0.7 0.5", "0.8 0.5", "0.9 0.5", "1 0.5"),cex=0.7, font=1, pos=1)
  abline(h=20000)
  abline(h=30000)
  points(output$ICER[8:13], col = 2)
  text(output$ICER[8:13], labels=c("0.5 0.6", "0.6 0.6", "0.7 0.6", "0.8 0.6", "0.9 0.6", "1 0.6"),cex=0.7, font=1, pos=1)
  points(output$ICER[14:19], col = 3)
  text(output$ICER[14:19], labels=c("0.5 0.7", "0.6 0.7", "0.7 0.7", "0.8 0.7", "0.9 0.7", "1 0.7"),cex=0.7, font=1, pos=1)
  points(output$ICER[20:25], col = 4)
  text(output$ICER[20:25], labels=c("0.5 0.8", "0.6 0.8", "0.7 0.8", "0.8 0.8", "0.9 0.8", "1 0.8"),cex=0.7, font=1, pos=1)
  points(output$ICER[26:31], col = 5)
  text(output$ICER[26:31], labels=c("0.5 0.9", "0.6 0.9", "0.7 0.9", "0.8 0.9", "0.9 0.9", "1 0.9"),cex=0.7, font=1, pos=1)
  points(output$ICER[32:37], col = 6)
  text(output$ICER[32:37], labels=c("0.5 1", "0.6 1", "0.7 1", "0.8 1", "0.9 1", "1 1"),cex=0.7, font=1, pos=1)
  
  plot(output$ICER[38:43],pch=19, ylim=c(0,75000), ylab = "ICER", main = "IGA TTG plus EMA", )
  text(output$ICER[38:43], labels=c("0.5 0.5", "0.6 0.5", "0.7 0.5", "0.8 0.5", "0.9 0.5", "1 0.5"),cex=0.7, font=1, pos=1)
  abline(h=20000)
  points(output$ICER[44:49], col = 2)
  text(output$ICER[44:49], labels=c("0.5 0.6", "0.6 0.6", "0.7 0.6", "0.8 0.6", "0.9 0.6", "1 0.6"),cex=0.7, font=1, pos=1)
  points(output$ICER[50:55], col = 3)
  text(output$ICER[50:55], labels=c("0.5 0.7", "0.6 0.7", "0.7 0.7", "0.8 0.7", "0.9 0.7", "1 0.7"),cex=0.7, font=1, pos=1)
  points(output$ICER[56:61], col = 4)
  text(output$ICER[56:61], labels=c("0.5 0.8", "0.6 0.8", "0.7 0.8", "0.8 0.8", "0.9 0.8", "1 0.8"),cex=0.7, font=1, pos=1)
  points(output$ICER[62:67], col = 5)
  text(output$ICER[62:67], labels=c("0.5 0.9", "0.6 0.9", "0.7 0.9", "0.8 0.9", "0.9 0.9", "1 0.9"),cex=0.7, font=1, pos=1)
  points(output$ICER[68:73], col = 6)
  text(output$ICER[68:73], labels=c("0.5 1", "0.6 1", "0.7 1", "0.8 1", "0.9 1", "1 1"),cex=0.7, font=1, pos=1)
  
  plot(output$ICER[74:79],pch=19, ylim=c(0,75000), ylab = "ICER", main = "IGA TTG", )
  text(output$ICER[74:79], labels=c("0.5 0.5", "0.6 0.5", "0.7 0.5", "0.8 0.5", "0.9 0.5", "1 0.5"),cex=0.7, font=1, pos=1)
  abline(h=20000)
  points(output$ICER[80:85], col = 2)
  text(output$ICER[80:85], labels=c("0.5 0.6", "0.6 0.6", "0.7 0.6", "0.8 0.6", "0.9 0.6", "1 0.6"),cex=0.7, font=1, pos=1)
  points(output$ICER[86:91], col = 3)
  text(output$ICER[86:91], labels=c("0.5 0.7", "0.6 0.7", "0.7 0.7", "0.8 0.7", "0.9 0.7", "1 0.7"),cex=0.7, font=1, pos=1)
  points(output$ICER[92:97], col = 4)
  text(output$ICER[92:97], labels=c("0.5 0.8", "0.6 0.8", "0.7 0.8", "0.8 0.8", "0.9 0.8", "1 0.8"),cex=0.7, font=1, pos=1)
  points(output$ICER[98:103], col = 5)
  text(output$ICER[98:103], labels=c("0.5 0.9", "0.6 0.9", "0.7 0.9", "0.8 0.9", "0.9 0.9", "1 0.9"),cex=0.7, font=1, pos=1)
  points(output$ICER[104:109], col = 6)
  text(output$ICER[104:109], labels=c("0.5 1", "0.6 1", "0.7 1", "0.8 1", "0.9 1", "1 1"),cex=0.7, font=1, pos=1)
 
  plot(output$ICER[110:115],pch=19, ylim=c(800, 1100), ylab = "ICER", main = "Tests plus HLA")
  legend(1, 940, legend=c("IgA TTG", "IgA EMA", "IgA TTGplusEMA"),
       pch=c(2, 1,4))
  text(output$ICER[110:115], labels=c("0.5 0.5", "0.6 0.5", "0.7 0.5", "0.8 0.5", "0.9 0.5", "1 0.5"),cex=0.7, font=1, pos=1)
  abline(h=20000)
  points(output$ICER[116:121], col = 2)
  text(output$ICER[116:121], labels=c("0.5 0.6", "0.6 0.6", "0.7 0.6", "0.8 0.6", "0.9 0.6", "1 0.6"),cex=0.7, font=1, pos=1)
  points(output$ICER[122:127], col = 3)
  text(output$ICER[122:127], labels=c("0.5 0.7", "0.6 0.7", "0.7 0.7", "0.8 0.7", "0.9 0.7", "1 0.7"),cex=0.7, font=1, pos=1)
  points(output$ICER[128:133], col = 4)
  text(output$ICER[128:133], labels=c("0.5 0.8", "0.6 0.8", "0.7 0.8", "0.8 0.8", "0.9 0.8", "1 0.8"),cex=0.7, font=1, pos=1)
  points(output$ICER[134:139], col = 5)
  text(output$ICER[134:139], labels=c("0.5 0.9", "0.6 0.9", "0.7 0.9", "0.8 0.9", "0.9 0.9", "1 0.9"),cex=0.7, font=1, pos=1)
  points(output$ICER[140:145], col = 6)
  text(output$ICER[140:145], labels=c("0.5 1", "0.6 1", "0.7 1", "0.8 1", "0.9 1", "1 1"),cex=0.7, font=1, pos=1)
  
  points(output$ICER[146:151],pch=4)
  text(output$ICER[146:151], labels=c("0.5 0.5", "0.6 0.5", "0.7 0.5", "0.8 0.5", "0.9 0.5", "1 0.5"),cex=0.7, font=1, pos=1)
  abline(h=20000)
  points(output$ICER[152:157], pch=4, col = 2)
  text(output$ICER[152:157], labels=c("0.5 0.6", "0.6 0.6", "0.7 0.6", "0.8 0.6", "0.9 0.6", "1 0.6"),cex=0.7, font=1, pos=1)
  points(output$ICER[158:163], pch=4, col = 3)
  text(output$ICER[158:163], labels=c("0.5 0.7", "0.6 0.7", "0.7 0.7", "0.8 0.7", "0.9 0.7", "1 0.7"),cex=0.7, font=1, pos=1)
  points(output$ICER[164:169], pch=4, col = 4)
  text(output$ICER[164:169], labels=c("0.5 0.8", "0.6 0.8", "0.7 0.8", "0.8 0.8", "0.9 0.8", "1 0.8"),cex=0.7, font=1, pos=1)
  points(output$ICER[170:175], pch=4, col = 5)
  text(output$ICER[170:175], labels=c("0.5 0.9", "0.6 0.9", "0.7 0.9", "0.8 0.9", "0.9 0.9", "1 0.9"),cex=0.7, font=1, pos=1)
  points(output$ICER[176:181], pch=4, col = 6)
  text(output$ICER[176:181], labels=c("0.5 1", "0.6 1", "0.7 1", "0.8 1", "0.9 1", "1 1"),cex=0.7, font=1, pos=1)
  
  points(output$ICER[182:187],pch=2)
  text(output$ICER[182:187], labels=c("0.5 0.5", "0.6 0.5", "0.7 0.5", "0.8 0.5", "0.9 0.5", "1 0.5"),cex=0.7, font=1, pos=1)
  abline(h=20000)
  points(output$ICER[188:193], pch=2, col = 2)
  text(output$ICER[188:193], labels=c("0.5 0.6", "0.6 0.6", "0.7 0.6", "0.8 0.6", "0.9 0.6", "1 0.6"),cex=0.7, font=1, pos=1)
  points(output$ICER[194:199], pch=2, col = 3)
  text(output$ICER[194:199], labels=c("0.5 0.7", "0.6 0.7", "0.7 0.7", "0.8 0.7", "0.9 0.7", "1 0.7"),cex=0.7, font=1, pos=1)
  points(output$ICER[200:205], pch=2, col = 4)
  text(output$ICER[200:205], labels=c("0.5 0.8", "0.6 0.8", "0.7 0.8", "0.8 0.8", "0.9 0.8", "1 0.8"),cex=0.7, font=1, pos=1)
  points(output$ICER[206:201], pch=2, col = 5)
  text(output$ICER[206:211], labels=c("0.5 0.9", "0.6 0.9", "0.7 0.9", "0.8 0.9", "0.9 0.9", "1 0.9"),cex=0.7, font=1, pos=1)
  points(output$ICER[212:217], pch=2, col = 6)
  text(output$ICER[212:217], labels=c("0.5 1", "0.6 1", "0.7 1", "0.8 1", "0.9 1", "1 1"),cex=0.7, font=1, pos=1)
  
  #plots of incremental net benefit
  par(mfrow=c(3,2))
  
  #IGA EMA
  plot(output$incremental_net_benefit[2:7],pch=19, ylim=c(-1100,15000), ylab = "Incremental net benefit", xlab = "sensitivity", main = "IGA EMA", xaxt="n" )
 # text(output$incremental_net_benefit[2:7], labels=c("0.5 0.5", "0.6 0.5", "0.7 0.5", "0.8 0.5", "0.9 0.5", "1 0.5"),cex=0.7, font=1, pos=1)
 lines(output$incremental_net_benefit[2:7], lwd=2)
   abline(h=0)
  axis(1,                         # Define x-axis manually
       at = 1:6,
       labels = c(0.5, 0.6, 0.7, 0.8, 0.9, 1))
  lines(output$incremental_net_benefit[8:13], col = 2, lwd=2, lty=2)
  #text(output$incremental_net_benefit[8:13], labels=c("0.5 0.6", "0.6 0.6", "0.7 0.6", "0.8 0.6", "0.9 0.6", "1 0.6"),cex=0.7, font=1, pos=1)
  lines(output$incremental_net_benefit[14:19], col = 3, lwd=4, lty=3)
 # text(output$incremental_net_benefit[14:19], labels=c("0.5 0.7", "0.6 0.7", "0.7 0.7", "0.8 0.7", "0.9 0.7", "1 0.7"),cex=0.7, font=1, pos=1)
  lines(output$incremental_net_benefit[20:25], col = 4, lwd=2, lty=4)
  #text(output$incremental_net_benefit[20:25], labels=c("0.5 0.8", "0.6 0.8", "0.7 0.8", "0.8 0.8", "0.9 0.8", "1 0.8"),cex=0.7, font=1, pos=1)
  lines(output$incremental_net_benefit[26:31], col = 5, lwd=2, lty=5)
 # text(output$incremental_net_benefit[26:31], labels=c("0.5 0.9", "0.6 0.9", "0.7 0.9", "0.8 0.9", "0.9 0.9", "1 0.9"),cex=0.7, font=1, pos=1)
 # lines(output$incremental_net_benefit[32:37], col = 6, lwd=2, lty=6)
  #text(output$incremental_net_benefit[32:37], labels=c("0.5 1", "0.6 1", "0.7 1", "0.8 1", "0.9 1", "1 1"),cex=0.7, font=1, pos=1)
  legend(0.9, 15000, legend=c("Specificity 0.5", "0.6", "0.7", "0.8", "0.9"),
           lty=c(1:5), col=c(1:5))
  
  
  #IGA TTG plus EMA
  plot(output$incremental_net_benefit[38:43],pch=19, ylim=c(-1100,15000), ylab = "incremental_net_benefit", main = "IGA TTG plus EMA", )
  text(output$incremental_net_benefit[38:43], labels=c("0.5 0.5", "0.6 0.5", "0.7 0.5", "0.8 0.5", "0.9 0.5", "1 0.5"),cex=0.7, font=1, pos=1)
  abline(h=0)
  points(output$incremental_net_benefit[44:49], col = 2)
  text(output$incremental_net_benefit[44:49], labels=c("0.5 0.6", "0.6 0.6", "0.7 0.6", "0.8 0.6", "0.9 0.6", "1 0.6"),cex=0.7, font=1, pos=1)
  points(output$incremental_net_benefit[50:55], col = 3)
  text(output$incremental_net_benefit[50:55], labels=c("0.5 0.7", "0.6 0.7", "0.7 0.7", "0.8 0.7", "0.9 0.7", "1 0.7"),cex=0.7, font=1, pos=1)
  points(output$incremental_net_benefit[56:61], col = 4)
  text(output$incremental_net_benefit[56:61], labels=c("0.5 0.8", "0.6 0.8", "0.7 0.8", "0.8 0.8", "0.9 0.8", "1 0.8"),cex=0.7, font=1, pos=1)
  points(output$incremental_net_benefit[62:67], col = 5)
  text(output$incremental_net_benefit[62:67], labels=c("0.5 0.9", "0.6 0.9", "0.7 0.9", "0.8 0.9", "0.9 0.9", "1 0.9"),cex=0.7, font=1, pos=1)
  points(output$incremental_net_benefit[68:73], col = 6)
  text(output$incremental_net_benefit[68:73], labels=c("0.5 1", "0.6 1", "0.7 1", "0.8 1", "0.9 1", "1 1"),cex=0.7, font=1, pos=1)
  
  #IGA TTG
  plot(output$incremental_net_benefit[74:79],pch=19, ylim=c(-1100,15000), ylab = "incremental_net_benefit", main = "IGA TTG", )
  text(output$incremental_net_benefit[74:79], labels=c("0.5 0.5", "0.6 0.5", "0.7 0.5", "0.8 0.5", "0.9 0.5", "1 0.5"),cex=0.7, font=1, pos=1)
  abline(h=0)
  points(output$incremental_net_benefit[80:85], col = 2)
  text(output$incremental_net_benefit[80:85], labels=c("0.5 0.6", "0.6 0.6", "0.7 0.6", "0.8 0.6", "0.9 0.6", "1 0.6"),cex=0.7, font=1, pos=1)
  points(output$incremental_net_benefit[86:91], col = 3)
  text(output$incremental_net_benefit[86:91], labels=c("0.5 0.7", "0.6 0.7", "0.7 0.7", "0.8 0.7", "0.9 0.7", "1 0.7"),cex=0.7, font=1, pos=1)
  points(output$incremental_net_benefit[92:97], col = 4)
  text(output$incremental_net_benefit[92:97], labels=c("0.5 0.8", "0.6 0.8", "0.7 0.8", "0.8 0.8", "0.9 0.8", "1 0.8"),cex=0.7, font=1, pos=1)
  points(output$incremental_net_benefit[98:103], col = 5)
  text(output$incremental_net_benefit[98:103], labels=c("0.5 0.9", "0.6 0.9", "0.7 0.9", "0.8 0.9", "0.9 0.9", "1 0.9"),cex=0.7, font=1, pos=1)
  points(output$incremental_net_benefit[104:109], col = 6)
  text(output$incremental_net_benefit[104:109], labels=c("0.5 1", "0.6 1", "0.7 1", "0.8 1", "0.9 1", "1 1"),cex=0.7, font=1, pos=1)
  
 
   plot(output$incremental_net_benefit[110:115],pch=19, ylim=c(-1100, 15000),ylab = "incremental_net_benefit", main = "IgA EMA plus HLA")
  #legend(1, 940, legend=c("IgA TTG", "IgA EMA", "IgA TTGplusEMA"),
      #   pch=c(2, 1,4))
   abline(h=0)
  text(output$incremental_net_benefit[110:115], labels=c("0.5 0.5", "0.6 0.5", "0.7 0.5", "0.8 0.5", "0.9 0.5", "1 0.5"),cex=0.7, font=1, pos=1)
  points(output$incremental_net_benefit[116:121], col = 2)
  text(output$incremental_net_benefit[116:121], labels=c("0.5 0.6", "0.6 0.6", "0.7 0.6", "0.8 0.6", "0.9 0.6", "1 0.6"),cex=0.7, font=1, pos=1)
  points(output$incremental_net_benefit[122:127], col = 3)
  text(output$incremental_net_benefit[122:127], labels=c("0.5 0.7", "0.6 0.7", "0.7 0.7", "0.8 0.7", "0.9 0.7", "1 0.7"),cex=0.7, font=1, pos=1)
  points(output$incremental_net_benefit[128:133], col = 4)
  text(output$incremental_net_benefit[128:133], labels=c("0.5 0.8", "0.6 0.8", "0.7 0.8", "0.8 0.8", "0.9 0.8", "1 0.8"),cex=0.7, font=1, pos=1)
  points(output$incremental_net_benefit[134:139], col = 5)
  text(output$incremental_net_benefit[134:139], labels=c("0.5 0.9", "0.6 0.9", "0.7 0.9", "0.8 0.9", "0.9 0.9", "1 0.9"),cex=0.7, font=1, pos=1)
  points(output$incremental_net_benefit[140:145], col = 6)
  text(output$incremental_net_benefit[140:145], labels=c("0.5 1", "0.6 1", "0.7 1", "0.8 1", "0.9 1", "1 1"),cex=0.7, font=1, pos=1)
  
  plot(output$incremental_net_benefit[146:151],pch=4, ylim=c(-1100, 15000),ylab = "incremental_net_benefit", main = "IgATTG plus EMA plus HLA")
  text(output$incremental_net_benefit[146:151], labels=c("0.5 0.5", "0.6 0.5", "0.7 0.5", "0.8 0.5", "0.9 0.5", "1 0.5"),cex=0.7, font=1, pos=1)
  abline(h=20000)
  points(output$incremental_net_benefit[152:157], pch=4, col = 2)
  text(output$incremental_net_benefit[152:157], labels=c("0.5 0.6", "0.6 0.6", "0.7 0.6", "0.8 0.6", "0.9 0.6", "1 0.6"),cex=0.7, font=1, pos=1)
  points(output$incremental_net_benefit[158:163], pch=4, col = 3)
  text(output$incremental_net_benefit[158:163], labels=c("0.5 0.7", "0.6 0.7", "0.7 0.7", "0.8 0.7", "0.9 0.7", "1 0.7"),cex=0.7, font=1, pos=1)
  points(output$incremental_net_benefit[164:169], pch=4, col = 4)
  text(output$incremental_net_benefit[164:169], labels=c("0.5 0.8", "0.6 0.8", "0.7 0.8", "0.8 0.8", "0.9 0.8", "1 0.8"),cex=0.7, font=1, pos=1)
  points(output$incremental_net_benefit[170:175], pch=4, col = 5)
  text(output$incremental_net_benefit[170:175], labels=c("0.5 0.9", "0.6 0.9", "0.7 0.9", "0.8 0.9", "0.9 0.9", "1 0.9"),cex=0.7, font=1, pos=1)
  points(output$incremental_net_benefit[176:181], pch=4, col = 6)
  text(output$incremental_net_benefit[176:181], labels=c("0.5 1", "0.6 1", "0.7 1", "0.8 1", "0.9 1", "1 1"),cex=0.7, font=1, pos=1)
  
  plot(output$incremental_net_benefit[182:187],pch=2, ylim=c(-1100, 15000),ylab = "incremental_net_benefit", main = "IgA TTG plus HLA")
 abline(h=0)
   text(output$incremental_net_benefit[182:187], labels=c("0.5 0.5", "0.6 0.5", "0.7 0.5", "0.8 0.5", "0.9 0.5", "1 0.5"),cex=0.7, font=1, pos=1)
  points(output$incremental_net_benefit[188:193], pch=2, col = 2)
  text(output$incremental_net_benefit[188:193], labels=c("0.5 0.6", "0.6 0.6", "0.7 0.6", "0.8 0.6", "0.9 0.6", "1 0.6"),cex=0.7, font=1, pos=1)
  points(output$incremental_net_benefit[194:199], pch=2, col = 3)
  text(output$incremental_net_benefit[194:199], labels=c("0.5 0.7", "0.6 0.7", "0.7 0.7", "0.8 0.7", "0.9 0.7", "1 0.7"),cex=0.7, font=1, pos=1)
  points(output$incremental_net_benefit[200:205], pch=2, col = 4)
  text(output$incremental_net_benefit[200:205], labels=c("0.5 0.8", "0.6 0.8", "0.7 0.8", "0.8 0.8", "0.9 0.8", "1 0.8"),cex=0.7, font=1, pos=1)
  points(output$incremental_net_benefit[206:211], pch=2, col = 5)
  text(output$incremental_net_benefit[206:211], labels=c("0.5 0.9", "0.6 0.9", "0.7 0.9", "0.8 0.9", "0.9 0.9", "1 0.9"),cex=0.7, font=1, pos=1)
  points(output$incremental_net_benefit[212:217], pch=2, col = 6)
  text(output$incremental_net_benefit[212:217], labels=c("0.5 1", "0.6 1", "0.7 1", "0.8 1", "0.9 1", "1 1"),cex=0.7, font=1, pos=1)
  
  
  write.csv(t(output$total_costs), "costs.csv")
  write.csv(t(output$total_qalys), "qalys.csv")
  
  write.csv(t(output$total_qalys), "qalys.csv")
  
  results <- data.frame(output$average_costs, output$average_effects)
  results <- results[order(results$output.average_costs),]
  for (i in 1:216) {
  results$x[i+1] <- ifelse(results$output.average_effects[i+1] > results$output.average_effects[i], 1,0) }
  results <- subset(results, x== 1)
  for (i in 1:138) {
    results$x[i+1] <- ifelse(results$output.average_effects[i+1] > results$output.average_effects[i], 1,0) }
  results <- subset(results, x== 1)
  for (i in 1:101) {
    results$x[i+1] <- ifelse(results$output.average_effects[i+1] > results$output.average_effects[i], 1,0) }
  results <- subset(results, x== 1)
  nrow(results)
  for (i in 1:77) {
    results$x[i+1] <- ifelse(results$output.average_effects[i+1] > results$output.average_effects[i], 1,0) }
  results <- subset(results, x== 1)
  nrow(results)
  for (i in 1:58) {
    results$x[i+1] <- ifelse(results$output.average_effects[i+1] > results$output.average_effects[i], 1,0) }
  results <- subset(results, x== 1)
  nrow(results)
  for (i in 1:41) {
    results$x[i+1] <- ifelse(results$output.average_effects[i+1] > results$output.average_effects[i], 1,0) }
  results <- subset(results, x== 1)
  nrow(results)
  for (i in 1:30) {
    results$x[i+1] <- ifelse(results$output.average_effects[i+1] > results$output.average_effects[i], 1,0) }
  results <- subset(results, x== 1)
  nrow(results)
  for (i in 1:28) {
    results$x[i+1] <- ifelse(results$output.average_effects[i+1] > results$output.average_effects[i], 1,0) }
  results <- subset(results, x== 1)
  nrow(results)
  for (i in 1:27) {
    results$x[i+1] <- ifelse(results$output.average_effects[i+1] > results$output.average_effects[i], 1,0) }
  results <- subset(results, x== 1)
  nrow(results)
  
  results$ICER <- 0
  for (i in 1:26) {
    results$cost_difference[i+1] <- results$output.average_costs[i+1] - results$output.average_costs[i]
   results$effect_difference[i+1] <- results$output.average_effects[i+1] - results$output.average_effects[i]
   # results$ICER[i+1] <- (results$output.average_costs[i+1] - results$output.average_costs[i])/ (results$output.average_effects[i+1] - results$output.average_effects[i]) }
  }
  # Now use the BCEA package to analyse the results
 # pkgs <- c("MASS","Rtools","devtools")
  #repos <- c("https://cran.rstudio.com", "https://www.math.ntnu.no/inla/R/stable") 
  #install.packages(pkgs,repos=repos,dependencies = "Depends")
  #devtools::install_github("giabaio/BCEA")
  
m <- bcea(e = t(output$total_qalys[1:37,]), c = t(output$total_costs[1:37,]), ref = 1, interventions = t_names[1:37])
summary(m)
m <- bcea(e = t(output$total_qalys[38:73,]), c = t(output$total_costs[38:73,]), ref = 1, interventions = t_names[38:73])
summary(m)
m <- bcea(e = t(output$total_qalys[74:109,]), c = t(output$total_costs[74:109,]), ref = 1, interventions = t_names[74:109])
summary(m)
m <- bcea(e = t(output$total_qalys[110:145,]), c = t(output$total_costs[110:145,]), ref = 1, interventions = t_names[110:145])
summary(m)

m <- bcea(e = t(output$total_qalys), c = t(output$total_costs), ref = 1, interventions = t_names)
summary(m)

eib.plot(m, comparison = NULL, pos =
           c(1, 0), size = NULL, plot.cri = NULL, graph
         = c("ggplot2"))
evi.plot(m, graph = c("base", "ggplot2",
                       "plotly"))
ceac.plot(m, comparison = NULL,
          pos = FALSE, graph = c("ggplot2"))

mce <- multi.ce(m)
ceaf.plot(mce, graph = c("ggplot2"))


ceplane.plot(m, comparison =
               NULL, pos = c(1, 0), graph = c("ggplot2"), point_colors = c(1:36))
sim.table(m)
toc()

