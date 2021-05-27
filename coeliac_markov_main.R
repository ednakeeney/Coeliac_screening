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
  
 
  
  treatments <- c("IgAEMA", "IgATTGplusEMA", "IgATTG", "IgAEMA plus HLA", "IgATTGplusEMA plus HLA", "IgATTG plus HLA")
  n_tests <- length(treatments)
  
  #pre-test probabilities of coeliac disease 
  sens_riskfactor <- c(0.5, 0.6, 0.7, 0.8, 0.9, 0.9999)
  spec_riskfactor <- c(0.5, 0.6, 0.7, 0.8, 0.9, 0.9999)
  combinations <- expand.grid(sens_riskfactor = sens_riskfactor, spec_riskfactor = spec_riskfactor)
  combinations$x <- paste(combinations$sens_riskfactor, combinations$spec_riskfactor)
  combinations_names <- combinations$x
  n_combinations <- length(combinations$x)
  
  # Define the number and names of treatments
  
  n_treatments <- (n_tests * n_combinations) + 1
  
  t_names <-  c("No screening", outer(combinations_names, treatments, FUN = "paste")[1:n_treatments-1])

  
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
  source("generate_net_benefit_hla.R")
  
  
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
  
  ICER_table <- as.data.frame(output$ICER[2:7])
  plot(output$ICER[2:7],pch=19, ylim=c(0,75000), ylab = "ICER", main = "IGA EMA")
  text(output$ICER[2:7], labels=c("0.5 0.5", "0.6 0.5", "0.7 0.5", "0.8 0.5", "0.9 0.5", "1 0.5"),cex=0.7, font=1, pos=3)
  abline(h=20000)
  points(output$ICER[8:13])
  text(output$ICER[8:13], labels=c("0.5 0.6", "0.6 0.6", "0.7 0.6", "0.8 0.6", "0.9 0.6", "1 0.6"),cex=0.7, font=1, pos=1)
  points(output$ICER[14:19], col = 2)
  text(output$ICER[14:19], labels=c("0.5 0.7", "0.6 0.7", "0.7 0.7", "0.8 0.7", "0.9 0.7", "1 0.7"),cex=0.7, font=1, pos=1)
  
  
  write.csv(t(output$total_costs), "costs.csv")
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
               NULL, pos = c(1, 0), graph = c("ggplot2"), point_colors = c(1:14))
sim.table(m)
toc()

