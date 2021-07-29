# Coeliac disease Markov model
# Edna Keeney
library(tictoc)
library(BCEA)
library(SimDesign)
library(BCEA)
library(dplyr)
library(EnvStats)
library(Rmisc)

setwd("C:/Users/ek14588/Downloads/Coeliac_screening")

tic()
rm(list=ls())
set.seed(14143)
  
 
  
  tests <- c("IgAEMA", "IgATTGplusEMA", "IgATTG", "IgAEMA plus HLA", "IgATTGplusEMA plus HLA", "IgATTG plus HLA")
  n_sero_tests <- length(tests)
  
  #pre-test probabilities of coeliac disease 
  sens_riskfactor <- c(0.5, 0.6, 0.7, 0.8, 0.99, 0.9999)
  spec_riskfactor <- c(0.5, 0.6, 0.7, 0.8, 0.99, 0.9999)
  combinations <- expand.grid(sens_riskfactor = sens_riskfactor, spec_riskfactor = spec_riskfactor)
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
  
  
  # Define simulation parameters
  # This is the number of PSA samples to use
  n_samples <- 1000
 
  perspective <- "NHS" #Options are "NHS" or "NHS+OOP" if out-of-pocket costs for iron supplements and gluten free products are to be included
  
  population <- "adults" #Options are "adults" or "children"
  
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
  output <- generate_net_benefit(input_parameters)

  
  
  strategies_excluded <- names(subset(output$incremental_net_benefit,output$incremental_net_benefit < 0)) #strategies with ENB less than no screening
  strategies_included <- names(subset(output$incremental_net_benefit,output$incremental_net_benefit > 0)) #strategies with ENB greater than no screening

  
  output$percentage_biopsy_IgAEMA 
  output$percentage_biopsy_IgATTGplusEMA 
  output$percentage_biopsy_IgATTG 
  output$percentage_biopsy_IgAEMAplusHLA 
  output$percentage_biopsy_IgATTGplusEMAplusHLA 
  output$percentage_biopsy_IgATTGplusHLA 
  
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
 jpeg("results/inbplot.jpeg")
  par(mfrow=c(3,2))
  par(mar = c(2, 1, 1, 1))
  
  #IGA EMA

  plot(output$incremental_net_benefit[2:7],pch=19, ylim=c(-15000,50000), ylab = "Incremental net benefit", xlab = "sensitivity", main = "IGA EMA", xaxt="n" )
 lines(output$incremental_net_benefit[2:7], lwd=2)
   abline(h=0)
  axis(1,                         # Define x-axis manually
       at = 1:6,
       labels = c(0.5, 0.6, 0.7, 0.8, 0.9, 1))
  lines(output$incremental_net_benefit[8:13], col = 2, lwd=2, lty=2)
  lines(output$incremental_net_benefit[14:19], col = 3, lwd=2, lty=3)
  lines(output$incremental_net_benefit[20:25], col = 4, lwd=2, lty=4)
  lines(output$incremental_net_benefit[26:31], col = 5, lwd=2, lty=5)
  
  lines(output$inb_lci[2:7], lwd=1)
  lines(output$inb_lci[8:13], col = 2, lwd=1, lty=2)
  lines(output$inb_lci[14:19], col = 3, lwd=1, lty=3)
  lines(output$inb_lci[20:25], col = 4, lwd=1, lty=4)
  lines(output$inb_lci[26:31], col = 5, lwd=1, lty=5)
  
  lines(output$inb_uci[2:7], lwd=1)
  lines(output$inb_uci[8:13], col = 2, lwd=1, lty=2)
  lines(output$inb_uci[14:19], col = 3, lwd=1, lty=3)
  lines(output$inb_uci[20:25], col = 4, lwd=1, lty=4)
  lines(output$inb_uci[26:31], col = 5, lwd=1, lty=5)

  #IGA TTG plus EMA
  plot(output$incremental_net_benefit[38:43],pch=19, ylim=c(-15000,50000), ylab = "Incremental net benefit", main = "IGA TTG plus EMA",xlab = "sensitivity", xaxt="n" )
  lines(output$incremental_net_benefit[38:43], lwd=2)
  abline(h=0)
  axis(1,                         # Define x-axis manually
       at = 1:6,
       labels = c(0.5, 0.6, 0.7, 0.8, 0.9, 1))
  lines(output$incremental_net_benefit[44:49], col = 2, lwd = 2, lty = 2)
  lines(output$incremental_net_benefit[50:55], col = 3, lwd = 2, lty = 3)
  lines(output$incremental_net_benefit[56:61], col = 4, lwd = 2, lty = 4)
  lines(output$incremental_net_benefit[62:67], col = 5, lwd = 2, lty = 5)
  
  lines(output$inb_lci[38:43], lwd=1)
  lines(output$inb_lci[44:49], col = 2, lwd=1, lty=2)
  lines(output$inb_lci[50:55], col = 3, lwd=1, lty=3)
  lines(output$inb_lci[56:61], col = 4, lwd=1, lty=4)
  lines(output$inb_lci[62:67], col = 5, lwd=1, lty=5)
  
  lines(output$inb_uci[38:43], lwd=1)
  lines(output$inb_uci[44:49], col = 2, lwd=1, lty=2)
  lines(output$inb_uci[50:55], col = 3, lwd=1, lty=3)
  lines(output$inb_uci[56:61], col = 4, lwd=1, lty=4)
  lines(output$inb_uci[62:67], col = 5, lwd=1, lty=5)

  #IGA TTG
  plot(output$incremental_net_benefit[74:79],pch=19, ylim=c(-15000,50000), ylab = "Incremental net benefit", main = "IGA TTG", xlab = "sensitivity", xaxt="n" )
lines(output$incremental_net_benefit[74:79], lwd = 2)
abline(h=0)
axis(1,                         # Define x-axis manually
     at = 1:6,
     labels = c(0.5, 0.6, 0.7, 0.8, 0.9, 1))
  lines(output$incremental_net_benefit[80:85], col = 2, lwd = 2, lty = 2)
  lines(output$incremental_net_benefit[86:91], col = 3, lwd = 2, lty = 3)
  lines(output$incremental_net_benefit[92:97], col = 4, lwd = 2, lty = 4)
  lines(output$incremental_net_benefit[98:103], col = 5, lwd = 2, lty = 5)
  
  lines(output$inb_lci[74:79], lwd=1)
  lines(output$inb_lci[80:85], col = 2, lwd=1, lty=2)
  lines(output$inb_lci[86:91], col = 3, lwd=1, lty=3)
  lines(output$inb_lci[92:97], col = 4, lwd=1, lty=4)
  lines(output$inb_lci[98:103], col = 5, lwd=1, lty=5)
  
  lines(output$inb_uci[74:79], lwd=1)
  lines(output$inb_uci[80:85], col = 2, lwd=1, lty=2)
  lines(output$inb_uci[86:91], col = 3, lwd=1, lty=3)
  lines(output$inb_uci[92:97], col = 4, lwd=1, lty=4)
  lines(output$inb_uci[98:103], col = 5, lwd=1, lty=5)

#IgA EMA plus HLA 
   plot(output$incremental_net_benefit[110:115],pch=19, ylim=c(-15000, 50000),ylab = "Incremental net benefit", main = "IgA EMA plus HLA", xlab = "sensitivity", xaxt="n" )
   lines(output$incremental_net_benefit[110:115], lwd = 2)
   abline(h=0)
   axis(1,                         # Define x-axis manually
        at = 1:6,
        labels = c(0.5, 0.6, 0.7, 0.8, 0.9, 1))
  lines(output$incremental_net_benefit[116:121], col = 2, lwd = 2, lty = 2)
  lines(output$incremental_net_benefit[122:127], col = 3, lwd = 2, lty = 3)
  lines(output$incremental_net_benefit[128:133], col = 4, lwd = 2, lty = 4)
  lines(output$incremental_net_benefit[134:139], col = 5, lwd = 2, lty = 5)
  
  
  lines(output$inb_lci[110:115], lwd=1)
  lines(output$inb_lci[116:121], col = 2, lwd=1, lty=2)
  lines(output$inb_lci[122:127], col = 3, lwd=1, lty=3)
  lines(output$inb_lci[128:133], col = 4, lwd=1, lty=4)
  lines(output$inb_lci[134:139], col = 5, lwd=1, lty=5)
  
  lines(output$inb_uci[110:115], lwd=1)
  lines(output$inb_uci[116:121], col = 2, lwd=1, lty=2)
  lines(output$inb_uci[122:127], col = 3, lwd=1, lty=3)
  lines(output$inb_uci[128:133], col = 4, lwd=1, lty=4)
  lines(output$inb_uci[134:139], col = 5, lwd=1, lty=5)

#IgA TTG plus EMA plus HLA
   plot(output$incremental_net_benefit[146:151],pch=19, ylim=c(-15000, 50000),ylab = "Incremental net benefit", main = "IgATTG plus EMA plus HLA", xlab = "sensitivity", xaxt="n" )
   lines(output$incremental_net_benefit[146:151], lwd = 2)
   abline(h=0)
   axis(1,                         # Define x-axis manually
        at = 1:6,
        labels = c(0.5, 0.6, 0.7, 0.8, 0.9, 1))
  lines(output$incremental_net_benefit[152:157], col = 2, lwd = 2, lty = 2)
  lines(output$incremental_net_benefit[158:163], col = 3, lwd = 2, lty = 3)
  lines(output$incremental_net_benefit[164:169], col = 4, lwd = 2, lty = 4)
  lines(output$incremental_net_benefit[170:175], col = 5, lwd = 2, lty = 5)
  
  lines(output$inb_lci[146:151], lwd=1)
  lines(output$inb_lci[152:157], col = 2, lwd=1, lty=2)
  lines(output$inb_lci[158:163], col = 3, lwd=1, lty=3)
  lines(output$inb_lci[164:169], col = 4, lwd=1, lty=4)
  lines(output$inb_lci[170:175], col = 5, lwd=1, lty=5)
  
  
  lines(output$inb_uci[145:151], lwd=1)
  lines(output$inb_uci[152:157], col = 2, lwd=1, lty=2)
  lines(output$inb_uci[158:163], col = 3, lwd=1, lty=3)
  lines(output$inb_uci[164:169], col = 4, lwd=1, lty=4)
  lines(output$inb_uci[170:175], col = 5, lwd=1, lty=5)
 
  #IgA TTG plus HLA
  plot(output$incremental_net_benefit[182:187],pch=2, ylim=c(-15000, 50000),ylab = "Incremental net benefit", main = "IgA TTG plus HLA", xlab = "sensitivity", xaxt="n" )
  lines(output$incremental_net_benefit[182:187], lwd = 2)
  abline(h=0)
  axis(1,                         # Define x-axis manually
       at = 1:6,
       labels = c(0.5, 0.6, 0.7, 0.8, 0.9, 1))
  lines(output$incremental_net_benefit[188:193], col = 2, lwd = 2, lty = 2)
  lines(output$incremental_net_benefit[194:199], col = 3, lwd = 2, lty = 3)
  lines(output$incremental_net_benefit[200:205], col = 4, lwd = 2, lty = 4)
  lines(output$incremental_net_benefit[206:211], col = 5, lwd = 2, lty = 5)
  
  lines(output$inb_lci[182:187], lwd=1)
  lines(output$inb_lci[188:193], col = 2, lwd=1, lty=2)
  lines(output$inb_lci[194:199], col = 3, lwd=1, lty=3)
  lines(output$inb_lci[200:205], col = 4, lwd=1, lty=4)
  lines(output$inb_lci[206:211], col = 5, lwd=1, lty=5)
  
  lines(output$inb_uci[182:187], lwd=1)
  lines(output$inb_uci[188:193], col = 2, lwd=1, lty=2)
  lines(output$inb_uci[194:199], col = 3, lwd=1, lty=3)
  lines(output$inb_uci[200:205], col = 4, lwd=1, lty=4)
  lines(output$inb_uci[206:211], col = 5, lwd=1, lty=5)

dev.off()

jpeg("results/legend_inbplot.jpeg")
par(mfrow=c(1,1))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c('Specificity 0.5', '0.6', '0.7',
                            '0.8', '0.9'),
       col = c(1:5), lty = (1:5))
dev.off

#plots of probability CE
jpeg("results/inbplot.jpeg")
par(mfrow=c(3,2))
par(mar = c(2, 1, 1, 1))

#IGA EMA
plot(output$probability_cost_effective[2:7],pch=19, ylim=c(0,1), ylab = "Incremental net benefit", xlab = "sensitivity", main = "IGA EMA", xaxt="n" )
lines(output$probability_cost_effective[2:7], lwd=2)
abline(h=0)
axis(1,                         # Define x-axis manually
     at = 1:6,
     labels = c(0.5, 0.6, 0.7, 0.8, 0.9, 1))
lines(output$probability_cost_effective[8:13], col = 2, lwd=2, lty=2)
lines(output$probability_cost_effective[14:19], col = 3, lwd=, lty=3)
lines(output$probability_cost_effective[20:25], col = 4, lwd=2, lty=4)
lines(output$probability_cost_effective[26:31], col = 5, lwd=2, lty=5)

#IGA TTG plus EMA
plot(output$probability_cost_effective[38:43],pch=19, ylim=c(0,1), ylab = "Incremental net benefit", main = "IGA TTG plus EMA",xlab = "sensitivity", xaxt="n" )
lines(output$probability_cost_effective[38:43], lwd=2)
abline(h=0)
axis(1,                         # Define x-axis manually
     at = 1:6,
     labels = c(0.5, 0.6, 0.7, 0.8, 0.9, 1))
lines(output$probability_cost_effective[44:49], col = 2, lwd = 2, lty = 2)
lines(output$probability_cost_effective[50:55], col = 3, lwd = 2, lty = 3)
lines(output$probability_cost_effective[56:61], col = 4, lwd = 2, lty = 4)
lines(output$probability_cost_effective[62:67], col = 5, lwd = 2, lty = 5)

#IGA TTG
plot(output$probability_cost_effective[74:79],pch=19, ylim=c(0,1), ylab = "Incremental net benefit", main = "IGA TTG", xlab = "sensitivity", xaxt="n" )
lines(output$probability_cost_effective[74:79], lwd = 2)
abline(h=0)
axis(1,                         # Define x-axis manually
     at = 1:6,
     labels = c(0.5, 0.6, 0.7, 0.8, 0.9, 1))
lines(output$probability_cost_effective[80:85], col = 2, lwd = 2, lty = 2)
lines(output$probability_cost_effective[86:91], col = 3, lwd = 2, lty = 3)
lines(output$probability_cost_effective[92:97], col = 4, lwd = 2, lty = 4)
lines(output$probability_cost_effective[98:103], col = 5, lwd = 2, lty = 5)

#IgA EMA plus HLA 
plot(output$probability_cost_effective[110:115],pch=19, ylim=c(0, 1),ylab = "Incremental net benefit", main = "IgA EMA plus HLA", xlab = "sensitivity", xaxt="n" )
lines(output$probability_cost_effective[110:115], lwd = 2)
abline(h=0)
axis(1,                         # Define x-axis manually
     at = 1:6,
     labels = c(0.5, 0.6, 0.7, 0.8, 0.9, 1))
lines(output$probability_cost_effective[116:121], col = 2, lwd = 2, lty = 2)
lines(output$probability_cost_effective[122:127], col = 3, lwd = 2, lty = 3)
lines(output$probability_cost_effective[128:133], col = 4, lwd = 2, lty = 4)
lines(output$probability_cost_effective[134:139], col = 5, lwd = 2, lty = 5)

#IgA TTG plus EMA plus HLA
plot(output$probability_cost_effective[146:151],pch=19, ylim=c(0, 1),ylab = "Incremental net benefit", main = "IgATTG plus EMA plus HLA", xlab = "sensitivity", xaxt="n" )
lines(output$probability_cost_effective[146:151], lwd = 2)
abline(h=0)
axis(1,                         # Define x-axis manually
     at = 1:6,
     labels = c(0.5, 0.6, 0.7, 0.8, 0.9, 1))
lines(output$probability_cost_effective[152:157], col = 2, lwd = 2, lty = 2)
lines(output$probability_cost_effective[158:163], col = 3, lwd = 2, lty = 3)
lines(output$probability_cost_effective[164:169], col = 4, lwd = 2, lty = 4)
lines(output$probability_cost_effective[170:175], col = 5, lwd = 2, lty = 5)

#IgA TTG plus HLA
plot(output$probability_cost_effective[182:187],pch=2, ylim=c(0, 1),ylab = "Incremental net benefit", main = "IgA TTG plus HLA", xlab = "sensitivity", xaxt="n" )
lines(output$probability_cost_effective[182:187], lwd = 2)
abline(h=0)
axis(1,                         # Define x-axis manually
     at = 1:6,
     labels = c(0.5, 0.6, 0.7, 0.8, 0.9, 1))
lines(output$probability_cost_effective[188:193], col = 2, lwd = 2, lty = 2)
lines(output$probability_cost_effective[194:199], col = 3, lwd = 2, lty = 3)
lines(output$probability_cost_effective[200:205], col = 4, lwd = 2, lty = 4)
lines(output$probability_cost_effective[206:211], col = 5, lwd = 2, lty = 5)

dev.off()

write.csv(t(output$total_qalys), "results/total qalys.csv")
write.csv(t(output$total_costs), "results/total costs.csv")


pCE <- as.data.frame(output$probability_cost_effective)
pCE$sens <- c(0,rep(sens_riskfactor, n_sero_tests*6))
pCE$spec <- c(0,rep(rep(spec_riskfactor, each=n_sero_tests),6))
  
#heat plot of cost-effectiveness
  ggplot(data = pCE, aes(x = spec,
             
             y = sens)) +
  
  geom_tile(aes(fill=output$probability_cost_effective)) +
  
  scale_fill_gradient2(high     = "red4",
                       
                       mid      = "white",
                       
                       low      = "midnightblue",
                       
                       midpoint = 0.5) +
  
  scale_x_reverse(expand = c(0, 0)) +
  
  scale_y_continuous(expand = c(0, 0)) 


  #table of costs and QALYs
  format.results<-function(x,n.digits=2)
  {
    paste(round(mean(x),digits=n.digits)," (",round(quantile(x,probs=0.025),digits=n.digits),", ",round(quantile(x,probs=0.975),digits=n.digits),")",sep="")
  }
  
  feasible.strategies.qalys <- output$total_qalys[c("No screening", "0.6 0.99 IgATTG", "0.6 0.99 IgAEMA", "0.6 0.99 IgATTGplusEMA", "0.6 0.99 IgATTG plus HLA", 
                                                    "0.6 0.99 IgAEMA plus HLA", "0.6 0.99 IgATTGplusEMA plus HLA", "0.99 0.99 IgATTG", 
                                                    "0.99 0.99 IgAEMA", "0.99 0.99 IgATTGplusEMA", "0.99 0.99 IgATTG plus HLA", 
                                                    "0.99 0.99 IgAEMA plus HLA", "0.99 0.99 IgATTGplusEMA plus HLA"),]
  feasible.strategies.costs <- output$total_costs[c("No screening", "0.6 0.99 IgATTG", "0.6 0.99 IgAEMA", "0.6 0.99 IgATTGplusEMA", "0.6 0.99 IgATTG plus HLA", 
                                                     "0.6 0.99 IgAEMA plus HLA", "0.6 0.99 IgATTGplusEMA plus HLA", "0.99 0.99 IgATTG", 
                                                     "0.99 0.99 IgAEMA", "0.99 0.99 IgATTGplusEMA", "0.99 0.99 IgATTG plus HLA", 
                                                     "0.99 0.99 IgAEMA plus HLA", "0.99 0.99 IgATTGplusEMA plus HLA"),]
# feasible.strategies.netbenefit <- output$all_net_benefit[c("No screening", "0.6 0.99 IgATTG", "0.6 0.99 IgAEMA", "0.6 0.99 IgATTGplusEMA", "0.6 0.99 IgATTG plus HLA", 
 #                                                           "0.6 0.99 IgAEMA plus HLA", "0.6 0.99 IgATTGplusEMA plus HLA", "0.99 0.99 IgATTG", 
  #                                                          "0.99 0.99 IgAEMA", "0.99 0.99 IgATTGplusEMA", "0.99 0.99 IgATTG plus HLA", 
   #                                                         "0.99 0.99 IgAEMA plus HLA", "0.99 0.99 IgATTGplusEMA plus HLA"),]
 feasible.strategies.inetbenefit <- output$all_incremental_net_benefit[c("No screening", "0.6 0.99 IgATTG", "0.6 0.99 IgAEMA", "0.6 0.99 IgATTGplusEMA", "0.6 0.99 IgATTG plus HLA", 
                                                            "0.6 0.99 IgAEMA plus HLA", "0.6 0.99 IgATTGplusEMA plus HLA", "0.99 0.99 IgATTG", 
                                                            "0.99 0.99 IgAEMA", "0.99 0.99 IgATTGplusEMA", "0.99 0.99 IgATTG plus HLA", 
                                                            "0.99 0.99 IgAEMA plus HLA", "0.99 0.99 IgATTGplusEMA plus HLA"),]
   feasible.strategies.qalys.table <- array(dim=c(1, 13), dimnames = list(NULL, NULL))
  feasible.strategies.costs.table <- array(dim=c(1, 13), dimnames = list(NULL, NULL))
 # feasible.strategies.netbenefit.table <- array(dim=c(1, 13), dimnames = list(NULL, NULL))
  feasible.strategies.inetbenefit.table <- array(dim=c(1, 13), dimnames = list(NULL, NULL))

 for (i in 1:13) { 
   feasible.strategies.qalys.table[,i] <- format.results(feasible.strategies.qalys[i,], n.digits = 4)
   feasible.strategies.costs.table[,i] <- format.results(feasible.strategies.costs[i,], n.digits = 0)
  # feasible.strategies.netbenefit.table[,i] <- format.results(feasible.strategies.netbenefit[i,], n.digits = 0)
   feasible.strategies.inetbenefit.table[,i] <- format.results(feasible.strategies.inetbenefit[i,], n.digits = 0)
 }
  
  table <- data.frame(t(feasible.strategies.costs.table), t(feasible.strategies.qalys.table), t(feasible.strategies.inetbenefit.table))
  colnames(table) <- c("Costs", "QALYs", "Incremental net benefit v no screening")
  rownames(table) <- c("No screening", "0.6 0.99 IgATTG", "0.6 0.99 IgAEMA", "0.6 0.99 IgATTGplusEMA", "0.6 0.99 IgATTG plus HLA", 
                       "0.6 0.99 IgAEMA plus HLA", "0.6 0.99 IgATTGplusEMA plus HLA", "0.99 0.99 IgATTG", 
                       "0.99 0.99 IgAEMA", "0.99 0.99 IgATTGplusEMA", "0.99 0.99 IgATTG plus HLA", 
                       "0.99 0.99 IgAEMA plus HLA", "0.99 0.99 IgATTGplusEMA plus HLA")
  table <-  table[order( table$`Incremental net benefit v no screening`),]
  write.csv(table, "results/table of results.csv")
 
  #for CEAC, if doing in excel. can also do using BCEA
   write.csv(t(output$total_qalys[c("No screening", "0.6 0.99 IgATTG", "0.6 0.99 IgAEMA", "0.6 0.99 IgATTGplusEMA", "0.6 0.99 IgATTG plus HLA", 
                                 "0.6 0.99 IgAEMA plus HLA", "0.6 0.99 IgATTGplusEMA plus HLA", "0.99 0.99 IgATTG", 
                                 "0.99 0.99 IgAEMA", "0.99 0.99 IgATTGplusEMA", "0.99 0.99 IgATTG plus HLA", 
                                 "0.99 0.99 IgAEMA plus HLA", "0.99 0.99 IgATTGplusEMA plus HLA"),]), "qalys_feasible.csv")
  write.csv(t(output$total_costs[c("No screening", "0.6 0.99 IgATTG", "0.6 0.99 IgAEMA", "0.6 0.99 IgATTGplusEMA", "0.6 0.99 IgATTG plus HLA", 
                                 "0.6 0.99 IgAEMA plus HLA", "0.6 0.99 IgATTGplusEMA plus HLA", "0.99 0.99 IgATTG", 
                                 "0.99 0.99 IgAEMA", "0.99 0.99 IgATTGplusEMA", "0.99 0.99 IgATTG plus HLA", 
                                 "0.99 0.99 IgAEMA plus HLA", "0.99 0.99 IgATTGplusEMA plus HLA"),]), "costs_feasible.csv")
  #cost breakdown
  feasible.strategies.testcosts <- output$test_costs_applied[,c("No screening", "0.6 0.99 IgATTG", "0.6 0.99 IgAEMA", "0.6 0.99 IgATTGplusEMA", "0.6 0.99 IgATTG plus HLA", 
                                                    "0.6 0.99 IgAEMA plus HLA", "0.6 0.99 IgATTGplusEMA plus HLA", "0.99 0.99 IgATTG", 
                                                    "0.99 0.99 IgAEMA", "0.99 0.99 IgATTGplusEMA", "0.99 0.99 IgATTG plus HLA", 
                                                    "0.99 0.99 IgAEMA plus HLA", "0.99 0.99 IgATTGplusEMA plus HLA")]
  feasible.strategies.fpcosts <- output$fp_costs[,c("No screening", "0.6 0.99 IgATTG", "0.6 0.99 IgAEMA", "0.6 0.99 IgATTGplusEMA", "0.6 0.99 IgATTG plus HLA", 
                                                                "0.6 0.99 IgAEMA plus HLA", "0.6 0.99 IgATTGplusEMA plus HLA", "0.99 0.99 IgATTG", 
                                                                "0.99 0.99 IgAEMA", "0.99 0.99 IgATTGplusEMA", "0.99 0.99 IgATTG plus HLA", 
                                                                "0.99 0.99 IgAEMA plus HLA", "0.99 0.99 IgATTGplusEMA plus HLA")]
  feasible.strategies.diagnosiscosts <- output$diagnosis_costs[,c("No screening", "0.6 0.99 IgATTG", "0.6 0.99 IgAEMA", "0.6 0.99 IgATTGplusEMA", "0.6 0.99 IgATTG plus HLA", 
                                                                        "0.6 0.99 IgAEMA plus HLA", "0.6 0.99 IgATTGplusEMA plus HLA", "0.99 0.99 IgATTG", 
                                                                        "0.99 0.99 IgAEMA", "0.99 0.99 IgATTGplusEMA", "0.99 0.99 IgATTG plus HLA", 
                                                                        "0.99 0.99 IgAEMA plus HLA", "0.99 0.99 IgATTGplusEMA plus HLA")]
 
  feasible.strategies.cyclecosts <- output$cycle_costs[c("No screening", "0.6 0.99 IgATTG", "0.6 0.99 IgAEMA", "0.6 0.99 IgATTGplusEMA", "0.6 0.99 IgATTG plus HLA", 
                                                                  "0.6 0.99 IgAEMA plus HLA", "0.6 0.99 IgATTGplusEMA plus HLA", "0.99 0.99 IgATTG", 
                                                                 "0.99 0.99 IgAEMA", "0.99 0.99 IgATTGplusEMA", "0.99 0.99 IgATTG plus HLA", 
                                                                  "0.99 0.99 IgAEMA plus HLA", "0.99 0.99 IgATTGplusEMA plus HLA")]
   feasible.strategies.testcosts.table <- array(dim=c(1, 13), dimnames = list(NULL, NULL))
  feasible.strategies.fpcosts.table <- array(dim=c(1, 13), dimnames = list(NULL, NULL))
  feasible.strategies.diagnosiscosts.table <- array(dim=c(1, 13), dimnames = list(NULL, NULL))
  feasible.strategies.cyclecosts.table <- array(dim=c(1, 13), dimnames = list(NULL, NULL))
  for (i in 1:13) { 
    feasible.strategies.testcosts.table[,i] <- format.results(feasible.strategies.testcosts[,i])
    feasible.strategies.fpcosts.table[,i] <- format.results(feasible.strategies.fpcosts[,i])
    feasible.strategies.diagnosiscosts.table[,i] <- format.results(feasible.strategies.diagnosiscosts[,i])
   feasible.strategies.cyclecosts.table[i] <- format.results(feasible.strategies.cyclecosts[i])
  }
  
  cost_breakdown_table <- data.frame(t(feasible.strategies.testcosts.table), t(feasible.strategies.fpcosts.table), t(feasible.strategies.diagnosiscosts.table), t(feasible.strategies.cyclecosts.table))
  colnames(cost_breakdown_table) <- c("Test costs", "False positive costs", "Diagnosis costs", "Cycle costs")
  rownames(cost_breakdown_table) <- c("No screening", "0.6 0.99 IgATTG", "0.6 0.99 IgAEMA", "0.6 0.99 IgATTGplusEMA", "0.6 0.99 IgATTG plus HLA", 
                       "0.6 0.99 IgAEMA plus HLA", "0.6 0.99 IgATTGplusEMA plus HLA", "0.99 0.99 IgATTG", 
                       "0.99 0.99 IgAEMA", "0.99 0.99 IgATTGplusEMA", "0.99 0.99 IgATTG plus HLA", 
                       "0.99 0.99 IgAEMA plus HLA", "0.99 0.99 IgATTGplusEMA plus HLA")
  
  write.csv(cost_breakdown_table, "results/cost breakdown table.csv")
 
  
  #utility breakdown
  feasible.strategies.disutilitybiopsy <- output$disutility_biopsy[,c("No screening", "0.6 0.99 IgATTG", "0.6 0.99 IgAEMA", "0.6 0.99 IgATTGplusEMA", "0.6 0.99 IgATTG plus HLA", 
                                                                "0.6 0.99 IgAEMA plus HLA", "0.6 0.99 IgATTGplusEMA plus HLA", "0.99 0.99 IgATTG", 
                                                                "0.99 0.99 IgAEMA", "0.99 0.99 IgATTGplusEMA", "0.99 0.99 IgATTG plus HLA", 
                                                                "0.99 0.99 IgAEMA plus HLA", "0.99 0.99 IgATTGplusEMA plus HLA")]
  feasible.strategies.disutilitybiopsywait <- output$disutility_biopsy_wait[,c("No screening", "0.6 0.99 IgATTG", "0.6 0.99 IgAEMA", "0.6 0.99 IgATTGplusEMA", "0.6 0.99 IgATTG plus HLA", 
                                                    "0.6 0.99 IgAEMA plus HLA", "0.6 0.99 IgATTGplusEMA plus HLA", "0.99 0.99 IgATTG", 
                                                    "0.99 0.99 IgAEMA", "0.99 0.99 IgATTGplusEMA", "0.99 0.99 IgATTG plus HLA", 
                                                    "0.99 0.99 IgAEMA plus HLA", "0.99 0.99 IgATTGplusEMA plus HLA")]
  feasible.strategies.disutilityfp <- output$disutility_fp[,c("No screening", "0.6 0.99 IgATTG", "0.6 0.99 IgAEMA", "0.6 0.99 IgATTGplusEMA", "0.6 0.99 IgATTG plus HLA", 
                                                                  "0.6 0.99 IgAEMA plus HLA", "0.6 0.99 IgATTGplusEMA plus HLA", "0.99 0.99 IgATTG", 
                                                                  "0.99 0.99 IgAEMA", "0.99 0.99 IgATTGplusEMA", "0.99 0.99 IgATTG plus HLA", 
                                                                  "0.99 0.99 IgAEMA plus HLA", "0.99 0.99 IgATTGplusEMA plus HLA")]
  
  feasible.strategies.cycleqalys <- output$cycle_qalys[c("No screening", "0.6 0.99 IgATTG", "0.6 0.99 IgAEMA", "0.6 0.99 IgATTGplusEMA", "0.6 0.99 IgATTG plus HLA", 
                                                         "0.6 0.99 IgAEMA plus HLA", "0.6 0.99 IgATTGplusEMA plus HLA", "0.99 0.99 IgATTG", 
                                                         "0.99 0.99 IgAEMA", "0.99 0.99 IgATTGplusEMA", "0.99 0.99 IgATTG plus HLA", 
                                                         "0.99 0.99 IgAEMA plus HLA", "0.99 0.99 IgATTGplusEMA plus HLA")]
  feasible.strategies.disutilitybiopsy.table <- array(dim=c(1, 13), dimnames = list(NULL, NULL))
  feasible.strategies.disutilitybiopsywait.table <- array(dim=c(1, 13), dimnames = list(NULL, NULL))
  feasible.strategies.disutilityfp.table <- array(dim=c(1, 13), dimnames = list(NULL, NULL))
  feasible.strategies.cycleqalys.table <- array(dim=c(1, 13), dimnames = list(NULL, NULL))
  for (i in 1:13) { 
    feasible.strategies.disutilitybiopsy.table[,i] <- format.results(feasible.strategies.disutilitybiopsy[,i], n.digits = 4)
    feasible.strategies.disutilitybiopsywait.table[,i] <- format.results(feasible.strategies.disutilitybiopsywait[,i], n.digits = 4)
    feasible.strategies.disutilityfp.table[,i] <- format.results(feasible.strategies.disutilityfp[,i], n.digits = 4)
    feasible.strategies.cycleqalys.table[i] <- format.results(feasible.strategies.cycleqalys[i], n.digits = 4)
  }
  
  qaly_breakdown_table <- data.frame(t(feasible.strategies.cycleqalys.table), t(feasible.strategies.disutilityfp.table), t(feasible.strategies.disutilitybiopsywait.table), t(feasible.strategies.disutilitybiopsy.table))
  colnames(qaly_breakdown_table) <- c("Cycle QALYs", "Disutility FP", "Disutility waiting for biopsy", "Disutility biopsy")
  rownames(qaly_breakdown_table) <- c("No screening", "0.6 0.99 IgATTG", "0.6 0.99 IgAEMA", "0.6 0.99 IgATTGplusEMA", "0.6 0.99 IgATTG plus HLA", 
                                      "0.6 0.99 IgAEMA plus HLA", "0.6 0.99 IgATTGplusEMA plus HLA", "0.99 0.99 IgATTG", 
                                      "0.99 0.99 IgAEMA", "0.99 0.99 IgATTGplusEMA", "0.99 0.99 IgATTG plus HLA", 
                                      "0.99 0.99 IgAEMA plus HLA", "0.99 0.99 IgATTGplusEMA plus HLA")
  
  write.csv(qaly_breakdown_table, "results/qaly breakdown table.csv")
  
  #time in states
  
  time_in_states <- output$time_in_states
  
  write.csv(output$time_in_states, "results/time in states.csv")
  
   # Now use the BCEA package to analyse the results
 # pkgs <- c("MASS","Rtools","devtools")
  #repos <- c("https://cran.rstudio.com", "https://www.math.ntnu.no/inla/R/stable") 
  #install.packages(pkgs,repos=repos,dependencies = "Depends")
  #devtools::install_github("giabaio/BCEA")


m <- bcea(e = t(output$total_qalys), c = t(output$total_costs), ref = 1, interventions = t_names)
summary(m)

m_feasible <- bcea(e = t(output$total_qalys[c("No screening", "0.6 0.99 IgATTG", "0.6 0.99 IgAEMA", "0.6 0.99 IgATTGplusEMA", "0.6 0.99 IgATTG plus HLA", 
                     "0.6 0.99 IgAEMA plus HLA", "0.6 0.99 IgATTGplusEMA plus HLA", "0.99 0.99 IgATTG", 
                     "0.99 0.99 IgAEMA", "0.99 0.99 IgATTGplusEMA", "0.99 0.99 IgATTG plus HLA", 
                     "0.99 0.99 IgAEMA plus HLA", "0.99 0.99 IgATTGplusEMA plus HLA"),]), 
                   c = t(output$total_costs[c("No screening", "0.6 0.99 IgATTG", "0.6 0.99 IgAEMA", "0.6 0.99 IgATTGplusEMA", 
                                              "0.6 0.99 IgATTG plus HLA", 
                                            "0.6 0.99 IgAEMA plus HLA", "0.6 0.99 IgATTGplusEMA plus HLA", "0.99 0.99 IgATTG", 
                                             "0.99 0.99 IgAEMA", "0.99 0.99 IgATTGplusEMA", "0.99 0.99 IgATTG plus HLA", 
                                           "0.99 0.99 IgAEMA plus HLA", "0.99 0.99 IgATTGplusEMA plus HLA"),]), ref=1, 
                   interventions = c("No screening", "0.6 0.99 IgATTG", "0.6 0.99 IgAEMA", "0.6 0.99 IgATTGplusEMA", "0.6 0.99 IgATTG plus HLA", 
                                    "0.6 0.99 IgAEMA plus HLA", "0.6 0.99 IgATTGplusEMA plus HLA", "0.99 0.99 IgATTG", 
                                  "0.99 0.99 IgAEMA", "0.99 0.99 IgATTGplusEMA", "0.99 0.99 IgATTG plus HLA", 
                                              "0.99 0.99 IgAEMA plus HLA", "0.99 0.99 IgATTGplusEMA plus HLA"))
summary(m_feasible) 


eib.plot(m_feasible, comparison = NULL, pos =
           c(1, 0), size = NULL, plot.cri = NULL, graph
         = c("ggplot2"))
evi.plot(m_feasible, graph = c("base", "ggplot2",
                       "plotly"))
ceac.plot(m_feasible, comparison = NULL,
          pos = FALSE, graph = c("ggplot2"))
par(mfrow = c(1,1))
mce <- multi.ce(m_feasible)
ceaf.plot(mce, graph = c("ggplot2"))
mce.plot(mce, color = c(1:13)) #CEAC 

ceplane.plot(m_feasible, comparison =
               NULL, pos = c(1, 0), graph = c("ggplot2"), point_colors = c(1:13))
sim.table(m_feasible)
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
utility_parameters <- data.frame(input_parameters$utility_GFD, input_parameters$utility_undiagnosedCD, input_parameters$disutility_subfertility, 
                                 input_parameters$disutility_osteoporosis, input_parameters$disutility_NHL, input_parameters$disutility_biopsy, 
                                 input_parameters$disutility_biopsy_wait, input_parameters$disutility_fp)
cost_parameters <- data.frame( input_parameters$cost_CDGFD, input_parameters$cost_osteoporosis, input_parameters$cost_undiagnosedCD, 
                               input_parameters$cost_IDA, input_parameters$cost_biopsy, 
                               input_parameters$cost_subfertility, input_parameters$cost_NHL,input_parameters$cost_diagnosis,
                               input_parameters$test_cost_IgAEMA, input_parameters$test_cost_IgATTG, input_parameters$test_cost_HLA)
info.rank(parameter = colnames(input_parameters), input = input_parameters, m_feasible, xlim = c(0,0.5), wtp=30000)
info.rank(parameter = colnames(utility_parameters), input = utility_parameters, m_feasible, xlim = c(0,0.5), wtp=30000)
info.rank(parameter = colnames(cost_parameters), input = cost_parameters, m_feasible, xlim = c(0,0.5), wtp=30000)



