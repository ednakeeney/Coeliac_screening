# Script to run MLMC for EVPPI
# Needs to be embedded within coeliac_markov_main script and runs after info.rank
# Howard Thom 10-November-2021

require(ggplot2)
require(grid)
require(Rcpp)
require(tictoc)
require(doRNG)

# Required for structure of net benefits in level l estimator

n_treatment <- length(strategies_of_interest)

# General repackaging of model functions
generate_NetB <- function(n_samples, input_parameters_mlmc) {
  # Generate transition matrices for given input parameters
  transition_matrices_mlmc <- generate_transition_matrices(input_parameters_mlmc, population = population)
  
  # Generate net benefit for strategies of interest
  output <- generate_net_benefit(input_parameters_mlmc, strategies_of_interest, 
                                 transition_matrices_mlmc,
                                 combinations = combinations,
                                 population = population)
  
  # Needs to be samples x strategies
  NetB <- t(output$all_net_benefit[strategies_of_interest, ])
  return(NetB)
}

###########################################################################

# Wrapper function to generate the net benefit holding the parameters of interest constant
# Repeats parameters of interest M times and generates random for remainder
# Calculates net benefit based on these
EVPPI_utilities_std_p<-function(M, N)
{
  # N is the number of outer samples, M=2^l is the number of inner samples
  # Total number of samples is NN
  NN <- M * N
  
  input_parameters <- generate_model_parameters(starting_age, population = population,
                                                combinations = combinations, n_samples = NN,
                                                hold_constant =   c("utility_undiagnosedCD", 
                                                                    "disutility_osteoporosis", "disutility_NHL",
                                                                    "disutility_biopsy", "disutility_biopsy_wait",
                                                                    "disutility_fp"))
  
  NetB <- generate_NetB(n_samples = NN, input_parameters_mlmc = input_parameters)

  return(NetB)
}

# Wrapper function for net benefit holding rates constant
EVPPI_rates_std_p<-function(M, N)
{
  # N is the number of outer samples, M=2^l is the number of inner samples
  # Total number of samples is NN
  NN <- M * N
  
  input_parameters <- generate_model_parameters(starting_age, population = population,
                                                combinations = combinations, n_samples = NN,
                                                hold_constant =  colnames(input_parameters)[grepl("rate", colnames(input_parameters))])
  
  NetB <- generate_NetB(n_samples = NN, input_parameters_mlmc = input_parameters)
  
  return(NetB)
}
 

# Wrapper function for net benefit holding GFD effects constant
EVPPI_gfd_std_p<-function(M, N)
{
  # N is the number of outer samples, M=2^l is the number of inner samples
  # Total number of samples is NN
  NN <- M * N
  
  input_parameters <- generate_model_parameters(starting_age, population = population,
                                                combinations = combinations, n_samples = NN,
                                                hold_constant =  c("utility_GFD",
                                                                   "log_or_osteoporosis_GFD",
                                                                   "log_or_osteoporosis_noGFD", 
                                                                   "log_rr_NHL_GFD" ,
                                                                   "log_rr_NHL_noGFD"))
  
  NetB <- generate_NetB(n_samples = NN, input_parameters_mlmc = input_parameters)
  
  return(NetB)
}

# Wrapper function for net benefit holding test accuracies constant
EVPPI_test_accuracy_std_p<-function(M, N)
{
  # N is the number of outer samples, M=2^l is the number of inner samples
  # Total number of samples is NN
  NN <- M * N
  
  input_parameters <- generate_model_parameters(starting_age, population = population,
                                                combinations = combinations, n_samples = NN,
                                                hold_constant =   colnames(input_parameters)[grepl("sens", colnames(input_parameters)) | grepl("spec", colnames(input_parameters))])
  
  NetB <- generate_NetB(n_samples = NN, input_parameters_mlmc = input_parameters)
  
  return(NetB)
}

# Wrapper function for net benefit holding probability of late diagnosis constant
EVPPI_late_diagnosis_std_p<-function(M, N)
{
  # N is the number of outer samples, M=2^l is the number of inner samples
  # Total number of samples is NN
  NN <- M * N
  
  input_parameters <- generate_model_parameters(starting_age, population = population,
                                                combinations = combinations, n_samples = NN,
                                                hold_constant =   c("probability_late_diagnosis"))
  
  NetB <- generate_NetB(n_samples = NN, input_parameters_mlmc = input_parameters)
  
  return(NetB)
}
###########################################################################
# Wrapper function for level l estimator of utilities
EVPPI_utilities_l_p <- function(l = l,N = N) {
  return(EVPPI_l_p(l, N, EVPPI_std_p = EVPPI_utilities_std_p))
}

# Wrapper function for level l estimator of rates of osteoporosis and NHL
EVPPI_rates_l_p <- function(l = l,N = N) {
  return(EVPPI_l_p(l, N, EVPPI_std_p = EVPPI_rates_std_p))
}
# Wrapper function for level l estimator of GFD effect
EVPPI_gfd_l_p <- function(l = l,N = N) {
  return(EVPPI_l_p(l, N, EVPPI_std_p = EVPPI_gfd_std_p))
}
# Wrapper function for level l estimator of test accuracies
EVPPI_test_accuracy_l_p <- function(l = l,N = N) {
  return(EVPPI_l_p(l, N, EVPPI_std_p = EVPPI_test_accuracy_std_p))
}
# Wrapper function for level l estimator of test accuracies
EVPPI_late_diagnosis_l_p <- function(l = l,N = N) {
  return(EVPPI_l_p(l, N, EVPPI_std_p = EVPPI_late_diagnosis_std_p))
}

# Tests of the function
#EVPPI_utilities_std_p(10, 10)
#EVPPI_utilities_l_p(l = 1, N = 10)
#EVPPI_rates_l_p(l = 1, N = 10)
#EVPPI_gfd_l_p(l = 1, N = 10)
#EVPPI_test_accuracy_l_p(l = 1, N = 10)
#EVPPI_late_diagnosis_l_p(l = 1, N = 10)

###########################################################################
# EVPPI of utilities
print("Utilities.........................................")
tic()
tst_utilities <- mlmc.test(EVPPI_utilities_l_p, M=2, N=16,
                   L=5, N0=16,
                   eps.v=c(60, 30),
                   Lmin=2, Lmax=10)
toc()
evpi_table["Utilities and disutilities", c("Per person", "Population")] <- (m$evi - tst_utilities$P[length(tst_utilities$P)]) * c(1, discounted_population_size)
write.csv(evpi_table, paste0("results/", population, "/", population, "_evppi_summary.csv"))

# EVPPI of rates
print("Rates.........................................")
tic()
tst_rates <- mlmc.test(EVPPI_rates_l_p, M=2, N=128,
                           L=2, N0=128,
                           eps.v=c(60,30,15),
                           Lmin=2, Lmax=10)
toc()
evpi_table["Rates of osteoporosis and NHL", c("Per person", "Population")] <- (m$evi - tst_rates$P[length(tst_rates$P)]) * c(1, discounted_population_size)
write.csv(evpi_table, paste0("results/", population, "/", population, "_evppi_summary.csv"))

# EVPPI of GFD effect
print("GFD.........................................")
tic()
tst_gfd <- mlmc.test(EVPPI_gfd_l_p, M=2, N=128,
                       L=2, N0=128,
                       eps.v=c(60,30,15),
                       Lmin=2, Lmax=10)
toc()
evpi_table["GFD effect", c("Per person", "Population")] <- (m$evi - tst_gfd$P[length(tst_gfd$P)]) * c(1, discounted_population_size)
write.csv(evpi_table, paste0("results/", population, "/", population, "_evppi_summary.csv"))


# EVPPI of test accuracies
print("Test accuracies.........................................")
tic()
tst_test_accuracy <- mlmc.test(EVPPI_test_accuracy_l_p, M=2, N=128,
                     L=2, N0=128,
                     eps.v=c(60,30,15),
                     Lmin=2, Lmax=10)
toc()
evpi_table["Test accuracies", c("Per person", "Population")] <- (m$evi - tst_test_accuracy$P[length(tst_test_accuracy$P)]) * c(1, discounted_population_size)
write.csv(evpi_table, paste0("results/", population, "/", population, "_evppi_summary.csv"))


# EVPPI of probability of late diagnosis
print("Late diagnosis.........................................")
tic()
tst_late_diagnosis <- mlmc.test(EVPPI_late_diagnosis_l_p, M=2, N=128,
                               L=2, N0=128,
                               eps.v=c(60,30,15),
                               Lmin=2, Lmax=10)
toc()
evpi_table["Probability of late diagnosis", c("Per person", "Population")] <- (m$evi - tst_late_diagnosis$P[length(tst_late_diagnosis$P)]) * c(1, discounted_population_size)
write.csv(evpi_table, paste0("results/", population, "/", population, "_evppi_summary.csv"))



# Tests of code
#tst_utilities$P[length(tst_utilities$P)]


###########################################################################



# Export results
write.csv(evpi_table, paste0("results/", population, "/", population, "_evppi_summary.csv"))

