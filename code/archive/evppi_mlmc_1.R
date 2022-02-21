# Script to run MLMC for EVPPI
# Needs to be embedded within coeliac_markov_main script and runs after info.rank
# Howard Thom 10-November-2021

require(ggplot2)
require(grid)
require(Rcpp)
require(tictoc)
require(doRNG)
#
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

n_treatment <- length(strategies_of_interest)

# Wrapper function for level l estimator of utilities
EVPPI_utilities_l_p <- function(l = l,N = N) {
  return(EVPPI_l_p(l, N, EVPPI_std_p = EVPPI_utilities_std_p))
}

EVPPI_utilities_std_p(10, 10)
EVPPI_utilities_l_p(l = 1, N = 10)

tic()
tst_utilities <- mlmc.test(EVPPI_utilities_l_p, M=2, N=128,
                   L=2, N0=128,
                   eps.v=c(60,30,15),
                   Lmin=2, Lmax=10)
toc()

# Error message
EVPPI_utilities_std_p(8, 16)


# DIFF with
tst_utilities$P[length(tst_utilities$P)]
old <- tst_utilities
