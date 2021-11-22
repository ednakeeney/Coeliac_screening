# This code is derived and adapted from the original GPL-2 Matlab version by
# Mike Giles.  See http://people.maths.ox.ac.uk/~gilesm/mlmc/

#' Multi-level Monte Carlo estimation
#'
#' This function is the Multi-level Monte Carlo driver which will sample from the levels of user specified function.
#'
#' Multilevel Monte Carlo Method method originated in works Giles (2008) and Heinrich (1998).
#'
#' Consider  a sequence \eqn{P_0, P_1, \ldots,} which approximates \eqn{P_L} with
#' increasing accuracy, but also increasing cost, we have the simple identity
#' \deqn{E[P_L] = E[P_0] + \sum_{l=1}^L E[P_l-P_{l-1}],}
#' and therefore we can use the following unbiased estimator for \eqn{E[P_L]},
#' \deqn{N_0^{-1} \sum_{n=1}^{N_0} P_0^{(0,n)} +
#' \sum_{l=1}^L \{
#' N_l^{-1} \sum_{n=1}^{N_l} (P_l^{(l,n)} - P_{l-1}^{(l,n)})
#' \}}
#' with the inclusion of the level \eqn{l} in the superscript \eqn{(l,n)} indicating
#' that the samples used at each level of correction are independent.
#'
#' Set \eqn{C_0}, and \eqn{V_0} to be the cost and variance of one sample of \eqn{P_0},
#' and \eqn{C_l, V_l} to be the cost and variance of one sample of \eqn{P_l - P_{l-1}},
#' then the overall cost and variance of the multilevel estimator is
#' \eqn{\sum_{l=0}^L N_l C_l\}}
#' and
#' \eqn{\sum_{l=0}^L N_l^{-1} V_l},
#' respectively.
#'
#' The idea begind the method, is that provided that the product \eqn{V_l C_l}
#' decreases with \eqn{l}, i.e. the cost increases with level slower than the
#' variance decreases, then one can achieve significant computational savings,
#' which can be formalised as in Theorem 1 of Giles (2015).
#'
#' For further information on multilevel Monte Carlo methods, see the webpage
#' \url{http://people.maths.ox.ac.uk/gilesm/mlmc_community.html}
#' which lists the research groups working in the area, and their
#' main publications.
#'
#' This function is based on GPL-2 'Matlab' code by Mike Giles.
#'
#' @author Louis Aslett <aslett@stats.ox.ac.uk>
#' @author Mike Giles <Mike.Giles@maths.ox.ac.uk>
#' @author Tigran Nagapetyan <nagapetyan@stats.ox.ac.uk>
#'
#' @references
#' M.B. Giles. Multilevel Monte Carlo path simulation. \emph{Operations Research}, 56(3):607-617, 2008.
#'
#' M.B. Giles. Multilevel Monte Carlo methods. \emph{Acta Numerica}, 24:259-328, 2015.
#'
#' S. Heinrich. Monte Carlo complexity of global solution of integral equations. \emph{Journal of Complexity}, 14(2):151-175, 1998.
#'
#' @param Lmin the minimum level of refinement.  Must be \eqn{\ge 2}.
#' @param Lmax the maximum level of refinement.  Must be \eqn{\ge} Lmin.
#' @param N0 initial number of samples which are used for the first 3 levels and
#'   for any subsequent levels which are automatically added.  Must be
#'   \eqn{> 0}.
#' @param eps the target accuracy of the estimate.  Must be \eqn{> 0}.
#' @param mlmc_l a user supplied function which provides the estimate for level
#'   l
#' @param alpha the weak error, \eqn{O(2^{-alpha*l})}.  If \code{NA} then
#'   \code{alpha} will be estimated.
#' @param beta the variance, \eqn{O(2^{-beta*l})}.  If \code{NA} then
#'   \code{beta} will be estimated.
#' @param gamma the sample cost, \eqn{O(2^{gamma*l})}.  Must be \eqn{> 0}.
#' @param parallel if an integer is supplied, R will fork \code{parallel} parallel
#'   processes an compute each level estimate in parallel.
#' @param ... additional arguments which are passed on when the user supplied
#'   \code{mlmc_l} function is called
#'
#' @return A list containing: \describe{
#'   \item{\code{P}}{The MLMC estimate;}
#'   \item{\code{Nl}}{A vector of the number of samples performed on each level.}
#' }
#'
#' @examples
#' mlmc(2, 6, 1000, 0.01, opre_l, gamma=1, option=1)
#'
#' mlmc(2, 10, 1000, 0.01, mcqmc06_l, gamma=1, option=1)
#'
#' @importFrom parallel mcmapply
#' @importFrom stats lm
#' @export
mlmc <- function(Lmin, Lmax, N0, eps, mlmc_l, alpha=NA, beta=NA, gamma, parallel=NA, ...) {
  # check parameters are acceptable
  if(Lmin<2) {
    stop("Lmin must be >= 2.")
  }
  if(Lmax<Lmin) {
    stop("must have Lmax >= Lmin.")
  }
  if(N0<=0 || eps<=0 || gamma <= 0){
    stop("N0, eps and gamma must all be greater than zero.")
  }
  if(!is.na(alpha) && alpha<=0) {
    stop("if specified, alpha must be greater than zero.  Set alpha to NA to automatically estimate.")
  }
  if(!is.na(beta) && beta<=0) {
    stop("if specified, beta must be greater than zero.  Set beta to NA to automatically estimate.")
  }

  # initialise the MLMC run
  est.alpha <- is.na(alpha)
  alpha <- ifelse(is.na(alpha), 0, alpha)
  est.beta  <- is.na(beta)
  beta <- ifelse(is.na(beta), 0, beta)

  theta <- 0.25

  L <- Lmin

  Nl <- rep(0, L+1)
  suml <- matrix(0, nrow=2, ncol=L+1)
  dNl <- rep(N0, L+1)

  while(sum(dNl) > 0) {
    # update sample sums from each level
    if(is.na(parallel)) {
      for(l in 0:L) {
        if(dNl[l+1] > 0) {
          sums        <- mlmc_l(l, dNl[l+1], ...)
          Nl[l+1]     <- Nl[l+1] + dNl[l+1]
          suml[1,l+1] <- suml[1,l+1] + sums[1]
          suml[2,l+1] <- suml[2,l+1] + sums[2]
        }
      }
    } else if(is.numeric(parallel)) {
      par.vars <- data.frame(l=0:L, dNl=dNl)[!(dNl==0),]
      sums <- mcmapply(function(l, dNl, ...) {
        mlmc_l(l, dNl, ...)
      }, l = par.vars$l, dNl = par.vars$dNl, ..., mc.preschedule = FALSE, mc.cores = parallel)
      Nl <- Nl+dNl
      suml[,!(dNl==0)] <- suml[,!(dNl==0)] + sums[1:2,]
    }

    # compute absolute average and variance
    ml <- abs(    suml[1,]/Nl)
    Vl <- pmax(0, suml[2,]/Nl - ml^2)

    # fix to cope with possible zero values for ml and Vl
    # (can happen in some applications when there are few samples)
    for(l in 3:(L+1)) {
      ml[l] <- max(ml[l], 0.5*ml[l-1]/2^alpha)
      Vl[l] <- max(Vl[l], 0.5*Vl[l-1]/2^beta)
    }

    # estimate alpha and beta by regression if specified
    if(est.alpha) {
      alpha <- max(0.5, -lm(y~x, data.frame(y=log2(ml[-1]), x=1:L))$coefficients["x"])
    }
    if(est.beta) {
      beta <- max(0.5, -lm(y~x, data.frame(y=log2(Vl[-1]), x=1:L))$coefficients["x"])
    }
    # choose optimal number of additional samples
    Cl  <- 2^(gamma*(0:L))
    Ns  <- ceiling(sqrt(Vl/Cl) * sum(sqrt(Vl*Cl)) / ((1-theta)*eps^2))
    dNl <- pmax(0, Ns-Nl)
    # if (almost) converged, estimate remaining error and decide
    # whether a new level is required
    if( !any(dNl > 0.01*Nl) ) {
      rem <- ml[L+1] / (2^alpha - 1)

      if(rem > sqrt(theta)*eps) {
        if(L==Lmax) {
          warning("failed to achieve weak convergence.")
        } else {
          L       <- L+1
          Vl[L+1] <- Vl[L] / 2^beta
          Nl[L+1] <- 0
          suml <- cbind(suml, 0)

          Cl <- 2^(gamma*(0:L))
          Ns  <- ceiling(sqrt(Vl/Cl) * sum(sqrt(Vl*Cl)) / ((1-theta)*eps^2))
          dNl <- pmax(0, Ns-Nl)
        }
      }
    }
  }

  # finally, evaluate multilevel estimator
  P <- sum(suml[1,]/Nl)
  list(P=P, Nl=Nl)
}



# Function to provide level l estimate based on parameters and net benefit function
EVPPI_l_p<-function(l,N, EVPPI_std_p = NULL)
{
  if(is.null(EVPPI_std_p)) return("EVPPI_std_p must be supplied")
  # this function generates all the random samples for the EVPPI calculation
  # then calculates the net benefit function
  
  # N is the number of outer samples, M=2^l is the number of inner samples
  sum1 <- rep(0, 7)
  # Need to store results of 2^0, 2^1,..., M= 2^l
  M = 2^(l + 1) 
  Np = max(M, 128)
  
  inputs = 1:ceiling(M * N / Np)
  
  # Iterate over the inputs
  # Can add parallel computation here
  results <- foreach(i=inputs,.export=c("n_samples", "n_treatment",
                                        "generate_net_benefit", "generate_input_parameters",
                                        'EVPPI_std_p'),.packages='MASS') %dorng% {
                                          NN=min(Np, N*M-(i-1)*Np)
                                          #####################################################
                                          EVPPI_std_p(M,NN/M)
                                          
                                          # if you want to calculate EVPPI for other parameters, you need to change this line using other corresponding function
                                          # and change the function name in the last element of the condition of foreach to the same one.
                                          ###########################################33########
                                        }   
  # Convert results of foreach to NetB matrix
  # If not interested in parallelisation could merge the foreach and for loop to simplify
  NetB = matrix(NA, M*N, n_treatment)
  for(i in inputs){
    NN <- min(Np, N*M-(i-1)*Np)
    nn = min(i*Np,N*M)
    NetB[(nn-NN+1):nn,] = matrix(unlist(results[i]),NN,n_treatment)
  }
  
  # Net benefit based on partial perfect information for every sample
  NetB_max = apply(NetB,1,max)
  # Expected value of max over each set of inner samples
  NetB_max_sample = apply(matrix(NetB_max,N,M,byrow=TRUE),1,mean)
  # Matrices needed for antithetic variable variance reduction
  NetB_low_c_1 = matrix(NA,N,n_treatment)
  NetB_low_c_2 = matrix(NA,N,n_treatment)
  NetB_low_f = matrix(NA,N,n_treatment)
  for(n in 1:n_treatment){
    # Net benefts for treatment n in matrix of size outer by inner samples
    temp = matrix(NetB[,n],N,M,byrow=TRUE)
    # Average net benefit over each set of inner samples
    NetB_low_f[,n] = apply(temp,1,mean)
    # Antithetic variable construction splits samples into first and second halves
    # Split formula into cases l=0 and l>0
    if(l==0){
      NetB_low_c_1[,n] = temp[,1]
      NetB_low_c_2[,n] = temp[,2]
    }else{
      NetB_low_c_1[,n] = apply(temp[,1:max(M/2,2)],1,mean)
      NetB_low_c_2[,n] = apply(temp[,min(M/2+1,M-1):M],1,mean)
    }
  }
  NetB_low_f_sample = apply(NetB_low_f,1,max)
  # Put antithetic variable construction together
  NetB_low_c_sample = (apply(NetB_low_c_1,1,max)+apply(NetB_low_c_2,1,max))/2
  Pf = NetB_max_sample - NetB_low_f_sample
  Pc = NetB_max_sample - NetB_low_c_sample
  
  sum1[1] = sum1[1] + sum(Pf-Pc);
  sum1[2] = sum1[2] + sum((Pf-Pc)^2);
  sum1[3] = sum1[3] + sum((Pf-Pc)^3);
  sum1[4] = sum1[4] + sum((Pf-Pc)^4);
  sum1[5] = sum1[5] + sum(Pf);
  sum1[6] = sum1[6] + sum(Pf^2);
  sum1[7] = sum1[7] + M*N;
  
  return(sum1)
}

