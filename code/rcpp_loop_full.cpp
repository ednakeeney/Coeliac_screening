#include <Rcpp.h>
using namespace Rcpp;

// You can source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


// [[Rcpp::export]]
// Arguments are the transition matrices in a column of n_cycles stacked matrices 
// and starting cohort vector used to find initial distribution
NumericMatrix rcpp_loop_full(NumericMatrix cohort_vectors_in, NumericMatrix transition_matrices,
                               int n_cycles, int n_tests, int n_samples, int n_states) {
  NumericMatrix cohort_vectors_out(n_cycles * n_tests * n_samples, 3 + n_states);
  // Initialise the output cohort vector
  cohort_vectors_out = cohort_vectors_in;
  NumericVector rm1, cm2;
  int current_cohort_row;
  NumericMatrix transition_matrices_temp(n_states, n_states);
  for(int i_test = 1; i_test <= n_tests; i_test++) {
    for(int i_sample = 1; i_sample <= n_samples; i_sample++) {
      for (int i_cycle = 1; i_cycle < n_cycles; i_cycle++) {
        // Cohort vector at previous cycle
        
        current_cohort_row = ((i_cycle - 1) * n_tests * n_samples) + 
          (i_test - 1) * (n_samples) +
          i_sample - 1;
        rm1 = cohort_vectors_out(Range(current_cohort_row, current_cohort_row), 
                                 Range(3, 3 + n_states - 1));

        transition_matrices_temp = transition_matrices(Range((i_cycle - 1) * (n_samples * n_states) +
          (i_sample - 1) * n_states,
          (i_cycle - 1) * (n_samples * n_states) +
            (i_sample - 1) * n_states + n_states - 1), 
            Range(3, 3 + n_states - 1));
        for (int j_state = 0; j_state < n_states; ++j_state) {
          //cm2 = transition_matrices(_,j_state);
          // Use the correct rows for this cycle
          cm2 = transition_matrices_temp(_, j_state);
          cohort_vectors_out((i_cycle) * (n_tests * n_samples) +
                                     (i_test - 1) * (n_samples) + 
                                     i_sample - 1,
                            3 + j_state) = std::inner_product(rm1.begin(), rm1.end(), cm2.begin(), 0.);              
        }
      }
    }
  }
  return(cohort_vectors_out);
}

