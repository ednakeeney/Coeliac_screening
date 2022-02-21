# Switched to using a data frame for transition matrices to ease implementation in C/C++

# Packages needed to easily manipulate data frames
require(reshape)
require(dplyr)


convert_cohort_vectors_to_df <- function(cohort_vectors) {
  
  # Convert transition matrices array to a dataframe
  # Starts in sample, cycle, from order but will be in cycle, sample, "from state", "to states" 
  # order to match C/C++
  # Also remove dimnames so data.frame can index using numbers (important if going to C/C++)
  cohort_vectors_temp <- cohort_vectors
  dimnames(cohort_vectors_temp) <- list(NULL, NULL, NULL, state_names)
  # List of data frames for transition from each state
  cohort_vectors_temp_df <- list()
  for(i_state in 1:n_states) {
    # Each stores the transition probabilities from i_state
    cohort_vectors_temp_df[[i_state]] <- cohort_vectors_temp[, , , i_state]
    # Convert the multidimensional array to a data frame
    cohort_vectors_temp_df[[i_state]] <- 
      melt(cohort_vectors_temp_df[[i_state]], varnames = c("test", "sample", "cycle"))
    # Name the state to which you're transiting
    colnames(cohort_vectors_temp_df[[i_state]])[4] <- state_names[i_state]
  }
  # Combine the data frames
  cohort_vectors_df <- do.call(cbind, cohort_vectors_temp_df)
  # Only keep unique columns
  cohort_vectors_df <- cohort_vectors_df[, unique(colnames(cohort_vectors_df))]
  # Sort 
  cohort_vectors_df  <- cohort_vectors_df %>% arrange(cycle, test, sample)
  # Ensure that first columns are cycle, sample, from
  cohort_vectors_df <- cohort_vectors_df[, c("cycle", "test", "sample", state_names)]
  
  return(cohort_vectors_df)
}

# Could test to see if they're the same
#cohort_vectors[1, 1, 1, ]
#cohort_vectors_df[which(with(cohort_vectors_df, test == 1 & sample == 1 & cycle ==1)), ]


