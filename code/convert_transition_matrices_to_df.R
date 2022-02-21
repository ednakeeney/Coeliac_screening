# Switched to using a data frame for transition matrices to ease implementation in C/C++

# Packages needed to easily manipulate data frames
require(reshape)
require(dplyr)


convert_transition_matrices_to_df <- function(transition_matrices) {
  
  # Convert transition matrices array to a dataframe
  # Starts in sample, cycle, from order but will be in cycle, sample, "from state", "to states" 
  # order to match C/C++
  # Also remove dimnames so data.frame can index using numbers (important if going to C/C++)
  transition_matrices_temp <- transition_matrices
  dimnames(transition_matrices_temp) <- list(NULL, NULL, NULL, state_names)
  # List of data frames for transition from each state
  transition_matrices_temp_df <- list()
  for(i_state in 1:n_states) {
    # Each stores the transition probabilities from i_state
    transition_matrices_temp_df[[i_state]] <- transition_matrices_temp[, , , i_state]
    # Convert the multidimensional array to a data frame
    transition_matrices_temp_df[[i_state]] <- 
      melt(transition_matrices_temp_df[[i_state]], varnames = c("sample", "cycle", "from"))
    # Name the state to which you're transiting
    colnames(transition_matrices_temp_df[[i_state]])[4] <- state_names[i_state]
  }
  # Combine the data frames
  transition_matrices_df <- do.call(cbind, transition_matrices_temp_df)
  # Only keep unique columns
  transition_matrices_df <- transition_matrices_df[, unique(colnames(transition_matrices_df))]
  # Sort 
  transition_matrices_df  <- transition_matrices_df %>% arrange(cycle, sample, from)
  # Ensure that first columns are cycle, sample, from
  transition_matrices_df <- transition_matrices_df[, c("cycle", "sample", "from", state_names)]
  
  return(transition_matrices_df)
}

# Compare transition matrices
#i_sample <- 10; i_cycle <- 50
#transition_matrices[i_sample, i_cycle, , ]
#transition_matrices_df[with(transition_matrices_df,  cycle == i_cycle & sample == i_sample), ]
#transition_matrices_df[transition_matrices_df$cycle == i_cycle & transition_matrices_df$sample == i_sample, ]



