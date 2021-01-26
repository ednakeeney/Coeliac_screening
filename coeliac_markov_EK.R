# Coeliac disease Markov model
# Edna Keeney
setwd("C:/Users/ek14588/OneDrive - University of Bristol/Coeliac screening")
#'
#' @return Output
#' @export
#markov_expanded <- function() {
  set.seed(14143)
  
  # Define the number and names of treatments

  n_treatments <- 3
  t_names <- c("Test", "Test + biopsy", "Double test")
  
  # Define the number and names of states of the model
  n_states <- 13
  state_names <- c("CD GFD no complications",
                  "CD GFD subfertility",
                  "CD GFD osteoporosis",
                  "CD GFD NHL",
                  "CD no GFD no complications",
                  "CD no GFD subfertility",
                  "CD no GFD osteoporosis",
                  "CD no GFD NHL",
                  "Undiagnosed CD no complications",
                  "Undiagnosed CD subfertility",
                  "Undiagnosed CD osteoporosis",
                  "Undiagnosed CD NHL",
                  "Death")
 
  
  
   # Define the number of cycles
  n_cycles <- 50
  
  # Define simulation parameters
  # This is the number of PSA samples to use
  n_samples <- 100
  
  starting_age <- 30
  #############################################################################
  ## Input parameters #########################################################
  #############################################################################
  
  adherence	<- 0.75
  probability_late_diagnosis <- 0.132
  
  #Initial cohort at diagnosis - will depend on age at diagnosis
  probability_subfertility <- 0.0517
  probability_osteoporosis <- 0.0125
  probability_NHL <- 0.0008
  probability_nocomplications <- 1 - (probability_subfertility + probability_osteoporosis + probability_NHL)
  
  #On GFD
  osteoporosis_probability_GFD_30	<- 0.0017
  osteoporosis_probability_GFD_40	<- 0.0041
  osteoporosis_probability_GFD_50 <- 	0.0098
  osteoporosis_probability_GFD_60	<- 0.0178
  osteoporosis_probability_GFD_70	<- 0.0286
  osteoporosis_probability_GFD_80 <- 	0.0305
  osteoporosis_probability_GFD_90 <- 	0.0573
  subfertility_probability_GFD_30 <-	0.0017
  subfertility_probability_GFD_40 <-	0.0002
  NHL_probability_GFD	<- 0.00072

  #Not on GFD
  osteoporosis_probability_noGFD_30	<- 0.0017
  osteoporosis_probability_noGFD_40	<- 0.0041
  osteoporosis_probability_noGFD_50 <- 	0.0098
  osteoporosis_probability_noGFD_60	<- 0.0178
  osteoporosis_probability_noGFD_70	<- 0.0286
  osteoporosis_probability_noGFD_80 <- 	0.0305
  osteoporosis_probability_noGFD_90 <- 	0.0573
  subfertility_probability_noGFD_30 <-	0.0017
  subfertility_probability_noGFD_40 <-	0.0002
  NHL_probability_noGFD	<- 0.00072
  
  lifetables <- read.csv("lifetables.csv")
  percentage_male <- 0.5
  lifetables$Overall <- (percentage_male * lifetables$Males) + ((1-percentage_male) * lifetables$Females)
  lifetables30plus <- subset(lifetables, lifetables$Age > 29)
  death_probability_nocomplications	<- lifetables30plus$Overall
 
  mortality_NHL <- read.csv("NHL mortality.csv")
  
   death_probability_NHL <-	mortality_NHL$Mortality
  death_probability_osteoporosis <-	0.0003
  
 
  # There is one transition matrix for each cycle and each PSA sample
  # Store them in an array with (before filling in below) NA entries
  transition_matrices <- array(dim=c(n_samples, n_cycles, n_states, n_states), #should n_cycles or n_samples come first?
                             dimnames=list(NULL, NULL, state_names, state_names))
  #change to list of arrays
  #different array for each cycle
  #transition matrices list
  #transition matrix 1 is an array
  
  # Now define the QALYS associated with the states per cycle
  # There is one for each PSA sample and each state
  # Store in an NA array and then fill in below
  state_qalys <- array(dim=c(n_samples, n_states), dimnames = list(NULL, state_names))
  
  # And finally define the state costs
  # There is one for each PSA sample and each state
  # Store in an NA array and then fill in below
  state_costs <- array(dim=c(n_samples, n_states), dimnames = list(NULL, state_names))
  
  
  # Define transition matrices, state utilities and costs 

  #CD GFD 
    transition_matrices[, c(1:10), "CD GFD no complications", "CD GFD subfertility"] <- subfertility_probability_GFD_30
    transition_matrices[, c(11:20), "CD GFD no complications", "CD GFD subfertility"] <- subfertility_probability_GFD_40
    transition_matrices[, c(21:n_cycles), "CD GFD no complications", "CD GFD subfertility"] <- 0
    transition_matrices[, c(1:10), "CD GFD no complications", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_30
    transition_matrices[, c(11:20), "CD GFD no complications", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_40
    transition_matrices[, c(21:30), "CD GFD no complications", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_50
    transition_matrices[, c(31:40), "CD GFD no complications", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_60
    transition_matrices[, c(41:50), "CD GFD no complications", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_70
    transition_matrices[, ,"CD GFD no complications", "CD GFD NHL"] <- NHL_probability_GFD
    for (i_age in 1:50){
    transition_matrices[, i_age, "CD GFD no complications", "Death"] <- Death_probability_nocomplications[i_age]
    transition_matrices[, c(1:10), "CD GFD no complications", "CD GFD no complications"] <- 1 - subfertility_probability_GFD_30 - 
      osteoporosis_probability_GFD_30 - NHL_probability_GFD - death_probability_nocomplications[i_age]
    transition_matrices[, c(11:20), "CD GFD no complications", "CD GFD no complications"] <- 1 - subfertility_probability_GFD_40 - 
      osteoporosis_probability_GFD_40 - NHL_probability_GFD - death_probability_nocomplications[i_age]
    transition_matrices[, c(21:30), "CD GFD no complications", "CD GFD no complications"] <- 1 - 
      osteoporosis_probability_GFD_50 - NHL_probability_GFD - death_probability_nocomplications[i_age]
    transition_matrices[, c(31:40), "CD GFD no complications", "CD GFD no complications"] <- 1 - 
      osteoporosis_probability_GFD_60 - NHL_probability_GFD - death_probability_nocomplications[i_age]
    transition_matrices[, c(41:50), "CD GFD no complications", "CD GFD no complications"] <- 1 - 
      osteoporosis_probability_GFD_70 - NHL_probability_GFD - death_probability_nocomplications[i_age]
    }
   

    transition_matrices[, c(1:10),"CD GFD subfertility", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_30
    transition_matrices[, c(11:20), "CD GFD subfertility", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_40
    transition_matrices[, c(21:30), "CD GFD subfertility", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_50
    transition_matrices[, c(31:40), "CD GFD subfertility", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_60
    transition_matrices[, c(41:50), "CD GFD subfertility", "CD GFD osteoporosis"] <- osteoporosis_probability_GFD_70
    transition_matrices[, ,"CD GFD subfertility", "CD GFD NHL"] <- NHL_probability_GFD
    for (i_age in 1:50){
      transition_matrices[, i_age, "CD GFD subfertility", "Death"] <- death_probability_nocomplications[i_age]
      transition_matrices[, c(1:10), "CD GFD subfertility", "CD GFD subfertility"] <- 1 - 
        osteoporosis_probability_GFD_30 - NHL_probability_GFD - death_probability_nocomplications[i_age]
      transition_matrices[, c(11:20), "CD GFD subfertility", "CD GFD subfertility"] <- 1 - 
        osteoporosis_probability_GFD_40 - NHL_probability_GFD - death_probability_nocomplications[i_age]
      transition_matrices[, c(21:30), "CD GFD subfertility", "CD GFD subfertility"] <- 1 - 
        osteoporosis_probability_GFD_50 - NHL_probability_GFD - death_probability_nocomplications[i_age]
      transition_matrices[, c(31:40), "CD GFD subfertility", "CD GFD subfertility"] <- 1 - 
        osteoporosis_probability_GFD_60 - NHL_probability_GFD - death_probability_nocomplications[i_age]
      transition_matrices[, c(41:50), "CD GFD subfertility", "CD GFD subfertility"] <- 1 - 
        osteoporosis_probability_GFD_70 - NHL_probability_GFD - death_probability_nocomplications[i_age]
    }
  
    
    transition_matrices[, ,"CD GFD osteoporosis", "CD GFD NHL"] <- NHL_probability_GFD
    transition_matrices[, , "CD GFD osteoporosis", "Death"] <- death_probability_osteoporosis
      transition_matrices[, , "CD GFD osteoporosis", "CD GFD osteoporosis"] <- 1 - 
        NHL_probability_GFD - death_probability_osteoporosis
    
    
   
    for (i_cycle in 1:10){
    transition_matrices[, i_cycle, "CD GFD NHL", "Death"] <- death_probability_NHL[i_cycle]
    transition_matrices[, i_cycle,"CD GFD NHL", "CD GFD NHL"] <- 1 - death_probability_NHL[i_cycle]
    }
    transition_matrices[, 11:n_cycles, "CD GFD NHL", "Death"] <- death_probability_NHL[10]
    transition_matrices[, 11:n_cycles,"CD GFD NHL", "CD GFD NHL"] <- 1 - death_probability_NHL[10]
    

#CD no GFD    
    transition_matrices[, c(1:10),"CD no GFD no complications", "CD no GFD subfertility"] <- subfertility_probability_noGFD_30
    transition_matrices[, c(11:20),"CD no GFD no complications", "CD no GFD subfertility"] <- subfertility_probability_noGFD_40
    transition_matrices[, c(21:n_cycles),"CD no GFD no complications", "CD no GFD subfertility"] <- 0
    transition_matrices[, c(1:10), "CD no GFD no complications", "CD no GFD osteoporosis"] <- osteoporosis_probability_noGFD_30
    transition_matrices[, c(11:20), "CD no GFD no complications", "CD no GFD osteoporosis"] <- osteoporosis_probability_noGFD_40
    transition_matrices[, c(21:30), "CD no GFD no complications", "CD no GFD osteoporosis"] <- osteoporosis_probability_noGFD_50
    transition_matrices[, c(31:40), "CD no GFD no complications", "CD no GFD osteoporosis"] <- osteoporosis_probability_noGFD_60
    transition_matrices[, c(41:50), "CD no GFD no complications", "CD no GFD osteoporosis"] <- osteoporosis_probability_noGFD_70
    transition_matrices[, ,"CD no GFD no complications", "CD no GFD NHL"] <- NHL_probability_noGFD
    for (i_age in 1:50){
      transition_matrices[, i_age, "CD no GFD no complications", "Death"] <- death_probability_nocomplications[i_age]
      transition_matrices[, c(1:10), "CD no GFD no complications", "CD no GFD no complications"] <- 1 - subfertility_probability_noGFD_30 - 
        osteoporosis_probability_noGFD_30 - NHL_probability_noGFD - death_probability_nocomplications[i_age]
      transition_matrices[, c(11:20), "CD no GFD no complications", "CD no GFD no complications"] <- 1 - subfertility_probability_noGFD_40 - 
        osteoporosis_probability_noGFD_40 - NHL_probability_noGFD - death_probability_nocomplications[i_age]
      transition_matrices[, c(21:30), "CD no GFD no complications", "CD no GFD no complications"] <- 1 - 
        osteoporosis_probability_noGFD_50 - NHL_probability_noGFD - death_probability_nocomplications[i_age]
      transition_matrices[, c(31:40), "CD no GFD no complications", "CD no GFD no complications"] <- 1 - 
        osteoporosis_probability_noGFD_60 - NHL_probability_noGFD - death_probability_nocomplications[i_age]
      transition_matrices[, c(41:50), "CD no GFD no complications", "CD no GFD no complications"] <- 1 - 
        osteoporosis_probability_noGFD_70 - NHL_probability_noGFD - death_probability_nocomplications[i_age]
    }
    
    
    transition_matrices[, c(1:10),"CD no GFD subfertility", "CD no GFD osteoporosis"] <- osteoporosis_probability_noGFD_30
    transition_matrices[, c(11:20), "CD no GFD subfertility", "CD no GFD osteoporosis"] <- osteoporosis_probability_noGFD_40
    transition_matrices[, c(21:30), "CD no GFD subfertility", "CD no GFD osteoporosis"] <- osteoporosis_probability_noGFD_50
    transition_matrices[, c(31:40), "CD no GFD subfertility", "CD no GFD osteoporosis"] <- osteoporosis_probability_noGFD_60
    transition_matrices[, c(41:50), "CD no GFD subfertility", "CD no GFD osteoporosis"] <- osteoporosis_probability_noGFD_70
    transition_matrices[, ,"CD no GFD subfertility", "CD no GFD NHL"] <- NHL_probability_noGFD
    for (i_age in 1:50){
      transition_matrices[, i_age, "CD no GFD subfertility", "Death"] <- death_probability_nocomplications[i_age]
      transition_matrices[, c(1:10), "CD no GFD subfertility", "CD no GFD subfertility"] <- 1- 
        osteoporosis_probability_noGFD_30 - NHL_probability_noGFD - death_probability_nocomplications[i_age]
      transition_matrices[, c(11:20), "CD no GFD subfertility", "CD no GFD subfertility"] <- 1- 
        osteoporosis_probability_noGFD_40 - NHL_probability_noGFD - death_probability_nocomplications[i_age]
      transition_matrices[, c(21:30), "CD no GFD subfertility", "CD no GFD subfertility"] <- 1- 
        osteoporosis_probability_noGFD_50 - NHL_probability_noGFD - death_probability_nocomplications[i_age]
      transition_matrices[, c(31:40), "CD no GFD subfertility", "CD no GFD subfertility"] <- 1- 
        osteoporosis_probability_noGFD_60 - NHL_probability_noGFD - death_probability_nocomplications[i_age]
      transition_matrices[, c(41:50), "CD no GFD subfertility", "CD no GFD subfertility"] <- 1- 
        osteoporosis_probability_noGFD_70 - NHL_probability_noGFD - death_probability_nocomplications[i_age]
    }
    
    
    transition_matrices[, ,"CD no GFD osteoporosis", "CD no GFD NHL"] <- NHL_probability_noGFD
    transition_matrices[, , "CD no GFD osteoporosis", "Death"] <- death_probability_osteoporosis
    transition_matrices[, , "CD no GFD osteoporosis", "CD no GFD osteoporosis"] <- 1- 
      NHL_probability_noGFD - death_probability_osteoporosis
    

for (i_cycle in 1:10){
  transition_matrices[, i_cycle, "CD no GFD NHL", "Death"] <-death_probability_NHL[i_cycle]
  transition_matrices[, i_cycle,"CD no GFD NHL", "CD no GFD NHL"] <- 1- death_probability_NHL[i_cycle]
}
transition_matrices[, 11:n_cycles, "CD no GFD NHL", "Death"] <-death_probability_NHL[10]
transition_matrices[, 11:n_cycles,"CD no GFD NHL", "CD no GFD NHL"] <- 1- death_probability_NHL[10]


#Undiagnosed CD  
    transition_matrices[, ,"Undiagnosed CD no complications", "CD GFD no complications"] <- adherence * probability_late_diagnosis
    transition_matrices[, ,"Undiagnosed CD no complications", "CD no GFD no complications"] <- (1-adherence) * probability_late_diagnosis
   
     transition_matrices[, c(1:10),"Undiagnosed CD no complications", "Undiagnosed CD subfertility"] <- subfertility_probability_noGFD_30
    transition_matrices[, c(11:20),"Undiagnosed CD no complications", "Undiagnosed CD subfertility"] <- subfertility_probability_noGFD_40
    transition_matrices[, c(21:n_cycles),"Undiagnosed CD no complications", "Undiagnosed CD subfertility"] <- 0
    transition_matrices[, c(1:10), "Undiagnosed CD no complications", "Undiagnosed CD osteoporosis"] <- osteoporosis_probability_noGFD_30
    transition_matrices[, c(11:20), "Undiagnosed CD no complications", "Undiagnosed CD osteoporosis"] <- osteoporosis_probability_noGFD_40
    transition_matrices[, c(21:30), "Undiagnosed CD no complications", "Undiagnosed CD osteoporosis"] <- osteoporosis_probability_noGFD_50
    transition_matrices[, c(31:40), "Undiagnosed CD no complications", "Undiagnosed CD osteoporosis"] <- osteoporosis_probability_noGFD_60
    transition_matrices[, c(41:50), "Undiagnosed CD no complications", "Undiagnosed CD osteoporosis"] <- osteoporosis_probability_noGFD_70
    transition_matrices[, ,"Undiagnosed CD no complications", "Undiagnosed CD NHL"] <- NHL_probability_noGFD
    for (i_age in 1:50){
      transition_matrices[, i_age, "Undiagnosed CD no complications", "Death"] <- death_probability_nocomplications[i_age]
      transition_matrices[, c(1:10), "Undiagnosed CD no complications", "Undiagnosed CD no complications"] <- 1- (adherence * probability_late_diagnosis) - 
        ((1-adherence) * probability_late_diagnosis) - subfertility_probability_noGFD_30 - 
        osteoporosis_probability_noGFD_30 - NHL_probability_noGFD - death_probability_nocomplications[i_age]
      transition_matrices[, c(11:20), "Undiagnosed CD no complications", "Undiagnosed CD no complications"] <- 1- (adherence * probability_late_diagnosis) - 
        ((1-adherence) * probability_late_diagnosis) - subfertility_probability_noGFD_40 - 
        osteoporosis_probability_noGFD_40 - NHL_probability_noGFD - death_probability_nocomplications[i_age]
      transition_matrices[, c(21:30), "Undiagnosed CD no complications", "Undiagnosed CD no complications"] <- 1- (adherence * probability_late_diagnosis) - 
        ((1-adherence) * probability_late_diagnosis) - 
        osteoporosis_probability_noGFD_50 - NHL_probability_noGFD - death_probability_nocomplications[i_age]
      transition_matrices[, c(31:40), "Undiagnosed CD no complications", "Undiagnosed CD no complications"] <- 1- (adherence * probability_late_diagnosis) - 
        ((1-adherence) * probability_late_diagnosis) -
        osteoporosis_probability_noGFD_60 - NHL_probability_noGFD - death_probability_nocomplications[i_age]
      transition_matrices[, c(41:50), "Undiagnosed CD no complications", "Undiagnosed CD no complications"] <- 1- (adherence * probability_late_diagnosis) - 
        ((1-adherence) * probability_late_diagnosis) -
        osteoporosis_probability_noGFD_70 - NHL_probability_noGFD - death_probability_nocomplications[i_age]
    }
    
    transition_matrices[, ,"Undiagnosed CD subfertility", "CD GFD subfertility"] <- adherence * probability_late_diagnosis
    transition_matrices[, ,"Undiagnosed CD subfertility", "CD no GFD subfertility"] <- (1-adherence) * probability_late_diagnosis
    transition_matrices[, c(1:10),"Undiagnosed CD subfertility", "Undiagnosed CD osteoporosis"] <- osteoporosis_probability_noGFD_30
    transition_matrices[, c(11:20), "Undiagnosed CD subfertility", "Undiagnosed CD osteoporosis"] <- osteoporosis_probability_noGFD_40
    transition_matrices[, c(21:30), "Undiagnosed CD subfertility", "Undiagnosed CD osteoporosis"] <- osteoporosis_probability_noGFD_50
    transition_matrices[, c(31:40), "Undiagnosed CD subfertility", "Undiagnosed CD osteoporosis"] <- osteoporosis_probability_noGFD_60
    transition_matrices[, c(41:50), "Undiagnosed CD subfertility", "Undiagnosed CD osteoporosis"] <- osteoporosis_probability_noGFD_70
    transition_matrices[, ,"Undiagnosed CD subfertility", "Undiagnosed CD NHL"] <- NHL_probability_noGFD
    for (i_age in 1:50){
      transition_matrices[, i_age, "Undiagnosed CD subfertility", "Death"] <- death_probability_nocomplications[i_age]
      transition_matrices[, c(1:10), "Undiagnosed CD subfertility", "Undiagnosed CD subfertility"] <- 1- (adherence * probability_late_diagnosis) - 
        ((1-adherence) * probability_late_diagnosis)  - 
        osteoporosis_probability_noGFD_30 - NHL_probability_noGFD - death_probability_nocomplications[i_age]
      transition_matrices[, c(11:20), "Undiagnosed CD subfertility", "Undiagnosed CD subfertility"] <- 1- (adherence * probability_late_diagnosis) - 
        ((1-adherence) * probability_late_diagnosis) -
        osteoporosis_probability_noGFD_40 - NHL_probability_noGFD - death_probability_nocomplications[i_age]
      transition_matrices[, c(21:30), "Undiagnosed CD subfertility", "Undiagnosed CD subfertility"] <- 1- (adherence * probability_late_diagnosis) - 
        ((1-adherence) * probability_late_diagnosis) -
        osteoporosis_probability_noGFD_50 - NHL_probability_noGFD - death_probability_nocomplications[i_age]
      transition_matrices[, c(31:40), "Undiagnosed CD subfertility", "Undiagnosed CD subfertility"] <- 1- (adherence * probability_late_diagnosis) - 
        ((1-adherence) * probability_late_diagnosis) -
        osteoporosis_probability_noGFD_60 - NHL_probability_noGFD - death_probability_nocomplications[i_age]
      transition_matrices[, c(41:50), "Undiagnosed CD subfertility", "Undiagnosed CD subfertility"] <- 1- (adherence * probability_late_diagnosis) - 
        ((1-adherence) * probability_late_diagnosis) -
        osteoporosis_probability_noGFD_70 - NHL_probability_noGFD - death_probability_nocomplications[i_age]
    }
    
    
    transition_matrices[, ,"Undiagnosed CD osteoporosis", "CD GFD osteoporosis"] <- adherence * probability_late_diagnosis
    transition_matrices[, ,"Undiagnosed CD osteoporosis", "CD no GFD osteoporosis"] <- (1-adherence) * probability_late_diagnosis
    transition_matrices[, ,"Undiagnosed CD osteoporosis", "Undiagnosed CD NHL"] <- NHL_probability_noGFD
    transition_matrices[, ,"Undiagnosed CD osteoporosis", "Death"] <- death_probability_osteoporosis
    transition_matrices[, ,"Undiagnosed CD osteoporosis", "Undiagnosed CD osteoporosis"] <- 1 -  (adherence * probability_late_diagnosis) - 
      ((1-adherence) * probability_late_diagnosis) - death_probability_osteoporosis - NHL_probability_noGFD
    
    transition_matrices[, ,"Undiagnosed CD NHL", "CD GFD NHL"] <- adherence * probability_late_diagnosis
    transition_matrices[, ,"Undiagnosed CD NHL", "CD no GFD NHL"] <- (1-adherence) * probability_late_diagnosis
    for (i_cycle in 1:10){
      transition_matrices[, i_cycle, "Undiagnosed CD NHL", "Death"] <-death_probability_NHL[i_cycle]
      transition_matrices[, i_cycle,"Undiagnosed CD NHL", "Undiagnosed CD NHL"] <- 1- (adherence * probability_late_diagnosis) - 
        ((1-adherence) * probability_late_diagnosis) - death_probability_NHL[i_cycle]
    }
    transition_matrices[, 11:n_cycles, "Undiagnosed CD NHL", "Death"] <-death_probability_NHL[10]
    transition_matrices[, 11:n_cycles,"Undiagnosed CD NHL", "Undiagnosed CD NHL"] <- 1- death_probability_NHL[10]
    
    # State utilities
    # Anything between 0 and 1

    state_qalys[, "CD GFD no complications"] <- 0.85
    state_qalys[, "CD GFD subfertility"] <- 0.85 - 0.07
    state_qalys[,"CD GFD osteoporosis"] <- 0.7
    state_qalys[,"CD GFD NHL"] <- 0.618
    state_qalys[,"CD no GFD no complications"] <- 0.71
    state_qalys[,"CD no GFD subfertility"] <- 0.71 - 0.07
    state_qalys[, "CD no GFD osteoporosis"] <- 0.5 
    state_qalys[,"CD no GFD NHL"] <- 0.618
    state_qalys[,"Undiagnosed CD no complications"] <- 0.65
    state_qalys[,"Undiagnosed CD subfertility"] <- 0.65 - 0.07
    state_qalys[,"Undiagnosed CD osteoporosis"] <- 0.5 
    state_qalys[,"Undiagnosed CD NHL"] <- 0.618
    state_qalys[,"Death"] <- 0
    
    
    # State costs
    # Assumed normal with sd small enough to avoid negative values
    state_costs[, "CD GFD no complications"] <- 50
    state_costs[, "CD GFD subfertility"] <- 60
    state_costs[,"CD GFD osteoporosis"] <- 10
      state_costs[,"CD GFD NHL"] <- 30
    state_costs[,"CD no GFD no complications"] <- 40
    state_costs[,"CD no GFD subfertility"] <- 70
    state_costs[, "CD no GFD osteoporosis"] <-80
      state_costs[,"CD no GFD NHL"] <- 90
    state_costs[,"Undiagnosed CD no complications"] <- 100
    state_costs[,"Undiagnosed CD subfertility"] <- 110
    state_costs[,"Undiagnosed CD osteoporosis"] <- 120
      state_costs[,"Undiagnosed CD NHL"] <- 30
    state_costs[,"Death"] <- 0
    
  #How to add costs and disutilities for anemia?
  
  
  # Define the treatment costs
  # One for each PSA sample and each treatment
  # Treatment costs are actually fixed but this allows flexibility if we
  # want to include uncertainty/randomness in the cost
  treatment_costs<-array(dim=c(n_treatments,n_samples),dimnames=list(t_names,NULL))
  
  # Cost of the smoking cessation website is a one-off subscription fee of ?50
  treatment_costs["Test", ] <-50
  # Zero cost for standard of care
  treatment_costs["Test + biopsy", ] <-1000
  treatment_costs["Double test", ] <-500
  
  #############################################################################
  ## Simulation ###############################################################
  #############################################################################
  
  # Build an array to store the cohort vector at each cycle
  # Each cohort vector has n_states elements: probability of being in each state,
  # There is one cohort vector for each treatment, for each PSA sample, for each cycle.
  cohort_vectors<-array(dim=c(n_treatments,n_samples,n_cycles,n_states),  #do I need to add another n_states?
                        dimnames=list(t_names,NULL,NULL,state_names, state_names))
  
  # This will be related to decision tree accuracy
  cohort_vectors[1, , 1,"CD GFD no complications"] <- tp[, 1] * probability_nocomplications * adherence
  cohort_vectors[2, , 1,"CD GFD no complications"] <- tp[, 2] * probability_nocomplications * adherence
  cohort_vectors[3, , 1,"CD GFD no complications"] <- tp[, 3] * probability_nocomplications * adherence
  cohort_vectors[1, , 1,"CD GFD subfertility"] <- tp[, 1] * probability_subfertility * adherence
  cohort_vectors[2, , 1,"CD GFD subfertility"] <- tp[, 2] * probability_subfertility * adherence
  cohort_vectors[3, , 1,"CD GFD subfertility"] <- tp[, 3] * probability_subfertility * adherence
  cohort_vectors[1, , 1,"CD GFD osteoporosis"] <- tp[, 1] * probability_osteoporosis * adherence
  cohort_vectors[2, , 1,"CD GFD osteoporosis"] <- tp[, 2] * probability_osteoporosis * adherence
  cohort_vectors[3, , 1,"CD GFD osteoporosis"] <- tp[, 3] * probability_osteoporosis * adherence
  cohort_vectors[1, , 1,"CD GFD NHL"] <- tp[, 1] * probability_NHL * adherence
  cohort_vectors[2, , 1,"CD GFD NHL"] <- tp[, 2] * probability_NHL * adherence
  cohort_vectors[3, , 1,"CD GFD NHL"] <- tp[, 3] * probability_NHL * adherence
  
  cohort_vectors[1, , 1,"CD no GFD no complications"] <- tp[, 1] * probability_nocomplications * (1-adherence)
  cohort_vectors[2, , 1,"CD no GFD no complications"] <- tp[, 2] * probability_nocomplications * (1-adherence)
  cohort_vectors[3, , 1,"CD no GFD no complications"] <- tp[, 3] * probability_nocomplications * (1-adherence)
  cohort_vectors[1, , 1,"CD no GFD subfertility"] <- tp[, 1] * probability_subfertility * (1-adherence)
  cohort_vectors[2, , 1,"CD no GFD subfertility"] <- tp[, 2] * probability_subfertility * (1-adherence)
  cohort_vectors[3, , 1,"CD no GFD subfertility"] <- tp[, 3] * probability_subfertility * (1-adherence)
  cohort_vectors[1, , 1,"CD no GFD osteoporosis"] <- tp[, 1] * probability_osteoporosis * (1-adherence)
  cohort_vectors[2, , 1,"CD no GFD osteoporosis"] <- tp[, 2] * probability_osteoporosis * (1-adherence)
  cohort_vectors[3, , 1,"CD no GFD osteoporosis"] <- tp[, 3] * probability_osteoporosis * (1-adherence)
  cohort_vectors[1, , 1,"CD no GFD NHL"] <- tp[, 1] * probability_NHL * (1-adherence)
  cohort_vectors[2, , 1,"CD no GFD NHL"] <- tp[, 2] * probability_NHL * (1-adherence)
  cohort_vectors[3, , 1,"CD no GFD NHL"] <- tp[, 3] * probability_NHL * (1-adherence)
  
  cohort_vectors[1, , 1,"Undiagnosed CD no complications"] <- fn[, 1] * probability_nocomplications 
  cohort_vectors[2, , 1,"Undiagnosed CD no complications"] <- fn[, 2] * probability_nocomplications 
  cohort_vectors[3, , 1,"Undiagnosed CD no complications"] <- fn[, 3] * probability_nocomplications 
  cohort_vectors[1, , 1,"Undiagnosed CD subfertility"] <- fn[, 1] * probability_subfertility 
  cohort_vectors[2, , 1,"Undiagnosed CD subfertility"] <- fn[, 2] * probability_subfertility 
  cohort_vectors[3, , 1,"Undiagnosed CD subfertility"] <- fn[, 3] * probability_subfertility 
  cohort_vectors[1, , 1,"Undiagnosed CD osteoporosis"] <- fn[, 1] * probability_osteoporosis 
  cohort_vectors[2, , 1,"Undiagnosed CD osteoporosis"] <- fn[, 2] * probability_osteoporosis 
  cohort_vectors[3, , 1,"Undiagnosed CD osteoporosis"] <- fn[, 3] * probability_osteoporosis 
  cohort_vectors[1, , 1,"Undiagnosed CD NHL"] <- fn[, 1] * probability_NHL 
  cohort_vectors[2, , 1,"Undiagnosed CD NHL"] <- fn[, 2] * probability_NHL 
  cohort_vectors[3, , 1,"Undiagnosed CD NHL"] <- fn[, 3] * probability_NHL 
  
  # Build an array to store the costs and QALYs accrued per cycle
  # One for each treatment, for each PSA sample, for each cycle
  # These will be filled in below in the main model code
  # Then discounted and summed to contribute to total costs and total QALYs
  cycle_costs <- array(dim = c(n_treatments, n_samples, n_cycles),
                     dimnames = list(t_names, NULL, NULL))
  cycle_qalys <- array(dim = c(n_treatments, n_samples, n_cycles),
                     dimnames = list(t_names, NULL, NULL))
  
  # Build arrays to store the total costs and total QALYs
  # There is one for each treatment and each PSA sample
  # These are filled in below using cycle_costs, 
  # treatment_costs, and cycle_qalys
  total_costs <- array(dim = c(n_treatments, n_samples),
                     dimnames = list(t_names, NULL))
  total_qalys <- array(dim = c(n_treatments, n_samples),
                     dimnames = list(t_names, NULL))
  
  #i_treatment <- 1
  #i_sample <- 1
  #i_cycle <- 2
  
  disc_vec <- (1/1.035) ^ rep(c(0:(n_cycles/2-1)), each = 2)
  
  # The remainder of the cohort_vectors will be filled in by Markov updating below
  
  # Main model code
  # Loop over the treatment options
  

    
    # Loop over the PSA samples
    for(i_sample in 1:n_samples)
    {
      
      transition_matrices_sample <- transition_matrices[i_sample, , , ]
      # Loop over the cycles
      # Cycle 1 is already defined so only need to update cycles 2:n_cycles
      for(i_cycle in 2:n_cycles)
      {
        # Markov update
        # Multiply previous cycle's cohort vector by transition matrix
        # i.e. pi_j = pi_(j-1)*P
        cohort_vectors[i_treatment, i_sample, i_cycle, ] <-
          cohort_vectors[i_treatment, i_sample, i_cycle-1, ] %*%
          transition_matrices_sample
      }
      
      cohort_vectors_tr_sample <- cohort_vectors[i_treatment, i_sample, , ]
      
      # Now use the cohort vectors to calculate the 
      # total costs for each cycle
      cycle_costs[i_treatment, i_sample, ] <-
        cohort_vectors_tr_sample %*% state_costs[i_sample, ]
      # And total QALYs for each cycle
      cycle_qalys[i_treatment, i_sample, ] <-
        cohort_vectors_tr_sample %*% state_qalys[i_sample, ]
      
      # Combine the cycle_costs and treatment_costs to get total costs
      # Apply the discount factor 
      # (1 in first year, 1.035 in second, 1.035^2 in third, and so on)
      # Each year acounts for two cycles so need to repeat the discount values
      total_costs[i_treatment, i_sample] <- treatment_costs[i_treatment, i_sample] +
        cycle_costs[i_treatment, i_sample, ] %*%
        disc_vec
      
      # Combine the cycle_qalys to get total qalys
      # Apply the discount factor 
      # (1 in first year, 1.035 in second, 1.035^2 in third, and so on)
      # Each year acounts for two cycles so need to repeat the discount values
      total_qalys[i_treatment, i_sample] <- cycle_qalys[i_treatment, i_sample, ] %*%
        disc_vec
    }
  
  
  #############################################################################
  ## Analysis of results ######################################################
  #############################################################################
  output <- list()
  # Average costs
  # These are ?50 on the website and 0 on standard of care as there are no
  # costs other than the website subscription cost
  output$average_costs <- rowMeans(total_costs)
  # Average effects (in QALY units)
  # These are slightly higher on the website as higher probability of 
  # quitting smoking
  output$average_effects <- rowMeans(total_qalys)
  
  # Incremental costs and effects relative to standard of care
  # No uncertainty in the costs as the website cost is fixed at ?50
  output$incremental_costs <- total_costs["SoC with website", ] - total_costs["SoC", ]
  # In some samples the website leads to higher QALYs but in others it is negative
  # There is uncertainty as to whether the website is an improvement over SoC
  output$incremental_effects <- total_qalys["SoC with website", ] - total_qalys["SoC", ]
  
  # The ICER comparing Standard of care with website to standard of care
  # This is much lower than the ?20,000 willingness-to-pay threshold indicating
  # good value for money
  output$ICER <- mean(output$incremental_costs)/mean(output$incremental_effects)
  
  # Incremental net benefit at the ?20,000 willingness-to-pay
  # Sometimes positive (website more cost-effective) and sometimes negative (SoC more cost-effective)
  # Need to look at averages and consider probabilities of cost-effectiveness
  output$incremental_net_benefit <- 20000*output$incremental_effects - output$incremental_costs
  
  # Average incremental net benefit
  # This is positive indicating cost-effectiveness at the ?20,000 threshold
  output$average_inb <- mean(output$incremental_net_benefit)
  
  # Probability cost-effective
  # This is the proportion of samples for which the incremental net benefit is positive
  # It is clost to 72%, representing good degree of certainty
  # in recommendation to adopt the smoking cessation website
  output$probability_cost_effective <- sum(output$incremental_net_benefit>0)/n_samples
  
  # Now use the BCEA package to analyse the results___
  output
}

