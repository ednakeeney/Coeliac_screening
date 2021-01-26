# Coeliac disease Markov model
# Edna Keeney
#setwd("C:/Users/ek14588/OneDrive - University of Bristol/Coeliac screening")
#'
#' @return Output
#' @export
markov_expanded <- function() {
  set.seed(14143)
  
  # Define the number and names of treatments

  n.treatments<-3
  t.names<-c("Test","Test + biopsy", "Double test")
  
  # Define the number and names of states of the model
  n.states<-13
  state.names<- c("CD GFD no complications",
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
  n.cycles<-50
  
  # Define simulation parameters
  # This is the number of PSA samples to use
  n.samples<- 100
  
  starting_age <- 30
  #############################################################################
  ## Input parameters #########################################################
  #############################################################################
  
  Adherence	<- 0.75
  Probability_late_diagnosis	<- 0.132
  
  #Initial cohort at diagnosis - will depend on age at diagnosis
  Probability_subfertility <- 0.0517
  Probability_osteoporosis <- 0.0125
  Probability_NHL <- 0.0008
  Probability_nocomplications <- 1 - (Probability_subfertility + Probability_osteoporosis + Probability_NHL)
  
  #On GFD
  Osteoporosis_probability_GFD_30	<- 0.0017
  Osteoporosis_probability_GFD_40	<- 0.0041
  Osteoporosis_probability_GFD_50 <- 	0.0098
  Osteoporosis_probability_GFD_60	<- 0.0178
  Osteoporosis_probability_GFD_70	<- 0.0286
  Osteoporosis_probability_GFD_80 <- 	0.0305
  Osteoporosis_probability_GFD_90 <- 	0.0573
  Subfertility_probability_GFD_30 <-	0.0017
  Subfertility_probability_GFD_40 <-	0.0002
  NHL_probability_GFD	<- 0.00072

  #Not on GFD
  Osteoporosis_probability_noGFD_30	<- 0.0017
  Osteoporosis_probability_noGFD_40	<- 0.0041
  Osteoporosis_probability_noGFD_50 <- 	0.0098
  Osteoporosis_probability_noGFD_60	<- 0.0178
  Osteoporosis_probability_noGFD_70	<- 0.0286
  Osteoporosis_probability_noGFD_80 <- 	0.0305
  Osteoporosis_probability_noGFD_90 <- 	0.0573
  Subfertility_probability_noGFD_30 <-	0.0017
  Subfertility_probability_noGFD_40 <-	0.0002
  NHL_probability_noGFD	<- 0.00072
  
  lifetables <- read.csv("lifetables.csv")
  percentage_male <- 0.5
  lifetables$Overall <- (percentage_male * lifetables$Males) + ((1-percentage_male) * lifetables$Females)
  lifetables30plus <- subset(lifetables, lifetables$Age > 29)
  
  Death_probability_nocomplications	<- lifetables30plus$Overall
 
  Mortality_NHL <- read.csv("NHL mortality.csv")
  
   Death_probability_NHL <-	Mortality_NHL$Mortality
  Death_probability_osteoporosis <-	0.0003
  
 
  # There is one transition matrix for each treatment option and each PSA sample
  # Store them in an array with (before filling in below) NA entries
  transition.matrices<-array(dim=c(n.samples, n.cycles, n.states,n.states),
                             dimnames=list(NULL,NULL, state.names,state.names))
  #change to list of arrays
  #different array for each cycle
  #transition matrices list
  #transition matrix 1 is an array
  
  # Now define the QALYS associated with the states per cycle
  # There is one for each PSA sample and each state
  # Store in an NA array and then fill in below
  state.qalys<-array(dim=c(n.samples, n.states),dimnames=list(NULL,state.names))
  
  # And finally define the state costs
  # There is one for each PSA sample and each state
  # Store in an NA array and then fill in below
  state.costs<-array(dim=c(n.samples, n.states),dimnames=list(NULL,state.names))
  
  
  # Define transition matrices, state utilities and costs 

  #CD GFD 
    transition.matrices[ , c(1:10),"CD GFD no complications", "CD GFD subfertility"]<- Subfertility_probability_GFD_30
    transition.matrices[ , c(11:20),"CD GFD no complications", "CD GFD subfertility"]<- Subfertility_probability_GFD_40
    transition.matrices[ , c(21:n.cycles),"CD GFD no complications", "CD GFD subfertility"]<- 0
    transition.matrices[ , c(1:10), "CD GFD no complications", "CD GFD osteoporosis"]<- Osteoporosis_probability_GFD_30
    transition.matrices[ , c(11:20), "CD GFD no complications", "CD GFD osteoporosis"]<- Osteoporosis_probability_GFD_40
    transition.matrices[ , c(21:30), "CD GFD no complications", "CD GFD osteoporosis"]<- Osteoporosis_probability_GFD_50
    transition.matrices[ , c(31:40), "CD GFD no complications", "CD GFD osteoporosis"]<- Osteoporosis_probability_GFD_60
    transition.matrices[ , c(41:50), "CD GFD no complications", "CD GFD osteoporosis"]<- Osteoporosis_probability_GFD_70
    transition.matrices[ , ,"CD GFD no complications", "CD GFD NHL"]<- NHL_probability_GFD
    for (i.age in 1:50){
    transition.matrices[, i.age, "CD GFD no complications", "Death"]<- Death_probability_nocomplications[i.age]
    transition.matrices[ , c(1:10), "CD GFD no complications", "CD GFD no complications"]<- 1- Subfertility_probability_GFD_30 - 
      Osteoporosis_probability_GFD_30 - NHL_probability_GFD - Death_probability_nocomplications[i.age]
    transition.matrices[ , c(11:20), "CD GFD no complications", "CD GFD no complications"]<- 1- Subfertility_probability_GFD_40 - 
      Osteoporosis_probability_GFD_40 - NHL_probability_GFD - Death_probability_nocomplications[i.age]
    transition.matrices[ , c(21:30), "CD GFD no complications", "CD GFD no complications"]<- 1- 
      Osteoporosis_probability_GFD_50 - NHL_probability_GFD - Death_probability_nocomplications[i.age]
    transition.matrices[ , c(31:40), "CD GFD no complications", "CD GFD no complications"]<- 1- 
      Osteoporosis_probability_GFD_60 - NHL_probability_GFD - Death_probability_nocomplications[i.age]
    transition.matrices[ , c(41:50), "CD GFD no complications", "CD GFD no complications"]<- 1- 
      Osteoporosis_probability_GFD_70 - NHL_probability_GFD - Death_probability_nocomplications[i.age]
    }
   

    transition.matrices[ , c(1:10),"CD GFD subfertility", "CD GFD osteoporosis"]<- Osteoporosis_probability_GFD_30
    transition.matrices[ , c(11:20), "CD GFD subfertility", "CD GFD osteoporosis"]<- Osteoporosis_probability_GFD_40
    transition.matrices[ , c(21:30), "CD GFD subfertility", "CD GFD osteoporosis"]<- Osteoporosis_probability_GFD_50
    transition.matrices[ , c(31:40), "CD GFD subfertility", "CD GFD osteoporosis"]<- Osteoporosis_probability_GFD_60
    transition.matrices[ , c(41:50), "CD GFD subfertility", "CD GFD osteoporosis"]<- Osteoporosis_probability_GFD_70
    transition.matrices[ , ,"CD GFD subfertility", "CD GFD NHL"]<- NHL_probability_GFD
    for (i.age in 1:50){
      transition.matrices[, i.age, "CD GFD subfertility", "Death"]<- Death_probability_nocomplications[i.age]
      transition.matrices[ , c(1:10), "CD GFD subfertility", "CD GFD subfertility"]<- 1- 
        Osteoporosis_probability_GFD_30 - NHL_probability_GFD - Death_probability_nocomplications[i.age]
      transition.matrices[ , c(11:20), "CD GFD subfertility", "CD GFD subfertility"]<- 1- 
        Osteoporosis_probability_GFD_40 - NHL_probability_GFD - Death_probability_nocomplications[i.age]
      transition.matrices[ , c(21:30), "CD GFD subfertility", "CD GFD subfertility"]<- 1- 
        Osteoporosis_probability_GFD_50 - NHL_probability_GFD - Death_probability_nocomplications[i.age]
      transition.matrices[ , c(31:40), "CD GFD subfertility", "CD GFD subfertility"]<- 1- 
        Osteoporosis_probability_GFD_60 - NHL_probability_GFD - Death_probability_nocomplications[i.age]
      transition.matrices[ , c(41:50), "CD GFD subfertility", "CD GFD subfertility"]<- 1- 
        Osteoporosis_probability_GFD_70 - NHL_probability_GFD - Death_probability_nocomplications[i.age]
    }
  
    
    transition.matrices[ , ,"CD GFD osteoporosis", "CD GFD NHL"]<- NHL_probability_GFD
    transition.matrices[, , "CD GFD osteoporosis", "Death"]<- Death_probability_osteoporosis
      transition.matrices[ , , "CD GFD osteoporosis", "CD GFD osteoporosis"]<- 1- 
        NHL_probability_GFD - Death_probability_osteoporosis
    
    
   
    for (i.cycle in 1:10){
    transition.matrices[ , i.cycle, "CD GFD NHL", "Death"]<-Death_probability_NHL[i.cycle]
    transition.matrices[, i.cycle,"CD GFD NHL", "CD GFD NHL"]<- 1- Death_probability_NHL[i.cycle]
    }
    transition.matrices[ , 11:n.cycles, "CD GFD NHL", "Death"]<-Death_probability_NHL[10]
    transition.matrices[, 11:n.cycles,"CD GFD NHL", "CD GFD NHL"]<- 1- Death_probability_NHL[10]
    

#CD no GFD    
    transition.matrices[ , c(1:10),"CD no GFD no complications", "CD no GFD subfertility"]<- Subfertility_probability_noGFD_30
    transition.matrices[ , c(11:20),"CD no GFD no complications", "CD no GFD subfertility"]<- Subfertility_probability_noGFD_40
    transition.matrices[ , c(21:n.cycles),"CD no GFD no complications", "CD no GFD subfertility"]<- 0
    transition.matrices[ , c(1:10), "CD no GFD no complications", "CD no GFD osteoporosis"]<- Osteoporosis_probability_noGFD_30
    transition.matrices[ , c(11:20), "CD no GFD no complications", "CD no GFD osteoporosis"]<- Osteoporosis_probability_noGFD_40
    transition.matrices[ , c(21:30), "CD no GFD no complications", "CD no GFD osteoporosis"]<- Osteoporosis_probability_noGFD_50
    transition.matrices[ , c(31:40), "CD no GFD no complications", "CD no GFD osteoporosis"]<- Osteoporosis_probability_noGFD_60
    transition.matrices[ , c(41:50), "CD no GFD no complications", "CD no GFD osteoporosis"]<- Osteoporosis_probability_noGFD_70
    transition.matrices[ , ,"CD no GFD no complications", "CD no GFD NHL"]<- NHL_probability_noGFD
    for (i.age in 1:50){
      transition.matrices[, i.age, "CD no GFD no complications", "Death"]<- Death_probability_nocomplications[i.age]
      transition.matrices[ , c(1:10), "CD no GFD no complications", "CD no GFD no complications"]<- 1- Subfertility_probability_noGFD_30 - 
        Osteoporosis_probability_noGFD_30 - NHL_probability_noGFD - Death_probability_nocomplications[i.age]
      transition.matrices[ , c(11:20), "CD no GFD no complications", "CD no GFD no complications"]<- 1- Subfertility_probability_noGFD_40 - 
        Osteoporosis_probability_noGFD_40 - NHL_probability_noGFD - Death_probability_nocomplications[i.age]
      transition.matrices[ , c(21:30), "CD no GFD no complications", "CD no GFD no complications"]<- 1- 
        Osteoporosis_probability_noGFD_50 - NHL_probability_noGFD - Death_probability_nocomplications[i.age]
      transition.matrices[ , c(31:40), "CD no GFD no complications", "CD no GFD no complications"]<- 1- 
        Osteoporosis_probability_noGFD_60 - NHL_probability_noGFD - Death_probability_nocomplications[i.age]
      transition.matrices[ , c(41:50), "CD no GFD no complications", "CD no GFD no complications"]<- 1- 
        Osteoporosis_probability_noGFD_70 - NHL_probability_noGFD - Death_probability_nocomplications[i.age]
    }
    
    
    transition.matrices[ , c(1:10),"CD no GFD subfertility", "CD no GFD osteoporosis"]<- Osteoporosis_probability_noGFD_30
    transition.matrices[ , c(11:20), "CD no GFD subfertility", "CD no GFD osteoporosis"]<- Osteoporosis_probability_noGFD_40
    transition.matrices[ , c(21:30), "CD no GFD subfertility", "CD no GFD osteoporosis"]<- Osteoporosis_probability_noGFD_50
    transition.matrices[ , c(31:40), "CD no GFD subfertility", "CD no GFD osteoporosis"]<- Osteoporosis_probability_noGFD_60
    transition.matrices[ , c(41:50), "CD no GFD subfertility", "CD no GFD osteoporosis"]<- Osteoporosis_probability_noGFD_70
    transition.matrices[ , ,"CD no GFD subfertility", "CD no GFD NHL"]<- NHL_probability_noGFD
    for (i.age in 1:50){
      transition.matrices[, i.age, "CD no GFD subfertility", "Death"]<- Death_probability_nocomplications[i.age]
      transition.matrices[ , c(1:10), "CD no GFD subfertility", "CD no GFD subfertility"]<- 1- 
        Osteoporosis_probability_noGFD_30 - NHL_probability_noGFD - Death_probability_nocomplications[i.age]
      transition.matrices[ , c(11:20), "CD no GFD subfertility", "CD no GFD subfertility"]<- 1- 
        Osteoporosis_probability_noGFD_40 - NHL_probability_noGFD - Death_probability_nocomplications[i.age]
      transition.matrices[ , c(21:30), "CD no GFD subfertility", "CD no GFD subfertility"]<- 1- 
        Osteoporosis_probability_noGFD_50 - NHL_probability_noGFD - Death_probability_nocomplications[i.age]
      transition.matrices[ , c(31:40), "CD no GFD subfertility", "CD no GFD subfertility"]<- 1- 
        Osteoporosis_probability_noGFD_60 - NHL_probability_noGFD - Death_probability_nocomplications[i.age]
      transition.matrices[ , c(41:50), "CD no GFD subfertility", "CD no GFD subfertility"]<- 1- 
        Osteoporosis_probability_noGFD_70 - NHL_probability_noGFD - Death_probability_nocomplications[i.age]
    }
    
    
    transition.matrices[ , ,"CD no GFD osteoporosis", "CD no GFD NHL"]<- NHL_probability_noGFD
    transition.matrices[, , "CD no GFD osteoporosis", "Death"]<- Death_probability_osteoporosis
    transition.matrices[ , , "CD no GFD osteoporosis", "CD no GFD osteoporosis"]<- 1- 
      NHL_probability_noGFD - Death_probability_osteoporosis
    

for (i.cycle in 1:10){
  transition.matrices[ , i.cycle, "CD no GFD NHL", "Death"]<-Death_probability_NHL[i.cycle]
  transition.matrices[, i.cycle,"CD no GFD NHL", "CD no GFD NHL"]<- 1- Death_probability_NHL[i.cycle]
}
transition.matrices[ , 11:n.cycles, "CD no GFD NHL", "Death"]<-Death_probability_NHL[10]
transition.matrices[, 11:n.cycles,"CD no GFD NHL", "CD no GFD NHL"]<- 1- Death_probability_NHL[10]


#Undiagnosed CD  
    transition.matrices[ , ,"Undiagnosed CD no complications", "CD GFD no complications"]<- Adherence * Probability_late_diagnosis
    transition.matrices[ , ,"Undiagnosed CD no complications", "CD no GFD no complications"]<- (1-Adherence) * Probability_late_diagnosis
   
     transition.matrices[ , c(1:10),"Undiagnosed CD no complications", "Undiagnosed CD subfertility"]<- Subfertility_probability_noGFD_30
    transition.matrices[ , c(11:20),"Undiagnosed CD no complications", "Undiagnosed CD subfertility"]<- Subfertility_probability_noGFD_40
    transition.matrices[ , c(21:n.cycles),"Undiagnosed CD no complications", "Undiagnosed CD subfertility"]<- 0
    transition.matrices[ , c(1:10), "Undiagnosed CD no complications", "Undiagnosed CD osteoporosis"]<- Osteoporosis_probability_noGFD_30
    transition.matrices[ , c(11:20), "Undiagnosed CD no complications", "Undiagnosed CD osteoporosis"]<- Osteoporosis_probability_noGFD_40
    transition.matrices[ , c(21:30), "Undiagnosed CD no complications", "Undiagnosed CD osteoporosis"]<- Osteoporosis_probability_noGFD_50
    transition.matrices[ , c(31:40), "Undiagnosed CD no complications", "Undiagnosed CD osteoporosis"]<- Osteoporosis_probability_noGFD_60
    transition.matrices[ , c(41:50), "Undiagnosed CD no complications", "Undiagnosed CD osteoporosis"]<- Osteoporosis_probability_noGFD_70
    transition.matrices[ , ,"Undiagnosed CD no complications", "Undiagnosed CD NHL"]<- NHL_probability_noGFD
    for (i.age in 1:50){
      transition.matrices[, i.age, "Undiagnosed CD no complications", "Death"]<- Death_probability_nocomplications[i.age]
      transition.matrices[ , c(1:10), "Undiagnosed CD no complications", "Undiagnosed CD no complications"]<- 1- (Adherence * Probability_late_diagnosis) - 
        ((1-Adherence) * Probability_late_diagnosis) - Subfertility_probability_noGFD_30 - 
        Osteoporosis_probability_noGFD_30 - NHL_probability_noGFD - Death_probability_nocomplications[i.age]
      transition.matrices[ , c(11:20), "Undiagnosed CD no complications", "Undiagnosed CD no complications"]<- 1- (Adherence * Probability_late_diagnosis) - 
        ((1-Adherence) * Probability_late_diagnosis) - Subfertility_probability_noGFD_40 - 
        Osteoporosis_probability_noGFD_40 - NHL_probability_noGFD - Death_probability_nocomplications[i.age]
      transition.matrices[ , c(21:30), "Undiagnosed CD no complications", "Undiagnosed CD no complications"]<- 1- (Adherence * Probability_late_diagnosis) - 
        ((1-Adherence) * Probability_late_diagnosis) - 
        Osteoporosis_probability_noGFD_50 - NHL_probability_noGFD - Death_probability_nocomplications[i.age]
      transition.matrices[ , c(31:40), "Undiagnosed CD no complications", "Undiagnosed CD no complications"]<- 1- (Adherence * Probability_late_diagnosis) - 
        ((1-Adherence) * Probability_late_diagnosis) -
        Osteoporosis_probability_noGFD_60 - NHL_probability_noGFD - Death_probability_nocomplications[i.age]
      transition.matrices[ , c(41:50), "Undiagnosed CD no complications", "Undiagnosed CD no complications"]<- 1- (Adherence * Probability_late_diagnosis) - 
        ((1-Adherence) * Probability_late_diagnosis) -
        Osteoporosis_probability_noGFD_70 - NHL_probability_noGFD - Death_probability_nocomplications[i.age]
    }
    
    transition.matrices[ , ,"Undiagnosed CD subfertility", "CD GFD subfertility"]<- Adherence * Probability_late_diagnosis
    transition.matrices[ , ,"Undiagnosed CD subfertility", "CD no GFD subfertility"]<- (1-Adherence) * Probability_late_diagnosis
    transition.matrices[ , c(1:10),"Undiagnosed CD subfertility", "Undiagnosed CD osteoporosis"]<- Osteoporosis_probability_noGFD_30
    transition.matrices[ , c(11:20), "Undiagnosed CD subfertility", "Undiagnosed CD osteoporosis"]<- Osteoporosis_probability_noGFD_40
    transition.matrices[ , c(21:30), "Undiagnosed CD subfertility", "Undiagnosed CD osteoporosis"]<- Osteoporosis_probability_noGFD_50
    transition.matrices[ , c(31:40), "Undiagnosed CD subfertility", "Undiagnosed CD osteoporosis"]<- Osteoporosis_probability_noGFD_60
    transition.matrices[ , c(41:50), "Undiagnosed CD subfertility", "Undiagnosed CD osteoporosis"]<- Osteoporosis_probability_noGFD_70
    transition.matrices[ , ,"Undiagnosed CD subfertility", "Undiagnosed CD NHL"]<- NHL_probability_noGFD
    for (i.age in 1:50){
      transition.matrices[, i.age, "Undiagnosed CD subfertility", "Death"]<- Death_probability_nocomplications[i.age]
      transition.matrices[ , c(1:10), "Undiagnosed CD subfertility", "Undiagnosed CD subfertility"]<- 1- (Adherence * Probability_late_diagnosis) - 
        ((1-Adherence) * Probability_late_diagnosis)  - 
        Osteoporosis_probability_noGFD_30 - NHL_probability_noGFD - Death_probability_nocomplications[i.age]
      transition.matrices[ , c(11:20), "Undiagnosed CD subfertility", "Undiagnosed CD subfertility"]<- 1- (Adherence * Probability_late_diagnosis) - 
        ((1-Adherence) * Probability_late_diagnosis) -
        Osteoporosis_probability_noGFD_40 - NHL_probability_noGFD - Death_probability_nocomplications[i.age]
      transition.matrices[ , c(21:30), "Undiagnosed CD subfertility", "Undiagnosed CD subfertility"]<- 1- (Adherence * Probability_late_diagnosis) - 
        ((1-Adherence) * Probability_late_diagnosis) -
        Osteoporosis_probability_noGFD_50 - NHL_probability_noGFD - Death_probability_nocomplications[i.age]
      transition.matrices[ , c(31:40), "Undiagnosed CD subfertility", "Undiagnosed CD subfertility"]<- 1- (Adherence * Probability_late_diagnosis) - 
        ((1-Adherence) * Probability_late_diagnosis) -
        Osteoporosis_probability_noGFD_60 - NHL_probability_noGFD - Death_probability_nocomplications[i.age]
      transition.matrices[ , c(41:50), "Undiagnosed CD subfertility", "Undiagnosed CD subfertility"]<- 1- (Adherence * Probability_late_diagnosis) - 
        ((1-Adherence) * Probability_late_diagnosis) -
        Osteoporosis_probability_noGFD_70 - NHL_probability_noGFD - Death_probability_nocomplications[i.age]
    }
    
    
    transition.matrices[ , ,"Undiagnosed CD osteoporosis", "CD GFD osteoporosis"]<- Adherence * Probability_late_diagnosis
    transition.matrices[ , ,"Undiagnosed CD osteoporosis", "CD no GFD osteoporosis"]<- (1-Adherence) * Probability_late_diagnosis
    transition.matrices[ , ,"Undiagnosed CD osteoporosis", "Undiagnosed CD NHL"]<- NHL_probability_noGFD
    transition.matrices[ , ,"Undiagnosed CD osteoporosis", "Death"]<- Death_probability_osteoporosis
    transition.matrices[ , ,"Undiagnosed CD osteoporosis", "Undiagnosed CD osteoporosis"]<- 1 -  (Adherence * Probability_late_diagnosis) - 
      ((1-Adherence) * Probability_late_diagnosis) - Death_probability_osteoporosis - NHL_probability_noGFD
    
    transition.matrices[ , ,"Undiagnosed CD NHL", "CD GFD NHL"]<- Adherence * Probability_late_diagnosis
    transition.matrices[ , ,"Undiagnosed CD NHL", "CD no GFD NHL"]<- (1-Adherence) * Probability_late_diagnosis
    for (i.cycle in 1:10){
      transition.matrices[ , i.cycle, "Undiagnosed CD NHL", "Death"]<-Death_probability_NHL[i.cycle]
      transition.matrices[, i.cycle,"Undiagnosed CD NHL", "Undiagnosed CD NHL"]<- 1- (Adherence * Probability_late_diagnosis) - 
        ((1-Adherence) * Probability_late_diagnosis) - Death_probability_NHL[i.cycle]
    }
    transition.matrices[ , 11:n.cycles, "Undiagnosed CD NHL", "Death"]<-Death_probability_NHL[10]
    transition.matrices[, 11:n.cycles,"Undiagnosed CD NHL", "Undiagnosed CD NHL"]<- 1- Death_probability_NHL[10]
    
    # State utilities
    # Anything between 0 and 1

    state.qalys[, "CD GFD no complications"]<- 0.85
    state.qalys[, "CD GFD subfertility"] <- 0.85 - 0.07
    state.qalys[,"CD GFD osteoporosis"] <- 0.7
    state.qalys[,"CD GFD NHL"] <- 0.618
    state.qalys[,"CD no GFD no complications"] <- 0.71
    state.qalys[,"CD no GFD subfertility"] <- 0.71 - 0.07
    state.qalys[, "CD no GFD osteoporosis"] <- 0.5 
    state.qalys[,"CD no GFD NHL"] <- 0.618
    state.qalys[,"Undiagnosed CD no complications"] <- 0.65
    state.qalys[,"Undiagnosed CD subfertility"] <- 0.65 - 0.07
    state.qalys[,"Undiagnosed CD osteoporosis"] <- 0.5 
    state.qalys[,"Undiagnosed CD NHL"] <- 0.618
    state.qalys[,"Death"] <- 0
    
    
    # State costs
    # Assumed normal with sd small enough to avoid negative values
    state.costs[, "CD GFD no complications"] <- 50
    state.costs[, "CD GFD subfertility"] <- 60
    state.costs[,"CD GFD osteoporosis"] <- 10
      state.costs[,"CD GFD NHL"] <- 30
    state.costs[,"CD no GFD no complications"] <- 40
    state.costs[,"CD no GFD subfertility"] <- 70
    state.costs[, "CD no GFD osteoporosis"] <-80
      state.costs[,"CD no GFD NHL"] <- 90
    state.costs[,"Undiagnosed CD no complications"] <- 100
    state.costs[,"Undiagnosed CD subfertility"] <- 110
    state.costs[,"Undiagnosed CD osteoporosis"] <- 120
      state.costs[,"Undiagnosed CD NHL"] <- 30
    state.costs[,"Death"] <- 0
    
  #How to add costs and disutilities for anemia?
  
  
  # Define the treatment costs
  # One for each PSA sample and each treatment
  # Treatment costs are actually fixed but this allows flexibility if we
  # want to include uncertainty/randomness in the cost
  treatment.costs<-array(dim=c(n.treatments,n.samples),dimnames=list(t.names,NULL))
  
  # Cost of the smoking cessation website is a one-off subscription fee of ?50
  treatment.costs["Test",]<-50
  # Zero cost for standard of care
  treatment.costs["Test + biopsy",]<-1000
  treatment.costs["Double test",]<-500
  
  #############################################################################
  ## Simulation ###############################################################
  #############################################################################
  
  # Build an array to store the cohort vector at each cycle
  # Each cohort vector has n.states elements: probability of being in each state,
  # There is one cohort vector for each treatment, for each PSA sample, for each cycle.
  cohort.vectors<-array(dim=c(n.treatments,n.samples,n.cycles,n.states),
                        dimnames=list(t.names,NULL,NULL,state.names))
  
  # This will be related to decision tree accuracy
  cohort.vectors[1,,1,"CD GFD no complications"] <- tp[,1] * Probability_nocomplications * Adherence
  cohort.vectors[2,,1,"CD GFD no complications"] <- tp[,2] * Probability_nocomplications * Adherence
  cohort.vectors[3,,1,"CD GFD no complications"] <- tp[,3] * Probability_nocomplications * Adherence
  cohort.vectors[1,,1,"CD GFD subfertility"] <- tp[,1] * Probability_subfertility * Adherence
  cohort.vectors[2,,1,"CD GFD subfertility"] <- tp[,2] * Probability_subfertility * Adherence
  cohort.vectors[3,,1,"CD GFD subfertility"] <- tp[,3] * Probability_subfertility * Adherence
  cohort.vectors[1,,1,"CD GFD osteoporosis"] <- tp[,1] * Probability_osteoporosis * Adherence
  cohort.vectors[2,,1,"CD GFD osteoporosis"] <- tp[,2] * Probability_osteoporosis * Adherence
  cohort.vectors[3,,1,"CD GFD osteoporosis"] <- tp[,3] * Probability_osteoporosis * Adherence
  cohort.vectors[1,,1,"CD GFD NHL"] <- tp[,1] * Probability_NHL * Adherence
  cohort.vectors[2,,1,"CD GFD NHL"] <- tp[,2] * Probability_NHL * Adherence
  cohort.vectors[3,,1,"CD GFD NHL"] <- tp[,3] * Probability_NHL * Adherence
  
  cohort.vectors[1,,1,"CD no GFD no complications"] <- tp[,1] * Probability_nocomplications * (1-Adherence)
  cohort.vectors[2,,1,"CD no GFD no complications"] <- tp[,2] * Probability_nocomplications * (1-Adherence)
  cohort.vectors[3,,1,"CD no GFD no complications"] <- tp[,3] * Probability_nocomplications * (1-Adherence)
  cohort.vectors[1,,1,"CD no GFD subfertility"] <- tp[,1] * Probability_subfertility * (1-Adherence)
  cohort.vectors[2,,1,"CD no GFD subfertility"] <- tp[,2] * Probability_subfertility * (1-Adherence)
  cohort.vectors[3,,1,"CD no GFD subfertility"] <- tp[,3] * Probability_subfertility * (1-Adherence)
  cohort.vectors[1,,1,"CD no GFD osteoporosis"] <- tp[,1] * Probability_osteoporosis * (1-Adherence)
  cohort.vectors[2,,1,"CD no GFD osteoporosis"] <- tp[,2] * Probability_osteoporosis * (1-Adherence)
  cohort.vectors[3,,1,"CD no GFD osteoporosis"] <- tp[,3] * Probability_osteoporosis * (1-Adherence)
  cohort.vectors[1,,1,"CD no GFD NHL"] <- tp[,1] * Probability_NHL * (1-Adherence)
  cohort.vectors[2,,1,"CD no GFD NHL"] <- tp[,2] * Probability_NHL * (1-Adherence)
  cohort.vectors[3,,1,"CD no GFD NHL"] <- tp[,3] * Probability_NHL * (1-Adherence)
  
  cohort.vectors[1,,1,"Undiagnosed CD no complications"] <- fn[,1] * Probability_nocomplications 
  cohort.vectors[2,,1,"Undiagnosed CD no complications"] <- fn[,2] * Probability_nocomplications 
  cohort.vectors[3,,1,"Undiagnosed CD no complications"] <- fn[,3] * Probability_nocomplications 
  cohort.vectors[1,,1,"Undiagnosed CD subfertility"] <- fn[,1] * Probability_subfertility 
  cohort.vectors[2,,1,"Undiagnosed CD subfertility"] <- fn[,2] * Probability_subfertility 
  cohort.vectors[3,,1,"Undiagnosed CD subfertility"] <- fn[,3] * Probability_subfertility 
  cohort.vectors[1,,1,"Undiagnosed CD osteoporosis"] <- fn[,1] * Probability_osteoporosis 
  cohort.vectors[2,,1,"Undiagnosed CD osteoporosis"] <- fn[,2] * Probability_osteoporosis 
  cohort.vectors[3,,1,"Undiagnosed CD osteoporosis"] <- fn[,3] * Probability_osteoporosis 
  cohort.vectors[1,,1,"Undiagnosed CD NHL"] <- fn[,1] * Probability_NHL 
  cohort.vectors[2,,1,"Undiagnosed CD NHL"] <- fn[,2] * Probability_NHL 
  cohort.vectors[3,,1,"Undiagnosed CD NHL"] <- fn[,3] * Probability_NHL 
  
  # Build an array to store the costs and QALYs accrued per cycle
  # One for each treatment, for each PSA sample, for each cycle
  # These will be filled in below in the main model code
  # Then discounted and summed to contribute to total costs and total QALYs
  cycle.costs<-array(dim=c(n.treatments,n.samples,n.cycles),
                     dimnames=list(t.names,NULL,NULL))
  cycle.qalys<-array(dim=c(n.treatments,n.samples,n.cycles),
                     dimnames=list(t.names,NULL,NULL))
  
  # Build arrays to store the total costs and total QALYs
  # There is one for each treatment and each PSA sample
  # These are filled in below using cycle.costs, 
  # treatment.costs, and cycle.qalys
  total.costs<-array(dim=c(n.treatments,n.samples),
                     dimnames=list(t.names,NULL))
  total.qalys<-array(dim=c(n.treatments,n.samples),
                     dimnames=list(t.names,NULL))
  
  #i.treatment <- 1
  #i.sample <- 1
  #i.cycle <- 2
  
  disc_vec <- (1/1.035)^rep(c(0:(n.cycles/2-1)),each=2)
  
  # The remainder of the cohort.vectors will be filled in by Markov updating below
  
  # Main model code
  # Loop over the treatment options
  

    
    # Loop over the PSA samples
    for(i.sample in 1:n.samples)
    {
      
      transition.matrices_sample <- transition.matrices[i.sample,,,]
      # Loop over the cycles
      # Cycle 1 is already defined so only need to update cycles 2:n.cycles
      for(i.cycle in 2:n.cycles)
      {
        # Markov update
        # Multiply previous cycle's cohort vector by transition matrix
        # i.e. pi_j = pi_(j-1)*P
        cohort.vectors[i.treatment, i.sample,i.cycle,]<-
          cohort.vectors[i.treatment, i.sample,i.cycle-1,] %*%
          transition.matrices_sample
      }
      
      cohort.vectors_tr_sample <- cohort.vectors[i.treatment,i.sample,,]
      
      # Now use the cohort vectors to calculate the 
      # total costs for each cycle
      cycle.costs[i.treatment,i.sample,]<-
        cohort.vectors_tr_sample%*%state.costs[i.sample,]
      # And total QALYs for each cycle
      cycle.qalys[i.treatment,i.sample,]<-
        cohort.vectors_tr_sample%*%state.qalys[i.sample,]
      
      # Combine the cycle.costs and treatment.costs to get total costs
      # Apply the discount factor 
      # (1 in first year, 1.035 in second, 1.035^2 in third, and so on)
      # Each year acounts for two cycles so need to repeat the discount values
      total.costs[i.treatment,i.sample]<-treatment.costs[i.treatment,i.sample] +
        cycle.costs[i.treatment,i.sample,]%*%
        disc_vec
      
      # Combine the cycle.qalys to get total qalys
      # Apply the discount factor 
      # (1 in first year, 1.035 in second, 1.035^2 in third, and so on)
      # Each year acounts for two cycles so need to repeat the discount values
      total.qalys[i.treatment,i.sample]<-cycle.qalys[i.treatment,i.sample,]%*%
        disc_vec
    }
  
  
  #############################################################################
  ## Analysis of results ######################################################
  #############################################################################
  output <- list()
  # Average costs
  # These are ?50 on the website and 0 on standard of care as there are no
  # costs other than the website subscription cost
  output$average.costs<-rowMeans(total.costs)
  # Average effects (in QALY units)
  # These are slightly higher on the website as higher probability of 
  # quitting smoking
  output$average.effects<-rowMeans(total.qalys)
  
  # Incremental costs and effects relative to standard of care
  # No uncertainty in the costs as the website cost is fixed at ?50
  output$incremental.costs<-total.costs["SoC with website",]-total.costs["SoC",]
  # In some samples the website leads to higher QALYs but in others it is negative
  # There is uncertainty as to whether the website is an improvement over SoC
  output$incremental.effects<-total.qalys["SoC with website",]-total.qalys["SoC",]
  
  # The ICER comparing Standard of care with website to standard of care
  # This is much lower than the ?20,000 willingness-to-pay threshold indicating
  # good value for money
  output$ICER<-mean(output$incremental.costs)/mean(output$incremental.effects)
  
  # Incremental net benefit at the ?20,000 willingness-to-pay
  # Sometimes positive (website more cost-effective) and sometimes negative (SoC more cost-effective)
  # Need to look at averages and consider probabilities of cost-effectiveness
  output$incremental.net.benefit<-20000*output$incremental.effects-output$incremental.costs
  
  # Average incremental net benefit
  # This is positive indicating cost-effectiveness at the ?20,000 threshold
  output$average.inb<-mean(output$incremental.net.benefit)
  
  # Probability cost-effective
  # This is the proportion of samples for which the incremental net benefit is positive
  # It is clost to 72%, representing good degree of certainty
  # in recommendation to adopt the smoking cessation website
  output$probability.cost.effective<-sum(output$incremental.net.benefit>0)/n.samples
  
  # Now use the BCEA package to analyse the results...
  output
}

