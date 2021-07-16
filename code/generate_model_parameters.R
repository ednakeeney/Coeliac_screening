generate_model_parameters <- function(starting_age) {
  
 starting_age <- ifelse(population == "adults", 50, 10) #based on mean age in under and over 18s in CPRD cost data
 
 # Define the number of cycles
 n_cycles <- 90 - starting_age
 
  duration_of_symptoms_adults <- 10.93  #Violato et al 
  duration_of_symptoms_adults_sd <- 13.10
  duration_of_symptoms_adults_location <- log(duration_of_symptoms_adults ^ 2 / sqrt(duration_of_symptoms_adults_sd ^ 2 + duration_of_symptoms_adults ^ 2))
  duration_of_symptoms_adults_shape <- sqrt(log(1 + (duration_of_symptoms_adults_sd ^ 2 / duration_of_symptoms_adults ^ 2)))
  duration_of_symptoms_adults <- rlnorm(n = n_samples, duration_of_symptoms_adults_location,  duration_of_symptoms_adults_shape)     #calculated from Violato et al 2019
  rate_of_symptoms_adults <- 1 / duration_of_symptoms_adults
  probability_late_diagnosis_adults <- 1 - exp(-rate_of_symptoms_adults)
  
  duration_of_symptoms_children <- 3.34  #Violato et al
  duration_of_symptoms_children_sd <- 3.71
  duration_of_symptoms_children_location <- log(duration_of_symptoms_children ^ 2 / sqrt(duration_of_symptoms_children_sd ^ 2 + duration_of_symptoms_children ^ 2))
  duration_of_symptoms_children_shape <- sqrt(log(1 + (duration_of_symptoms_children_sd ^ 2 / duration_of_symptoms_children ^ 2)))
  duration_of_symptoms_children <- rlnorm(n = n_samples, duration_of_symptoms_children_location,  duration_of_symptoms_children_shape)     #calculated from Violato et al 2019
  rate_of_symptoms_children <- 1 / duration_of_symptoms_children
  probability_late_diagnosis_children <- 1 - exp(-rate_of_symptoms_children)
  
  probability_late_diagnosis <- rep(0, times = n_samples)
  for (i_sample in 1:n_samples) {
  probability_late_diagnosis[i_sample] <- ifelse(population == "adults", probability_late_diagnosis_adults[i_sample], probability_late_diagnosis_children[i_sample])
  }
  
  
  prevalence <- read.csv("data/CPRD prevalence.csv")
  
  #Initial cohort at diagnosis - depends on age at diagnosis
  probability_subfertility <- rbeta(n=n_samples, shape1 = prevalence[prevalence$Age.categories == starting_age, "Subfertility_r"], shape2 = prevalence[prevalence$Age.categories == starting_age, "N"] - prevalence[prevalence$Age.categories == starting_age, "Subfertility_r"])
  probability_osteoporosis <- rbeta(n=n_samples, shape1 = prevalence[prevalence$Age.categories == starting_age, "Osteoporosis_r"], shape2 = prevalence[prevalence$Age.categories == starting_age, "N"] - prevalence[prevalence$Age.categories == starting_age, "Osteoporosis_r"])
  probability_NHL <- rbeta(n=n_samples, shape1 = prevalence[prevalence$Age.categories == starting_age, "NHL_r"], shape2 = prevalence[prevalence$Age.categories == starting_age, "N"] - prevalence[prevalence$Age.categories == starting_age, "NHL_r"])
  probability_nocomplications <- 1 - probability_subfertility - probability_osteoporosis - probability_NHL
  
  probability_IDA_0 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[1], shape2 = (prevalence$N[1] - prevalence$IDA_r[1]))
  probability_IDA_10 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[2], shape2 = (prevalence$N[2] - prevalence$IDA_r[2]))
  probability_IDA_20 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[3], shape2 = (prevalence$N[3] - prevalence$IDA_r[3]))
  probability_IDA_30 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[4], shape2 = (prevalence$N[4] - prevalence$IDA_r[4]))
  probability_IDA_40 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[5], shape2 = (prevalence$N[5] - prevalence$IDA_r[5]))
  probability_IDA_50 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[6], shape2 = (prevalence$N[6] - prevalence$IDA_r[6]))
  probability_IDA_60 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[7], shape2 = (prevalence$N[7] - prevalence$IDA_r[7]))
  probability_IDA_70 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[8], shape2 = (prevalence$N[8] - prevalence$IDA_r[8]))
  probability_IDA_80 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[9], shape2 = (prevalence$N[9] - prevalence$IDA_r[9]))
  probability_IDA_90 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[10], shape2 = (prevalence$N[10] - prevalence$IDA_r[10]))
  probability_IDA <- data.frame(probability_IDA_0 , probability_IDA_10 , probability_IDA_20 , probability_IDA_30 , probability_IDA_40
                                ,probability_IDA_50 , probability_IDA_60 , probability_IDA_70 , probability_IDA_80 , probability_IDA_90)
  
  # Osteoporosis probabilities On GFD
  osteoporosis_probability <- read.csv("data/osteoporosis_rate_nice.csv")
  osteoporosis_probability_GFD_0	<- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_alpha[1]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_beta[1]))
  osteoporosis_probability_GFD_10	<- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_alpha[2]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_beta[2]))
  osteoporosis_probability_GFD_20	<- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_alpha[3]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_beta[3]))
  osteoporosis_probability_GFD_30	<- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_alpha[4]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_beta[4]))
  osteoporosis_probability_GFD_40	<- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_alpha[5]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_beta[5]))
  osteoporosis_probability_GFD_50 <- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_alpha[6]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_beta[6]))
  osteoporosis_probability_GFD_60	<- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_alpha[7]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_beta[7]))
  osteoporosis_probability_GFD_70	<- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_alpha[8]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_beta[8]))
  osteoporosis_probability_GFD_80 <- 	rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_alpha[9]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_beta[9]))
  osteoporosis_probability_GFD_90 <- 	rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_alpha[10]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_GFD_beta[10]))
  osteoporosis_probability_GFD_all <- data.frame(osteoporosis_probability_GFD_0, osteoporosis_probability_GFD_10, osteoporosis_probability_GFD_20, osteoporosis_probability_GFD_30,
                                                 osteoporosis_probability_GFD_40, osteoporosis_probability_GFD_50, osteoporosis_probability_GFD_60,
                                                 osteoporosis_probability_GFD_70, osteoporosis_probability_GFD_80, osteoporosis_probability_GFD_90)
  
  # Subfertility probabilities on GFD
  subfertility_probability <- read.csv("data/subfertility_rate_nice.csv")
  subfertility_probability_GFD_0	<- rbeta(n=n_samples, shape1 = as.numeric(subfertility_probability$subfertility_GFD_alpha[1]), shape2 = as.numeric(subfertility_probability$subfertility_GFD_beta[1]))
  subfertility_probability_GFD_10	<- rbeta(n=n_samples, shape1 = as.numeric(subfertility_probability$subfertility_GFD_alpha[2]), shape2 = as.numeric(subfertility_probability$subfertility_GFD_beta[2]))
  subfertility_probability_GFD_20	<- rbeta(n=n_samples, shape1 = as.numeric(subfertility_probability$subfertility_GFD_alpha[3]), shape2 = as.numeric(subfertility_probability$subfertility_GFD_beta[3]))
  subfertility_probability_GFD_30 <-	rbeta(n=n_samples, shape1 = as.numeric(subfertility_probability$subfertility_GFD_alpha[4]), shape2 = as.numeric(subfertility_probability$subfertility_GFD_beta[4]))
  subfertility_probability_GFD_40 <- rbeta(n=n_samples, shape1 = as.numeric(subfertility_probability$subfertility_GFD_alpha[5]), shape2 = as.numeric(subfertility_probability$subfertility_GFD_beta[5]))
  subfertility_probability_GFD_50 <- 0
  subfertility_probability_GFD_60	<- 0
  subfertility_probability_GFD_70	<- 0
  subfertility_probability_GFD_80 <- 0
  subfertility_probability_GFD_90 <- 0 
  subfertility_probability_GFD_all <- data.frame(subfertility_probability_GFD_0, subfertility_probability_GFD_10, subfertility_probability_GFD_20, subfertility_probability_GFD_30,
                                                 subfertility_probability_GFD_40, subfertility_probability_GFD_50, subfertility_probability_GFD_60,
                                                 subfertility_probability_GFD_70, subfertility_probability_GFD_80, subfertility_probability_GFD_90)
  
  # NHL probabilities on GFD
  NHL_probability <- read.csv("data/NHL_rate_nice.csv")
  NHL_probability_GFD_18orless	<- rbeta(n=n_samples, shape1 = as.numeric(NHL_probability$NHL_GFD_alpha[1]), shape2 = as.numeric(NHL_probability$NHL_GFD_beta[1]))
  NHL_probability_GFD_18plus	<- rbeta(n=n_samples, shape1 = as.numeric(NHL_probability$NHL_GFD_alpha[2]), shape2 = as.numeric(NHL_probability$NHL_GFD_beta[2]))
  NHL_probability_GFD <- rep(0, times = n_samples)
  for (i_sample in 1:n_samples) {
  NHL_probability_GFD[i_sample] <- ifelse(population == "adults", NHL_probability_GFD_18plus[i_sample], NHL_probability_GFD_18orless[i_sample])
  }
  
  # Corresponding probabilities not on GFD
  osteoporosis_probability_noGFD_0	<- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_alpha[1]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_beta[1]))
  osteoporosis_probability_noGFD_10	<- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_alpha[2]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_beta[2]))
  osteoporosis_probability_noGFD_20	<- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_alpha[3]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_beta[3]))
  osteoporosis_probability_noGFD_30	<- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_alpha[4]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_beta[4]))
  osteoporosis_probability_noGFD_40	<- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_alpha[5]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_beta[5]))
  osteoporosis_probability_noGFD_50 <- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_alpha[6]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_beta[6]))
  osteoporosis_probability_noGFD_60	<- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_alpha[7]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_beta[7]))
  osteoporosis_probability_noGFD_70	<- rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_alpha[8]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_beta[8]))
  osteoporosis_probability_noGFD_80 <- 	rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_alpha[9]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_beta[9]))
  osteoporosis_probability_noGFD_90 <- 	rbeta(n=n_samples, shape1 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_alpha[10]), shape2 = as.numeric(osteoporosis_probability$Osteoporosis_noGFD_beta[10]))
  osteoporosis_probability_noGFD_all <- data.frame(osteoporosis_probability_noGFD_0, osteoporosis_probability_noGFD_10, osteoporosis_probability_noGFD_20, osteoporosis_probability_noGFD_30,
                                                   osteoporosis_probability_noGFD_40, osteoporosis_probability_noGFD_50, osteoporosis_probability_noGFD_60,
                                                   osteoporosis_probability_noGFD_70, osteoporosis_probability_noGFD_80, osteoporosis_probability_noGFD_90)
  
  
  
  subfertility_probability_noGFD_0	<- rbeta(n=n_samples, shape1 = as.numeric(subfertility_probability$subfertility_noGFD_alpha[1]), shape2 = as.numeric(subfertility_probability$subfertility_noGFD_beta[1]))
  subfertility_probability_noGFD_10	<- rbeta(n=n_samples, shape1 = as.numeric(subfertility_probability$subfertility_noGFD_alpha[2]), shape2 = as.numeric(subfertility_probability$subfertility_noGFD_beta[2]))
  subfertility_probability_noGFD_20	<- rbeta(n=n_samples, shape1 = as.numeric(subfertility_probability$subfertility_noGFD_alpha[3]), shape2 = as.numeric(subfertility_probability$subfertility_noGFD_beta[3]))
  subfertility_probability_noGFD_30 <-	rbeta(n=n_samples, shape1 = as.numeric(subfertility_probability$subfertility_noGFD_alpha[4]), shape2 = as.numeric(subfertility_probability$subfertility_noGFD_beta[4]))
  subfertility_probability_noGFD_40 <- rbeta(n=n_samples, shape1 = as.numeric(subfertility_probability$subfertility_noGFD_alpha[5]), shape2 = as.numeric(subfertility_probability$subfertility_noGFD_beta[5]))
  subfertility_probability_noGFD_50 <- 0
  subfertility_probability_noGFD_60	<- 0
  subfertility_probability_noGFD_70	<- 0
  subfertility_probability_noGFD_80 <- 0
  subfertility_probability_noGFD_90 <- 0 
  subfertility_probability_noGFD_all <- data.frame(subfertility_probability_noGFD_0, subfertility_probability_noGFD_10, subfertility_probability_noGFD_20, subfertility_probability_noGFD_30,
                                                   subfertility_probability_noGFD_40, subfertility_probability_noGFD_50, subfertility_probability_noGFD_60,
                                                   subfertility_probability_noGFD_70, subfertility_probability_noGFD_80, subfertility_probability_noGFD_90)
  
  
  NHL_probability_noGFD_18orless	<- rbeta(n=n_samples, shape1 = as.numeric(NHL_probability$NHL_noGFD_alpha[1]), shape2 = as.numeric(NHL_probability$NHL_noGFD_beta[1]))
  NHL_probability_noGFD_18plus	<- rbeta(n=n_samples, shape1 = as.numeric(NHL_probability$NHL_noGFD_alpha[2]), shape2 = as.numeric(NHL_probability$NHL_noGFD_beta[2]))
  NHL_probability_noGFD <- rep(0, times = n_samples)
  for (i_sample in 1:n_samples) {
    NHL_probability_noGFD[i_sample] <- ifelse(population == "adults", NHL_probability_noGFD_18plus[i_sample], NHL_probability_noGFD_18orless[i_sample])
  }
  
  
  death_hazard_NHL <- rnorm(n = n_samples, mean = exp(-2.092), sd = 0.006378) #do I exponentiate the sd?
  death_probability_NHL <-	1-exp(-death_hazard_NHL) 

  #############################################################################
  ## Utilities ################################################################
  #############################################################################
  
  utility_GFD_adults <- 0.85
  utility_GFdse_adults <- ((0.86-0.84)/3.92)
  utility_GFDalpha_adults <- (utility_GFD_adults ^ 2 * (1 - utility_GFD_adults)/utility_GFdse_adults ^ 2) - utility_GFD_adults
  utility_GFDbeta_adults <- (utility_GFDalpha_adults / utility_GFD_adults) - utility_GFDalpha_adults
  utility_GFD_adults <- rbeta(n = n_samples, shape1 = utility_GFDalpha_adults, shape2 = utility_GFDbeta_adults)
  
  utility_GFD_children <- 0.88
  utility_GFdse_children <- ((0.92-0.85)/3.92)
  utility_GFDalpha_children <- (utility_GFD_children ^ 2 * (1 - utility_GFD_children)/utility_GFdse_children ^ 2) - utility_GFD_children
  utility_GFDbeta_children <- (utility_GFDalpha_children / utility_GFD_children) - utility_GFDalpha_children
  utility_GFD_children <- rbeta(n = n_samples, shape1 = utility_GFDalpha_children, shape2 = utility_GFDbeta_children)
  
  utility_GFD <- rep(0, times = n_samples)
  for (i_sample in 1:n_samples) {
    utility_GFD[i_sample] <- ifelse(population == "adults", utility_GFD_adults[i_sample], utility_GFD_children[i_sample])
  }
  
  utility_undiagnosedCD_adults <-  0.65 
  utility_undiagnosedCD_se_adults <- (0.67 - 0.63)/3.92
  utility_undiagnosedCD_alpha_adults <- (utility_undiagnosedCD_adults ^ 2 * (1 - utility_undiagnosedCD_adults)/utility_undiagnosedCD_se_adults ^ 2) - utility_undiagnosedCD_adults
  utility_undiagnosedCD_beta_adults <- (utility_undiagnosedCD_alpha_adults/utility_undiagnosedCD_adults) - utility_undiagnosedCD_alpha_adults
  utility_undiagnosedCD_adults <- rbeta(n = n_samples, shape1 = utility_undiagnosedCD_alpha_adults, shape2 = utility_undiagnosedCD_beta_adults)
  
  utility_undiagnosedCD_children <-  0.65 
  utility_undiagnosedCD_se_children <- (0.67 - 0.63)/3.92
  utility_undiagnosedCD_alpha_children <- (utility_undiagnosedCD_children ^ 2 * (1 - utility_undiagnosedCD_children)/utility_undiagnosedCD_se_children ^ 2) - utility_undiagnosedCD_children
  utility_undiagnosedCD_beta_children <- (utility_undiagnosedCD_alpha_children/utility_undiagnosedCD_children) - utility_undiagnosedCD_alpha_children
  utility_undiagnosedCD_children <- rbeta(n = n_samples, shape1 = utility_undiagnosedCD_alpha_children, shape2 = utility_undiagnosedCD_beta_children)
  
  utility_undiagnosedCD <- rep(0, times = n_samples)
  for (i_sample in 1:n_samples) {
    utility_undiagnosedCD[i_sample] <- ifelse(population == "adults", utility_undiagnosedCD_adults[i_sample], utility_undiagnosedCD_children[i_sample])
  }
  
  disutility_subfertility <-  0.158 
  disutility_subfertility_se <- (0.173 - 0.143)/3.92
  disutility_subfertility_alpha <- (disutility_subfertility ^ 2 * (1 - disutility_subfertility)/disutility_subfertility_se ^ 2) - disutility_subfertility
  disutility_subfertility_beta <- (disutility_subfertility_alpha/disutility_subfertility) - disutility_subfertility_alpha
  disutility_subfertility <- rbeta(n = n_samples, shape1 = disutility_subfertility_alpha, shape2 = disutility_subfertility_beta)
  
  probability_hipfracture <- 0.00196
  probability_vertebralfracture <- 0.00071
  probability_wristfracture <- 0.00125
  disutility_hipfracture <- 0.817 - 0.59
  disutility_hipfractureSE <- (((0.817 - 0.65) - (0.817 - 0.54))/3.92)
  disutility_hipfracture_alpha <- (disutility_hipfracture ^ 2 * (1 - disutility_hipfracture)/disutility_hipfractureSE ^ 2) - disutility_hipfracture
  disutility_hipfracture_beta <- (disutility_hipfracture_alpha/disutility_hipfracture) - disutility_hipfracture_alpha
  disutility_hipfracture <- rbeta(n = n_samples, shape1 = disutility_hipfracture_alpha, shape2 = disutility_hipfracture_beta)
  disutility_wristfracture <- 0.817 - 0.55
  disutility_wristfractureSE <- (((0.817 - 0.60) - (0.817 - 0.50))/3.92)
  disutility_wristfracture_alpha <- (disutility_wristfracture ^ 2 * (1 - disutility_wristfracture)/disutility_wristfractureSE ^ 2) - disutility_wristfracture
  disutility_wristfracture_beta <- (disutility_wristfracture_alpha/disutility_wristfracture) - disutility_wristfracture_alpha
  disutility_wristfracture <- rbeta(n = n_samples, shape1 = disutility_wristfracture_alpha, shape2 = disutility_wristfracture_beta)
  disutility_vertebralfracture <- 0.817 - 0.78
  disutility_vertebralfractureSE <- (((0.817 - 0.84) - (0.817 - 0.72))/3.92)
  disutility_vertebralfracture_alpha <- (disutility_vertebralfracture ^ 2 * (1 - disutility_vertebralfracture)/disutility_vertebralfractureSE ^ 2) - disutility_vertebralfracture
  disutility_vertebralfracture_beta <- (disutility_vertebralfracture_alpha/disutility_vertebralfracture) - disutility_vertebralfracture_alpha
  disutility_vertebralfracture <- rbeta(n = n_samples, shape1 = disutility_vertebralfracture_alpha, shape2 = disutility_vertebralfracture_beta)
  disutility_osteoporosis <- (probability_hipfracture * disutility_hipfracture) + (probability_wristfracture * disutility_wristfracture) + (probability_vertebralfracture * disutility_vertebralfracture)
  

  disutility_NHL <- runif(n = n_samples, min = 0.036, max = 0.136) 
  
  disutility_biopsy_adults <- rtri(n = n_samples, min = 0, max = 0.005, mode = 0.003)
  disutility_biopsy_children <- rtri(n = n_samples, min = 0, max = 0.010, mode = 0.006)
  disutility_biopsy <- rep(0, times = n_samples)
  for (i_sample in 1:n_samples) {
    disutility_biopsy[i_sample] <- ifelse(population == "adults", disutility_biopsy_adults[i_sample], disutility_biopsy_children[i_sample])
  }
  
  disutility_biopsy_wait <- (utility_GFD - utility_undiagnosedCD) * 6/52 
 
   disutility_fp <- ifelse(disutility_fp_diagnosis == "Yes", 0.009, 0) #from NICE guideline
  #############################################################################
  ## Costs ####################################################################
  #############################################################################
  
  cost_hipfracture <- 19073
  cost_hipfractureSE <- ((16515 * 1.17) - (16097 * 1.17)) / 3.92
  cost_hipfracture_alpha <- (cost_hipfracture / cost_hipfractureSE) ^ 2
  cost_hipfracture_beta <- (cost_hipfractureSE ^ 2) / cost_hipfracture
  cost_hipfracture <- rgamma(n = n_samples, shape = cost_hipfracture_alpha, scale = cost_hipfracture_beta)
  cost_wristfracture <- 842.40
  cost_wristfractureSE <- cost_wristfracture/10
  cost_wristfracture_alpha <- (cost_wristfracture / cost_wristfractureSE) ^ 2
  cost_wristfracture_beta <- (cost_wristfractureSE ^ 2) / cost_wristfracture
  cost_wristfracture <- rgamma(n = n_samples, shape = cost_wristfracture_alpha, scale = cost_wristfracture_beta)
  cost_vertebralfracture <- 862.20
  cost_vertebralfractureSE <- cost_vertebralfracture/10
  cost_vertebralfracture_alpha <- (cost_vertebralfracture / cost_vertebralfractureSE) ^ 2
  cost_vertebralfracture_beta <- (cost_vertebralfractureSE ^ 2) / cost_vertebralfracture
  cost_vertebralfracture <- rgamma(n = n_samples, shape = cost_vertebralfracture_alpha, scale = cost_vertebralfracture_beta)
  cost_osteoporosis <- (probability_hipfracture * cost_hipfracture) + (probability_wristfracture * cost_wristfracture) + (probability_vertebralfracture * cost_vertebralfracture)
  
  cost_IDA <- if(perspective == "NHS") 0 else 17.89
  cost_gfp <- if(perspective == "NHS") 0 else 100
  
  cost_undiagnosedCD_adults <- 421
  cost_undiagnosedCD_se_adults <- 3.34
  cost_undiagnosedCD_alpha_adults <- (cost_undiagnosedCD_adults/cost_undiagnosedCD_se_adults)^2
  cost_undiagnosedCD_beta_adults <- (cost_undiagnosedCD_se_adults^2)/cost_undiagnosedCD_adults
  cost_undiagnosedCD_adults <- rgamma(n = n_samples, shape = cost_undiagnosedCD_alpha_adults, scale = cost_undiagnosedCD_beta_adults) 
  
  cost_undiagnosedCD_children <- 248
  cost_undiagnosedCD_se_children <- 4.97
  cost_undiagnosedCD_alpha_children <- (cost_undiagnosedCD_children/cost_undiagnosedCD_se_children)^2
  cost_undiagnosedCD_beta_children <- (cost_undiagnosedCD_se_children^2)/cost_undiagnosedCD_children
  cost_undiagnosedCD_children <- rgamma(n = n_samples, shape = cost_undiagnosedCD_alpha_children, scale = cost_undiagnosedCD_beta_children) 
  
  cost_undiagnosedCD <- rep(0, times = n_samples)
  for (i_sample in 1:n_samples) {
    cost_undiagnosedCD[i_sample] <- ifelse(population == "adults", cost_undiagnosedCD_adults[i_sample], cost_undiagnosedCD_children[i_sample])
  }
  
  cost_CDGFD_adults <- 757
  cost_CDGFD_se_adults <- 5.3
  cost_CDGFD_alpha_adults <- (cost_CDGFD_adults/cost_CDGFD_se_adults)^2
  cost_CDGFD_beta_adults <- (cost_CDGFD_se_adults^2)/cost_CDGFD_adults
  cost_CDGFD_adults <- rgamma(n = n_samples, shape = cost_CDGFD_alpha_adults, scale = cost_CDGFD_beta_adults)
  
  cost_CDGFD_children <- 452
  cost_CDGFD_se_children <- 20.6
  cost_CDGFD_alpha_children <- (cost_CDGFD_children/cost_CDGFD_se_children)^2
  cost_CDGFD_beta_children <- (cost_CDGFD_se_children^2)/cost_CDGFD_children
  cost_CDGFD_children <- rgamma(n = n_samples, shape = cost_CDGFD_alpha_children, scale = cost_CDGFD_beta_children)
  
  cost_CDGFD <- rep(0, times = n_samples)
  for (i_sample in 1:n_samples) {
    cost_CDGFD[i_sample] <- ifelse(population == "adults", cost_CDGFD_adults[i_sample], cost_CDGFD_children[i_sample])
  }
  
  
  cost_biopsy_adults <- 530
  cost_biopsy_children <- 823
  cost_biopsy <- ifelse(population == "adults", cost_biopsy_adults, cost_biopsy_children)
  probability_biopsy <- runif(n = n_samples, min = 0.6, max = 0.8)
  
  cost_subfertility <- 8079.75 
  cost_subfertility_se <- (8742.96 - 7432.11)/3.92
  cost_subfertility_alpha <- (cost_subfertility/cost_subfertility_se)^2
  cost_subfertility_beta <- (cost_subfertility_se^2)/cost_subfertility
  cost_subfertility <- rgamma(n = n_samples, shape = cost_subfertility_alpha, scale = cost_subfertility_beta)
  
  cost_NHL <- 18396
  cost_NHL_sd <- sqrt(271) * ((18415 - 18377)/3.92)
  cost_NHL_location <- log(cost_NHL^2 / sqrt(cost_NHL_sd^2 + cost_NHL^2))
  cost_NHL_shape <- sqrt(log(1 + (cost_NHL_sd^2 / cost_NHL^2)))
  cost_NHL <- rlnorm(n = n_samples, cost_NHL_location,  cost_NHL_shape)
  
  cost_diagnosis <- 379.27
  cost_diagnosis_sd <- cost_diagnosis/10
  cost_diagnosis_location <- log(cost_diagnosis^2 / sqrt(cost_diagnosis_sd^2 + cost_diagnosis^2))
  cost_diagnosis_shape <- sqrt(log(1 + (cost_diagnosis_sd^2 / cost_diagnosis^2)))
  cost_diagnosis <- rlnorm(n = n_samples, cost_diagnosis_location,  cost_diagnosis_shape)
  

  #Need to include cost of GFD - Penny contacting Coeliac UK

  test_cost_IgAEMA <- 14.92 #based on correspondence from labs 
  test_cost_IgAEMA_se <- test_cost_IgAEMA/8  #this gives a min of 9.5 and max of 21 which is in line with lab estimates
  test_cost_IgAEMA_alpha <- (test_cost_IgAEMA/test_cost_IgAEMA_se)^2
  test_cost_IgAEMA_beta <- (test_cost_IgAEMA_se^2)/test_cost_IgAEMA
  test_cost_IgAEMA <- rgamma(n = n_samples, shape = test_cost_IgAEMA_alpha, scale = test_cost_IgAEMA_beta)
  #rgamma(n = n_samples, shape = 122.57, scale = 0.08) #from NICE model
  test_cost_IgATTG <- 10.77 #based on correspondence from labs
  test_cost_IgATTG_se <- test_cost_IgATTG/5
  test_cost_IgATTG_alpha <- (test_cost_IgATTG/test_cost_IgATTG_se)^2
  test_cost_IgATTG_beta <- (test_cost_IgATTG_se^2)/test_cost_IgATTG
  test_cost_IgATTG <- rgamma(n = n_samples, shape = test_cost_IgATTG_alpha, scale = test_cost_IgATTG_beta) #this gives a min of 5.04 and max of 18.21 which is in line with lab estimates
  #rgamma(n = n_samples, shape = 26.16, scale = 0.42) #from NICE model
  test_cost_HLA <- 122.34 #correspondence from labs
  test_cost_HLA_se <- test_cost_HLA/5  #this gives a min of 61.15 and max of 214.29 which is in line with lab estimates
  test_cost_HLA_alpha <- (test_cost_HLA/test_cost_HLA_se)^2
  test_cost_HLA_beta <- (test_cost_HLA_se^2)/test_cost_HLA
  test_cost_HLA <- rgamma(n = n_samples, shape = test_cost_HLA_alpha, scale = test_cost_HLA_beta)
  
  #capital_cost_double_test <- 0.44 #based on NICE guideline capital cost for IgATTG + IgAEMA inflated to 2021 prices
  #############################################################################
  ## Accuracy of tests ########################################################
  #############################################################################
  sens_biopsy <- 1 #assumption to be varied in sensitivity analysis
  spec_biopsy <- 1 #assumption to be varied in sensitivity analysis

  
  pre_test_probability_overall <- rbeta(n = n_samples, shape1 = 10872, shape2 = 4519128) #based on West 2014
  tp_riskfactor <- array(dim=c(n_samples, n_combinations), dimnames = list(NULL, paste(combinations_names, "tp_riskfactor")))
    fn_riskfactor <- array(dim=c(n_samples, n_combinations), dimnames = list(NULL, paste(combinations_names, "fn_riskfactor")))
      fp_riskfactor <- array(dim=c(n_samples, n_combinations), dimnames = list(NULL, paste(combinations_names, "fp_riskfactor")))
        tn_riskfactor <- array(dim=c(n_samples, n_combinations), dimnames = list(NULL, paste(combinations_names, "tn_riskfactor")))
  
  for (i in 1:n_combinations) {
    tp_riskfactor[,i] <- pre_test_probability_overall * combinations$sens_riskfactor[i]
    fn_riskfactor[,i] <- pre_test_probability_overall - tp_riskfactor[,i]  
    tn_riskfactor[,i] <- (1 - pre_test_probability_overall) * combinations$spec_riskfactor[i]
    fp_riskfactor[,i] <- (1 - pre_test_probability_overall) - tn_riskfactor[,i]
  }
  
 write.csv(data.frame(tp_riskfactor[1,], fn_riskfactor[1,], tn_riskfactor[1,], fp_riskfactor[1,]), "risk_factor.csv")
  pre_test_probability <- tp_riskfactor/(tp_riskfactor+fp_riskfactor)
  colnames(pre_test_probability) <- paste(combinations_names, "pre_test_probability")
  write.csv(colMeans(pre_test_probability), "pretestprob.csv")
  
  pre_test_odds <- array(0, dim=c(n_samples, n_combinations), dimnames = list(NULL, combinations_names))
  
  for (i in 1:n_combinations){
    pre_test_odds[,i] <- pre_test_probability[,i]/(1 - pre_test_probability[,i])
  }

  #################################################################################################################
  #Iga EMA adults
  #E(logitSE) coef = 1.993122, SE = 0.4508497
  #E(logitSP) coef = 5.54022, SE = 1.556019
  #Covariance = -0.2689103
  
  #random normal values with mean [1.993122, 5.54022] and SEs [0.4508497, 1.556019], and covariance -0.2689103
  sigma_IgAEMA_adults <- matrix(c((0.4508497^2),-0.2689103,-0.2689103,(1.556019^2)), 2, 2)
  mu_IgAEMA_adults <- c(1.993122, 5.54022)
  x_IgAEMA_adults <- rmvnorm(n_samples, mu_IgAEMA_adults, sigma_IgAEMA_adults)
  head(x_IgAEMA_adults)
  SensSpec_IgAEMA_adults <- exp(x_IgAEMA_adults)/(1+exp(x_IgAEMA_adults))
  sens_IgAEMA_adults <- SensSpec_IgAEMA_adults[,1]
  spec_IgAEMA_adults <- SensSpec_IgAEMA_adults[,2]
  LR_IgAEMA_adults <- SensSpec_IgAEMA_adults[,1]/ (1 - SensSpec_IgAEMA_adults[,2])
  
  post_test_odds_IgAEMA_adults <- array(dim=c(n_samples, n_combinations),dimnames=list(NULL, combinations_names))
  post_test_probability_IgAEMA_adults <- array(dim=c(n_samples, n_combinations),dimnames=list(NULL, combinations_names))
  
  for (i in 1:n_combinations){
  post_test_odds_IgAEMA_adults[,i] <- pre_test_odds[,i] * LR_IgAEMA_adults
  post_test_probability_IgAEMA_adults[,i] <- post_test_odds_IgAEMA_adults[,i]/(1 + post_test_odds_IgAEMA_adults[,i])
  }
  
 
  #Iga EMA children
  #E(logitSE) coef = 2.839716, SE = 0.3886658
  #E(logitSP) coef = 2.716697, SE = 0.4927015
  #Covariance = 0.1592634
  
  #random normal values with mean [2.839716, 2.716697] and variances [0.3886658, 0.4927015], and covariance 0.1592634
  sigma <- matrix(c((0.3886658^2),0.1592634,0.1592634,(0.4927015^2)), 2, 2)
  mu <- c(2.839716, 2.716697)
  x <- rmvnorm(n_samples, mu, sigma)
  head(x)
  SensSpec_IgAEMA_children <- exp(x)/(1+exp(x))
  sens_IgAEMA_children <- SensSpec_IgAEMA_children[,1]
  spec_IgAEMA_children <- SensSpec_IgAEMA_children[,2]
  LR_IgAEMA_children <- SensSpec_IgAEMA_children[,1]/ (1 - SensSpec_IgAEMA_children[,2])
  
  post_test_odds_IgAEMA_children <- array(dim=c(n_samples, n_combinations),dimnames=list(NULL, combinations_names))
  post_test_probability_IgAEMA_children <- array(dim=c(n_samples, n_combinations),dimnames=list(NULL, combinations_names))
  
  for (i in 1:n_combinations){
    post_test_odds_IgAEMA_children[,i] <- pre_test_odds[,i] * LR_IgAEMA_children
    post_test_probability_IgAEMA_children[,i] <- post_test_odds_IgAEMA_children[,i]/(1 + post_test_odds_IgAEMA_children[,i])
  }

  
  post_test_probability_IgAEMA <- array(dim=c(n_samples, n_combinations),dimnames=list(NULL, paste(combinations_names, "IgAEMA")))
  sens_IgAEMA <- rep(0, times = n_samples)
  spec_IgAEMA <- rep(0, times = n_samples)
  for(i_sample in 1:n_samples){
  for (i in 1:n_combinations){
  post_test_probability_IgAEMA[i_sample,i] <- ifelse(population == "adults", post_test_probability_IgAEMA_adults[i_sample,i], post_test_probability_IgAEMA_children[i_sample,i])
  }}

  for(i_sample in 1:n_samples){
  sens_IgAEMA[i_sample] <- ifelse(population == "adults", sens_IgAEMA_adults[i_sample], sens_IgAEMA_children[i_sample])
  spec_IgAEMA[i_sample] <- ifelse(population == "adults", spec_IgAEMA_adults[i_sample], spec_IgAEMA_children[i_sample])
  }
  ##################################################################################################
  
  #IgATTGplusEMA
  #E(logitSE) coef = 1.939019, SE = 0.3723165
  #E(logitSP) coef = 4.252873, SE = 0.5569404
  #Covariance = -0.358947
  
  #random normal values with mean [1.939019, 4.252873] and variances [0.3723165, 0.5569404], and covariance -0.358947
  sigma <- matrix(c((0.3723165^2),-0.358947,-0.358947, (0.5569404^2)), 2, 2)
  mu <- c(1.939019, 4.252873)
  x <- rmvnorm(n_samples, mu, sigma)
  head(x)
  SensSpec_IgATTGplusEMA <- exp(x)/(1+exp(x))
  sens_IgATTGplusEMA <- SensSpec_IgATTGplusEMA[,1]
  spec_IgATTGplusEMA <- SensSpec_IgATTGplusEMA[,2]
  LR_IgATTGplusEMA <- SensSpec_IgATTGplusEMA[,1]/ (1 - SensSpec_IgATTGplusEMA[,2])
  
  post_test_odds_IgATTGplusEMA <- array(dim=c(n_samples, n_combinations),dimnames=list(NULL, combinations_names))
  post_test_probability_IgATTGplusEMA <- array(dim=c(n_samples, n_combinations),dimnames=list(NULL, paste(combinations_names, "IgATTGplusEMA")))
  
  for (i in 1:n_combinations){
    post_test_odds_IgATTGplusEMA[,i] <- pre_test_odds[,i] * LR_IgATTGplusEMA[1:n_samples]
    post_test_probability_IgATTGplusEMA[,i] <- post_test_odds_IgATTGplusEMA[,i]/(1 + post_test_odds_IgATTGplusEMA[,i])
  }

  ######################################################################################################################
  
  #IgATTG
  #E(logitSE) coef = 2.272053, SE = 0.173953
  #E(logitSP) coef = 1.940514, se = 0.1290671
  #Covariance = 0.0019497
  sigma <- matrix(c((0.173953^2),0.0019497,0.0019497,(0.1290671^2)), 2, 2)
  mu <- c(2.272053, 1.940514)
  x <- rmvnorm(n_samples, mu, sigma)
  head(x)
  SensSpec_IgATTG <- exp(x)/(1+exp(x))
  sens_IgATTG <- SensSpec_IgATTG[,1]
  spec_IgATTG <- SensSpec_IgATTG[,2]
  LR_IgATTG <- SensSpec_IgATTG[,1]/ (1 - SensSpec_IgATTG[,2])
  
  post_test_odds_IgATTG <- array(dim=c(n_samples, n_combinations),dimnames=list(NULL, combinations_names))
  post_test_probability_IgATTG <- array(dim=c(n_samples, n_combinations),dimnames=list(NULL, paste(combinations_names, "IgATTG")))
  
  for (i in 1:n_combinations){
    post_test_odds_IgATTG[,i] <- pre_test_odds[,i] * LR_IgATTG
    post_test_probability_IgATTG[,i] <- post_test_odds_IgATTG[,i]/(1 + post_test_odds_IgATTG[,i])
  }

  ######################################################################################################################
  
  #HLA
  #E(logitSE) coef = 3.145430687, SE = 0.704324701
  #E(logitSP) coef = 0.568258353, 0.169650745

  #Covariance = -0.526523769

  
  #random normal values with mean [3.145430687, 0.568258353] and variances [0.704324701, 0.169650745], and covariance -0.526523769
  sigma_HLA <- matrix(c((0.704324701^2),(-0.526523769*0.704324701*0.169650745),(-0.526523769*0.704324701*0.169650745),(0.169650745^2)), 2, 2)
  mu_HLA <- c(3.145430687, 0.568258353)
  x_HLA <- rmvnorm(n_samples, mu_HLA, sigma_HLA)
  head(x_HLA)
  SensSpec_HLA <- exp(x_HLA)/(1+exp(x_HLA))
  sens_HLA <- SensSpec_HLA[,1]
  spec_HLA <- SensSpec_HLA[,2]
  LR_HLA <- SensSpec_HLA[,1]/ (1 - SensSpec_HLA[,2])

  
  
  return(data.frame(probability_late_diagnosis, probability_subfertility, probability_osteoporosis, probability_NHL, probability_nocomplications,
                    osteoporosis_probability_GFD_all, subfertility_probability_GFD_all, NHL_probability_GFD, osteoporosis_probability_noGFD_all, subfertility_probability_noGFD_all,
                    NHL_probability_noGFD, death_probability_NHL, 
                    utility_GFD, utility_undiagnosedCD, disutility_subfertility, disutility_osteoporosis, disutility_NHL, disutility_biopsy, disutility_biopsy_wait, disutility_fp,
                    cost_CDGFD, cost_osteoporosis, cost_undiagnosedCD, cost_IDA, cost_biopsy, probability_biopsy,
                    cost_subfertility, cost_NHL, probability_IDA, cost_diagnosis, test_cost_IgAEMA, test_cost_IgATTG, test_cost_HLA, 
                    sens_IgATTGplusEMA, spec_IgATTGplusEMA, sens_IgAEMA, spec_IgAEMA, sens_IgATTG, spec_IgATTG, cost_gfp, sens_biopsy, spec_biopsy, 
                    post_test_probability_IgAEMA, post_test_probability_IgATTGplusEMA, post_test_probability_IgATTG, pre_test_probability, pre_test_probability_overall,
                   tp_riskfactor, fn_riskfactor, fp_riskfactor, LR_HLA, sens_HLA, spec_HLA))
}

generate_model_parameters(starting_age)