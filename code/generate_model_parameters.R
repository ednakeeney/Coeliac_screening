# Utility functions
expit <- function(logO) {
  return(exp(logO)/(1 + exp(logO)))
}
logit <- function(p) {
  return(log(p/(1-p)))
}

generate_model_parameters <- function(starting_age, population = NULL,
                                      combinations, n_samples,
                                      hold_constant = c()) {
  
  # Derived variables
  n_combinations <- dim(combinations)[1]
  
 starting_age <- ifelse(population == "men" | population == "women", 18, 10) # child starting age based on mean age in under 18s in CPRD cost data
 
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
  probability_late_diagnosis[i_sample] <- ifelse(population == "men" | population == "women", probability_late_diagnosis_adults[i_sample], probability_late_diagnosis_children[i_sample])
  }
  
  # Use the appropriate prevalence
  prevalence_mixed <- as.data.frame(readxl::read_excel(path = "data/CPRD prevalence.xlsx", sheet = "mixed"))
  prevalence_men <- as.data.frame(readxl::read_excel(path = "data/CPRD prevalence.xlsx", sheet = "men"))
  prevalence_women <- as.data.frame(readxl::read_excel(path = "data/CPRD prevalence.xlsx", sheet = "women"))
  if(population == "children") prevalence = prevalence_mixed
  if(population == "men") prevalence = prevalence_men
  if(population == "women") prevalence = prevalence_women
  # For backward compatibility rename prevalence column names
  colnames(prevalence) <- gsub(" ", ".", colnames(prevalence))
  
  
  #Initial cohort at diagnosis - depends on age at diagnosis
  # e.g. 10 year old starts with 10-20 prevalence and 18 year old starts with 10-20 prevalence. 
  starting_age_category <- which(prevalence$Age.categories > starting_age)[1] - 1
  # prevalence$Age.categories == starting_age
  probability_osteoporosis <- rbeta(n=n_samples, shape1 = prevalence[starting_age_category, "Osteoporosis_r"], shape2 = prevalence[starting_age_category, "N"] - prevalence[starting_age_category, "Osteoporosis_r"])
  probability_NHL <- rbeta(n=n_samples, shape1 = prevalence[starting_age_category, "NHL_r"], shape2 = prevalence[starting_age_category, "N"] - prevalence[starting_age_category, "NHL_r"])
  probability_nocomplications <- 1 - probability_osteoporosis - probability_NHL
  
  # IDA prevalence changes with age of cohort so is age stratified
  # Prevalence of osteoporosis and NHL are combined with incidence rates to model prevalence changing with age
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
  if(population == "children") {
    osteoporosis_probability <- as.data.frame(readxl::read_excel(path = "data/osteoporosis_rate.xlsx", sheet = "mixed"))
  } else {
    osteoporosis_probability <- as.data.frame(readxl::read_excel(path = "data/osteoporosis_rate.xlsx", sheet = population))
  }
  # Add correction for zero rates
  osteoporosis_probability$Osteoporosis_rate_low[osteoporosis_probability$Osteoporosis_rate_low == 0] <- osteoporosis_probability$Osteoporosis_rate_low[osteoporosis_probability$Osteoporosis_rate_low == 0] + 0.000001
  
  # Log odds ratio for diagnosed CD
  log_or_osteoporosis_GFD <- rnorm(n_samples, mean = log(1.40), sd = ((log(1.50) - log(1.30))/(2*1.96)))
  log_or_osteoporosis_noGFD <- rnorm(n_samples, mean = log(2.59), sd = ((log(5.09) - log(1.32))/(2*1.96)))

  # Log rates in general population
  osteoporosis_lograte <- matrix(NA, nrow = n_samples, ncol = 10)
  colnames(osteoporosis_lograte) <- paste0("osteoporosis_lograte_", seq(0, 90, 10))
  for(i_age_category in 1:10) {
    osteoporosis_lograte[, i_age_category] <- rnorm(n_samples, mean = log(osteoporosis_probability$Osteoporosis_rate[i_age_category]),
                                                    sd = (log(osteoporosis_probability$Osteoporosis_rate_high[i_age_category]) - log(osteoporosis_probability$Osteoporosis_rate_low[i_age_category]))/(2*1.96))
  }
  
  
  if(population == "children") {
    NHL_probability <- as.data.frame(readxl::read_excel(path = "data/NHL_rate.xlsx", sheet = "mixed"))
  } else {
      NHL_probability <- as.data.frame(readxl::read_excel(path = "data/NHL_rate.xlsx", sheet = population))
  }
  # Ensure no rates are zero
  NHL_probability$NHL_rate_low[NHL_probability$NHL_rate_low == 0] <- 0.000001
  
  # Log incidence ratios for NHL
  log_rr_NHL_GFD <- rnorm(n_samples, mean = log(3.28), sd = (log(6.28) - log(1.49))/(2*1.96))
  log_rr_NHL_noGFD <- rnorm(n_samples, mean = log(4.7), sd = (log(7.3) - log(2.9))/(2*1.96))
  # NHL probabilities on GFD
  NHL_lograte <- as.data.frame(matrix(NA, nrow = n_samples, ncol = dim(NHL_probability)[1]))
  colnames(NHL_lograte) <- c("NHL_lograte_18orless", "NHL_lograte_18plus")
  
  for(i_age_category in 1:2) {
    NHL_lograte[, i_age_category] <- rnorm(n_samples, mean = log(NHL_probability$NHL_rate[i_age_category]),
                                                    sd = (log(NHL_probability$NHL_rate_high[i_age_category]) - log(NHL_probability$NHL_rate_low[i_age_category]))/(2*1.96))
  }
  
  
  
  death_log_hazard_NHL <- rnorm(n = n_samples, mean = exp(-2.092), sd = 0.006378)
  death_probability_NHL <-	1-exp(-exp(death_log_hazard_NHL))
  
  # HR 3.5 (95% CI: 3.28–3.74)
  death_log_hr_osteoporosis_male <- rnorm(n = n_samples, mean = log(3.5),
                                          sd = (log(3.74) - log(3.28))/(2*1.96))
  # 2.4 (95% CI: 2.31–2.50)
  death_log_hr_osteoporosis_female <- rnorm(n = n_samples, mean = log(2.4),
                                          sd = (log(2.50) - log(2.31))/(2*1.96))
  
  
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
    utility_GFD[i_sample] <- ifelse(population == "men" | population == "women", utility_GFD_adults[i_sample], utility_GFD_children[i_sample])
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
    utility_undiagnosedCD[i_sample] <- ifelse(population == "men" | population == "women", utility_undiagnosedCD_adults[i_sample], utility_undiagnosedCD_children[i_sample])
  }
  
  
  
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
  
  disutility_biopsy_adults <- rtriangle(n = n_samples, a = 0, b = 0.005, c = 0.003)
  disutility_biopsy_children <- rtriangle(n = n_samples, a = 0, b = 0.010, c = 0.006)
  disutility_biopsy <- rep(0, times = n_samples)
  for (i_sample in 1:n_samples) {
    disutility_biopsy[i_sample] <- ifelse(population == "men" | population == "women", disutility_biopsy_adults[i_sample], disutility_biopsy_children[i_sample])
  }
  
  disutility_biopsy_wait <- (utility_GFD - utility_undiagnosedCD) * 6/52 

  disutility_fp_temp <-rnorm(n_samples, mean = -8.3, sd = 3.83) * rnorm(n_samples, mean = 0.0011, sd = 0.0002) 
   disutility_fp <- (disutility_fp_diagnosis == "Yes") * disutility_fp_temp + 0 # NICE guidelines 2015
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
  cost_gfp <- if(perspective == "NHS") 0 else 0 # Not used
  
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
    cost_undiagnosedCD[i_sample] <- ifelse(population == "men" | population == "women", cost_undiagnosedCD_adults[i_sample], cost_undiagnosedCD_children[i_sample])
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
    cost_CDGFD[i_sample] <- ifelse(population == "men" | population == "women", cost_CDGFD_adults[i_sample], cost_CDGFD_children[i_sample])
  }
  
  
  cost_biopsy_adults <- 530
  cost_biopsy_children <- 823
  cost_biopsy <- ifelse(population == "men" | population == "women", cost_biopsy_adults, cost_biopsy_children)
  probability_biopsy <- runif(n = n_samples, min = 0.6, max = 0.8)
  
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
  sens_biopsy <- 1 # Assumed perfectly accurate
  spec_biopsy <- 1 # Assumed perfectly accurate

  
  pre_test_probability_overall <- rbeta(n = n_samples, shape1 = 10872, shape2 = 4519128) #based on West 2014
  


  #################################################################################################################
 
   #Iga EMA adults
  #E(logitSE) coef = 1.993122, SE = 0.4508497
  #E(logitSP) coef = 5.54022, SE = 1.556019
  #Covariance = -0.2689103

   #Sens 87% (77.7–92.8), Spec	 98% (97.4–98.6) Hopper 2008

  
  sens_igaema_adults <- logit(0.87)
  sens_igaema_adults_lci <- logit(0.777)
  sens_igaema_adults_uci <- logit(0.928)
 sd = (sens_igaema_adults_uci-sens_igaema_adults_lci)/3.92
 sens_igaema_adults <- rnorm(n=n_samples, sens_igaema_adults, sd=sd)
 sens_IgAEMA_adults <- expit(sens_igaema_adults)
 
 spec_igaema_adults <- logit(0.98)
 spec_igaema_adults_lci <- logit(0.974)
 spec_igaema_adults_uci <- logit(0.986)
 sd = (spec_igaema_adults_uci-spec_igaema_adults_lci)/3.92
 spec_igaema_adults <- rnorm(n=n_samples, spec_igaema_adults, sd=sd)
 spec_IgAEMA_adults <- expit(spec_igaema_adults)
 
  

 
  #Iga EMA children
  #E(logitSE) coef = 2.839716, SE = 0.3886658
  #E(logitSP) coef = 2.716697, SE = 0.4927015
  #Covariance = 0.1592634
  
 # 95.8% (93.7 - 97.2)	94% (91 - 96.1) Wolf 2017
  
  sens_igaema_children <- logit(0.958)
  sens_igaema_children_lci <- logit(0.937)
  sens_igaema_children_uci <- logit(0.972)
  sd = (sens_igaema_children_uci-sens_igaema_children_lci)/3.92
  sens_igaema_children <- rnorm(n=n_samples, sens_igaema_children, sd=sd)
  sens_IgAEMA_children <- expit(sens_igaema_children)
  
  spec_igaema_children <- logit(0.94)
  spec_igaema_children_lci <- logit(0.91)
  spec_igaema_children_uci <- logit(0.961)
  sd = (spec_igaema_children_uci-spec_igaema_children_lci)/3.92
  spec_igaema_children <- rnorm(n=n_samples, spec_igaema_children, sd=sd)
  spec_IgAEMA_children <- expit(spec_igaema_children)
  
  #random normal values with mean [2.839716, 2.716697] and variances [0.3886658, 0.4927015], and covariance 0.1592634
 # sigma <- matrix(c((0.3886658^2),0.1592634,0.1592634,(0.4927015^2)), 2, 2)
  #mu <- c(2.839716, 2.716697)
  #x <- rmvnorm(n_samples, mu, sigma)
  #head(x)
  #SensSpec_IgAEMA_children <- exp(x)/(1+exp(x))
  #sens_IgAEMA_children <- SensSpec_IgAEMA_children[,1]
  #spec_IgAEMA_children <- SensSpec_IgAEMA_children[,2]
  

  
  ##################################################################################################
  
  #IgATTGplusEMA_adults
  #E(logitSE) coef = 1.939019, SE = 0.3723165
  #E(logitSP) coef = 4.252873, SE = 0.5569404
  #Covariance = -0.358947
  
  #85.7% (76.2–91.8) 	98.6% (98 –99) Hopper 2008
  
  sens_igattgplusema_adults <- logit(0.857)
  sens_igattgplusema_adults_lci <- logit(0.762)
  sens_igattgplusema_adults_uci <- logit(0.918)
  sd = (sens_igattgplusema_adults_uci-sens_igattgplusema_adults_lci)/3.92
  sens_igattgplusema_adults <- rnorm(n=n_samples, sens_igattgplusema_adults, sd=sd)
  sens_IgATTGplusEMA_adults <- expit(sens_igattgplusema_adults)
  
  spec_igattgplusema_adults <- logit(0.986)
  spec_igattgplusema_adults_lci <- logit(0.98)
  spec_igattgplusema_adults_uci <- logit(0.99)
  sd = (spec_igattgplusema_adults_uci-spec_igattgplusema_adults_lci)/3.92
  spec_igattgplusema_adults <- rnorm(n=n_samples, spec_igattgplusema_adults, sd=sd)
  spec_IgATTGplusEMA_adults <- expit(spec_igattgplusema_adults)
  
  #random normal values with mean [1.939019, 4.252873] and variances [0.3723165, 0.5569404], and covariance -0.358947
  #sigma <- matrix(c((0.3723165^2),-0.358947,-0.358947, (0.5569404^2)), 2, 2)
  #mu <- c(1.939019, 4.252873)
  #x <- rmvnorm(n_samples, mu, sigma)
  #head(x)
  #SensSpec_IgATTGplusEMA <- exp(x)/(1+exp(x))
  #sens_IgATTGplusEMA <- SensSpec_IgATTGplusEMA[,1]
  #spec_IgATTGplusEMA <- SensSpec_IgATTGplusEMA[,2]
  
  #IgATTGplusEMA_children
 
  
  #95.1% (92.9 - 96.8)	94.5% (91.5 - 96.7) Wolf 2017
  
  sens_igattgplusema_children<- logit(0.951)
  sens_igattgplusema_children_lci <- logit(0.929)
  sens_igattgplusema_children_uci <- logit(0.968)
  sd = (sens_igattgplusema_children_uci-sens_igattgplusema_children_lci)/3.92
  sens_igattgplusema_children<- rnorm(n=n_samples, sens_igattgplusema_children, sd=sd)
  sens_IgATTGplusEMA_children <- expit(sens_igattgplusema_children)
  
  spec_igattgplusema_children <- logit(0.945)
  spec_igattgplusema_children_lci <- logit(0.915)
  spec_igattgplusema_children_uci <- logit(0.967)
  sd = (spec_igattgplusema_children_uci-spec_igattgplusema_children_lci)/3.92
  spec_igattgplusema_children <- rnorm(n=n_samples, spec_igattgplusema_children, sd=sd)
  spec_IgATTGplusEMA_children <- expit(spec_igattgplusema_children)
  
 
  

  ######################################################################################################################
  
  #Iga TTG adults

  
  #90.9% (82.4–94.5) 	90.9% (89.5–92.1) Hopper 2008

  sens_igattg_adults <- logit(0.909)
  sens_igattg_adults_lci <- logit(0.824)
  sens_igattg_adults_uci <- logit(0.945)
  sd = (sens_igattg_adults_uci-sens_igattg_adults_lci)/3.92
  sens_igattg_adults <- rnorm(n=n_samples, sens_igattg_adults, sd=sd)
  sens_IgATTG_adults <- expit(sens_igattg_adults)
  
  spec_igattg_adults <- logit(0.909)
  spec_igattg_adults_lci <- logit(0.895)
  spec_igattg_adults_uci <- logit(0.921)
  sd = (spec_igattg_adults_uci-spec_igattg_adults_lci)/3.92
  spec_igattg_adults <- rnorm(n=n_samples, spec_igattg_adults, sd=sd)
  spec_IgATTG_adults <- expit(spec_igattg_adults)
  
  
  
  #Iga TTG children
  
  
  # 97.1% (95.3 -98.3)	89.3% (85.5-92.1)Wolf 2017
  
  sens_igaTTG_children <- logit(0.971)
  sens_igaTTG_children_lci <- logit(0.953)
  sens_igaTTG_children_uci <- logit(0.983)
  sd = (sens_igaTTG_children_uci-sens_igaTTG_children_lci)/3.92
  sens_igaTTG_children <- rnorm(n=n_samples, sens_igaTTG_children, sd=sd)
  sens_IgATTG_children <- expit(sens_igaTTG_children)
  
  spec_igaTTG_children <- logit(0.893)
  spec_igaTTG_children_lci <- logit(0.855)
  spec_igaTTG_children_uci <- logit(0.921)
  sd = (spec_igaTTG_children_uci-spec_igaTTG_children_lci)/3.92
  spec_igaTTG_children <- rnorm(n=n_samples, spec_igaTTG_children, sd=sd)
  spec_IgATTG_children <- expit(spec_igaTTG_children)
  

  ######################################################################################################################
  
  #HLA
  # Sent by Martha Elwenspoek on 3rd November 2021
  # sens		  spec    	log_sens	log_spec	log_sens_se	log_spec_se	corrmean
  # 0.992431	0.556197	4.876031	0.225742	1.663271	  0.110675  	-0.60621
  mu_HLA <- c(4.876031, 0.22572)
  sigma_HLA <- matrix(c(1.663271^2, rep(1.663271 * 0.110675 * (-0.60621), 2), 0.110675^2), nrow = 2)
  x_HLA <- rmvnorm(n_samples, mu_HLA, sigma_HLA)
  SensSpec_HLA <- expit(x_HLA)
  sens_HLA <- SensSpec_HLA[,1]
  spec_HLA <- SensSpec_HLA[,2]
  
  
  
  input_parameters_temp <- data.frame(probability_late_diagnosis, probability_osteoporosis, probability_NHL, probability_nocomplications,
                    osteoporosis_lograte, log_or_osteoporosis_GFD, log_or_osteoporosis_noGFD,
                    NHL_lograte, log_rr_NHL_GFD, log_rr_NHL_noGFD,
                    death_probability_NHL, death_log_hr_osteoporosis_male, death_log_hr_osteoporosis_female,
                    utility_GFD, utility_undiagnosedCD, disutility_osteoporosis, disutility_NHL, disutility_biopsy, disutility_biopsy_wait, disutility_fp,
                    cost_CDGFD, cost_osteoporosis, cost_undiagnosedCD, cost_IDA, cost_biopsy, probability_biopsy,
                    cost_NHL, probability_IDA, cost_diagnosis, test_cost_IgAEMA, test_cost_IgATTG, test_cost_HLA, 
                    sens_IgATTGplusEMA_adults, spec_IgATTGplusEMA_adults, sens_IgAEMA_adults, 
                    spec_IgAEMA_adults, sens_IgATTG_adults, spec_IgATTG_adults, 
                    sens_IgATTGplusEMA_children, spec_IgATTGplusEMA_children, sens_IgAEMA_children, 
                    spec_IgAEMA_children, sens_IgATTG_children, spec_IgATTG_children, 
                    cost_gfp, sens_biopsy, spec_biopsy, 
                     sens_HLA, spec_HLA,
                    pre_test_probability_overall)
  
  if(length(hold_constant) != 0) {
    input_parameters_temp[, hold_constant] <- input_parameters_temp[1, hold_constant]
  }
  
  return(input_parameters_temp)
}

