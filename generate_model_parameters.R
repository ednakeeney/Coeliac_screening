generate_model_parameters <- function(starting_age) {
  
 
  
   #adherence	<- rbeta(n = n_samples, shape1 = 75, shape2 = 100)       #will come from targeted review
  
  duration_of_symptoms <- 10.93
  duration_of_symptoms_sd <- 13.10
  duration_of_symptoms_location <- log(duration_of_symptoms ^ 2 / sqrt(duration_of_symptoms_sd ^ 2 + duration_of_symptoms ^ 2))
  duration_of_symptoms_shape <- sqrt(log(1 + (duration_of_symptoms_sd ^ 2 / duration_of_symptoms ^ 2)))
  duration_of_symptoms <- rlnorm(n = n_samples, duration_of_symptoms_location,  duration_of_symptoms_shape)     #calculated from Violato et al 2019
  rate_of_symptoms <- 1 / duration_of_symptoms
  probability_late_diagnosis <- 1 - exp(-rate_of_symptoms)
  
  
  prevalence <- read.csv("CPRD prevalence.csv")
  
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
  osteoporosis_probability <- read.csv("osteoporosis_rate_nice.csv")
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
  subfertility_probability <- read.csv("subfertility_rate_nice.csv")
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
  NHL_probability <- read.csv("NHL_rate_nice.csv")
  NHL_probability_GFD_18orless	<- rbeta(n=n_samples, shape1 = as.numeric(NHL_probability$NHL_GFD_alpha[1]), shape2 = as.numeric(NHL_probability$NHL_GFD_beta[1]))
  NHL_probability_GFD_18plus	<- rbeta(n=n_samples, shape1 = as.numeric(NHL_probability$NHL_GFD_alpha[2]), shape2 = as.numeric(NHL_probability$NHL_GFD_beta[2]))
  NHL_probability_GFD <- NHL_probability_GFD_18plus
  
  
  
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
  NHL_probability_noGFD <- NHL_probability_noGFD_18plus
  
 
  #lifetables <- read.csv("lifetables.csv")
  #percentage_male <- 0.5
  #lifetables$Overall <- (percentage_male * lifetables$Males) + ((1-percentage_male) * lifetables$Females)
  #death_probability_nocomplications	<- data.frame(lifetables$Age, lifetables$Overall)
  
  death_hazard_NHL <- rnorm(n = n_samples, mean = exp(-2.092), sd = 0.006378) #do I exponentiate the sd?
  death_probability_NHL <-	1-exp(-death_hazard_NHL) #will relate to table above
  
  #During the total study period, the hazard ratio for one-year all-cause mortality was 3.5 times (95% CI: 3.28–3.74) 
  #greater for male hip fracture patients than control subjects and 2.4 times (95% CI: 2.31–2.50) greater than 
  #controls for females. 
  #death_probability_osteoporosis <-	lifetables
  #death_probability_osteoporosis$Males <- death_probability_osteoporosis$Males*3.5   #Currently  not probabilistic as lifetables are not probabilistic
  #death_probability_osteoporosis$Females <- death_probability_osteoporosis$Females*2.4
  #death_probability_osteoporosis$Overall <- (percentage_male * death_probability_osteoporosis$Males) + ((1-percentage_male) * death_probability_osteoporosis$Females)
  #death_probability_osteoporosis	<- data.frame(death_probability_osteoporosis$Age, death_probability_osteoporosis$Overall)
  
  
  
  
  utility_GFD <- 0.85
  utility_GFdse <- ((0.86-0.84)/3.92)
  utility_GFDalpha <- (utility_GFD ^ 2 * (1 - utility_GFD)/utility_GFdse ^ 2) - utility_GFD
  utility_GFDbeta <- (utility_GFDalpha / utility_GFD) - utility_GFDalpha
  utility_GFD <- rbeta(n = n_samples, shape1 = utility_GFDalpha, shape2 = utility_GFDbeta)
  
  utility_undiagnosedCD <-  0.65 
  utility_undiagnosedCD_se <- (0.67 - 0.63)/3.92
  utility_undiagnosedCD_alpha <- (utility_undiagnosedCD ^ 2 * (1 - utility_undiagnosedCD)/utility_undiagnosedCD_se ^ 2) - utility_undiagnosedCD
  utility_undiagnosedCD_beta <- (utility_undiagnosedCD_alpha/utility_undiagnosedCD) - utility_undiagnosedCD_alpha
  utility_undiagnosedCD <- rbeta(n = n_samples, shape1 = utility_undiagnosedCD_alpha, shape2 = utility_undiagnosedCD_beta)
  
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
  
  #disutility_noGFD <-  0.14 
  #disutility_noGFD_se <- (0.31 - (-0.03))/3.92
  #disutility_noGFD_alpha <- (disutility_noGFD ^ 2 * (1 - disutility_noGFD)/disutility_noGFD_se ^ 2) - disutility_noGFD
  #disutility_noGFD_beta <- (disutility_noGFD_alpha/disutility_noGFD) - disutility_noGFD_alpha
  #disutility_noGFD <- rbeta(n = n_samples, shape1 = disutility_noGFD_alpha, shape2 = disutility_noGFD_beta)
  
  
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
  
  cost_undiagnosedCD <- 340
  cost_undiagnosedCD_se <- 2.96
  cost_undiagnosedCD_alpha <- (cost_undiagnosedCD/cost_undiagnosedCD_se)^2
  cost_undiagnosedCD_beta <- (cost_undiagnosedCD_se^2)/cost_undiagnosedCD
  cost_undiagnosedCD <- rgamma(n = n_samples, shape = cost_undiagnosedCD_alpha, scale = cost_undiagnosedCD_beta) 
  
  cost_CDGFD <- 650
  cost_CDGFD_se <- 4.68
  cost_CDGFD_alpha <- (cost_CDGFD/cost_CDGFD_se)^2
  cost_CDGFD_beta <- (cost_CDGFD_se^2)/cost_CDGFD
  cost_CDGFD <- rgamma(n = n_samples, shape = cost_CDGFD_alpha, scale = cost_CDGFD_beta)
  
  cost_biopsy <- 530
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
  
  
  return(data.frame(probability_late_diagnosis, probability_subfertility, probability_osteoporosis, probability_NHL, probability_nocomplications,
                    osteoporosis_probability_GFD_all, subfertility_probability_GFD_all, NHL_probability_GFD, osteoporosis_probability_noGFD_all, subfertility_probability_noGFD_all,
                    NHL_probability_noGFD, death_probability_NHL, 
                    utility_GFD, utility_undiagnosedCD, disutility_subfertility, disutility_osteoporosis, disutility_NHL,
                    cost_CDGFD, cost_osteoporosis, cost_undiagnosedCD, cost_IDA, cost_biopsy, probability_biopsy,
                    cost_subfertility, cost_NHL, probability_IDA, cost_diagnosis))
}

generate_model_parameters(starting_age)