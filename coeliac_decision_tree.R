# Coeliac screening decision tree

# HT: I added extra spaces after commas and equal signs to align with Tidyverse
# HT: I added a set.seed() at the start. Ultimately need one at the start of your 'coeliac_main' script.
# HT: There are a lot of hard coded numbers in this script. Ultimately need them in a separate
# R script (e.g. and coming from a generate_input_parameters() function) or Excel file. 

set.seed(26452345)

# Number of PSA samples
n_samples <- 100

# Number and names of strategies
n_treat <- 3
t_names <- c("Test", "Test + biopsy", "Double test")


#prevalence of coeliac disease
p_cd <- 0.5 

# Accuracy of test or test + biopsy
tp <- fn <- fp <- tn <- matrix(nrow=n_samples, ncol=n_treat)

#IgA EMA sensitivity in adults: 88.0 (75.2, 94.7)
sens_test <- 0.88
sens_se <- (0.947 - 0.752) / 3.92
sens_alpha <- (sens_test ^ 2 * (1 - sens_test)/sens_se ^ 2) - sens_test
sens_beta <- (sens_alpha / sens_test) - sens_alpha
#IgA EMA specificity in adults: 99.6 (92.3, 100.0)
spec_test <- 0.996
spec_se <- (1 - 0.923) / 3.92 
spec_alpha <- (spec_test ^ 2 * (1 - spec_test) / spec_se ^ 2) - spec_test
spec_beta <- (spec_alpha / spec_test) - spec_alpha

sens_testbiopsy <- 1
spec_testbiopsy <- 1
sens_doubletest <- 1
spec_doubletest <- 1

# Probabilities for test 
tp[, 1] <- (n_samples * p_cd * rbeta(n = n_samples, shape1 = sens_alpha, shape2 = sens_beta)) / n_samples
fn[, 1] <- ((n_samples * p_cd) - tp[, 1]) / n_samples  
tn[, 1] <- (n_samples * p_cd * rbeta(n = n_samples, shape1 = spec_alpha, shape2 = spec_beta)) / n_samples
fp[, 1] <- ((n_samples * p_cd) - tn[, 1]) / n_samples

# Probabilities for test + biopsy
tp[, 2] <- (n_samples * p_cd * sens_testbiopsy) / n_samples
fn[, 2] <- ((n_samples * p_cd) - tp[, 2]) / n_samples
tn[, 2] <- (n_samples * p_cd * spec_testbiopsy)/n_samples
fp[, 2] <- ((n_samples * p_cd) - tn[, 2]) / n_samples

# Probabilities for Double test
tp[, 3] <- (n_samples * p_cd * sens_doubletest) / n_samples
fn[, 3] <- ((n_samples * p_cd) - tp[, 3]) / n_samples
tn[, 3] <- (n_samples * p_cd * spec_doubletest) / n_samples
fp[, 3] <- ((n_samples * p_cd) - tn[, 3]) / n_samples


# Costs for FP - Cost of GFP for one year
c_fp <- rnorm(n = n_samples, mean = 1000, sd = 50)

# Cost of test, test + biopsy or double test
c_test <- t(matrix(rep(c(300, 500, 400), n_samples), ncol = n_samples, nrow = 3))


