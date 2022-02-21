# Script to run GAM or GP for EVPPI
# Needs to be embedded within coeliac_markov_main script and runs after info.rank
# Howard Thom 12-November-2021

require(ggplot2)
require(BCEA)

evpi_table_bcea <- evpi_table_bcea


bcea_evppi_utilities <- evppi(param_idx =  c("utility_undiagnosedCD", 
                                     "disutility_osteoporosis", "disutility_NHL",
                                     "disutility_biopsy", #"disutility_biopsy_wait",
                                     "disutility_fp"), 
                      input = input_parameters_info_rank, he = m_of_interest, method = 'gp')
evpi_table_bcea["Utilities and disutilities", c("Per person", "Population")] <- bcea_evppi_utilities$evppi[201] * c(1, discounted_population_size)
write.csv(evpi_table_bcea, paste0("results/", population, "/", population, "_evppi_bcea_summary.csv"))

print("Rate...")
bcea_evppi_rate <- evppi(param_idx =  colnames(input_parameters_info_rank)[grepl("rate", colnames(input_parameters_info_rank))], 
                              input = input_parameters_info_rank, he = m_of_interest, method = 'gp')
evpi_table_bcea["Rates of osteoporosis and NHL", c("Per person", "Population")] <- bcea_evppi_rate$evppi * c(1, discounted_population_size)
write.csv(evpi_table_bcea, paste0("results/", population, "/", population, "_evppi_bcea_summary.csv"))

print("gfd...")
bcea_evppi_gfd <- evppi(param_idx =  
                          c("utility_GFD",
                            "log_or_osteoporosis_GFD",
                            "log_or_osteoporosis_noGFD", 
                            "log_rr_NHL_GFD" ,
                            "log_rr_NHL_noGFD"), 
                         input = input_parameters_info_rank, he = m_of_interest, method = 'gp')
evpi_table_bcea["GFD effect", c("Per person", "Population")] <- bcea_evppi_gfd$evppi * c(1, discounted_population_size)
write.csv(evpi_table_bcea, paste0("results/", population, "/", population, "_evppi_bcea_summary.csv"))

print("sens...")
bcea_evppi_sens <- evppi(param_idx =  
                           colnames(input_parameters_info_rank)[grepl("sens", colnames(input_parameters_info_rank)) | grepl("spec", colnames(input_parameters))],
                        input = input_parameters_info_rank, he = m_of_interest, method = 'gp')
evpi_table_bcea["Test accuracies", c("Per person", "Population")] <- bcea_evppi_sens$evppi * c(1, discounted_population_size)
write.csv(evpi_table_bcea, paste0("results/", population, "/", population, "_evppi_bcea_summary.csv"))

print("late...")
bcea_evppi_late <- evppi(param_idx =  c("probability_late_diagnosis"),
                         input = input_parameters_info_rank, he = m_of_interest, method = 'gp')
evpi_table_bcea["Probability of late diagnosis", c("Per person", "Population")] <- bcea_evppi_late$evppi * c(1, discounted_population_size)
write.csv(evpi_table_bcea, paste0("results/", population, "/", population, "_evppi_bcea_summary.csv"))



# Export results
write.csv(evpi_table_bcea, paste0("results/", population, "/", population, "_evppi_bcea_summary.csv"))

