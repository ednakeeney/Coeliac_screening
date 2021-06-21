results <- data.frame(output$average_costs, output$average_effects)
results <- results[order(results$output.average_costs),]
for (i in 1:216) {
  results$x[i+1] <- ifelse(results$output.average_effects[i+1] > results$output.average_effects[i], 1,0) }
results <- subset(results, x== 1)
for (i in 1:138) {
  results$x[i+1] <- ifelse(results$output.average_effects[i+1] > results$output.average_effects[i], 1,0) }
results <- subset(results, x== 1)
for (i in 1:101) {
  results$x[i+1] <- ifelse(results$output.average_effects[i+1] > results$output.average_effects[i], 1,0) }
results <- subset(results, x== 1)
nrow(results)
for (i in 1:77) {
  results$x[i+1] <- ifelse(results$output.average_effects[i+1] > results$output.average_effects[i], 1,0) }
results <- subset(results, x== 1)
nrow(results)
for (i in 1:58) {
  results$x[i+1] <- ifelse(results$output.average_effects[i+1] > results$output.average_effects[i], 1,0) }
results <- subset(results, x== 1)
nrow(results)
for (i in 1:41) {
  results$x[i+1] <- ifelse(results$output.average_effects[i+1] > results$output.average_effects[i], 1,0) }
results <- subset(results, x== 1)
nrow(results)
for (i in 1:30) {
  results$x[i+1] <- ifelse(results$output.average_effects[i+1] > results$output.average_effects[i], 1,0) }
results <- subset(results, x== 1)
nrow(results)
for (i in 1:28) {
  results$x[i+1] <- ifelse(results$output.average_effects[i+1] > results$output.average_effects[i], 1,0) }
results <- subset(results, x== 1)
nrow(results)
for (i in 1:27) {
  results$x[i+1] <- ifelse(results$output.average_effects[i+1] > results$output.average_effects[i], 1,0) }
results <- subset(results, x== 1)
nrow(results)

results$ICER <- 0
for (i in 1:26) {
  results$cost_difference[i+1] <- results$output.average_costs[i+1] - results$output.average_costs[i]
  results$effect_difference[i+1] <- results$output.average_effects[i+1] - results$output.average_effects[i]
  # results$ICER[i+1] <- (results$output.average_costs[i+1] - results$output.average_costs[i])/ (results$output.average_effects[i+1] - results$output.average_effects[i]) }
}