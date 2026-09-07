### ==============================================================================
### SCRIPT 3: Calculation of PTI and stAUGC Using the Best Model (Logistic)
### Aim: Extracting the best parameters and calculating the profitability area
### ==============================================================================

library(dplyr)
library(tidyr)

# 1. Load the Results Generated from Script 2
df_results <- read.csv("02_Bayesian_7_Models_Full_Stats.csv")

# 2. Filter the Best Model (Logistic)
# Extracting only the data of the Logistic model with successfully calculated Parameter A
df_champ <- df_results %>% filter(Model == "Logistic" & !is.na(Asymptote_A_kg))

# Mathematical Formulas (Specifically for the Logistic Model)
logistic_func <- function(t, A, B, k) { A / (1 + B * exp(-k * t)) }
find_t <- function(target_W, A, B, k) { -log(((A / target_W) - 1) / B) / k }

results_temp <- data.frame()
cat("Calculating PTI and stAUGC Values...\n")

# 3. Find the Time to Insemination (PTI) for Each Animal
for(i in 1:nrow(df_champ)) {
  h <- df_champ$ID[i]
  A_est <- df_champ$Asymptote_A_kg[i]
  B_est <- df_champ$B_Param[i]
  k_est <- df_champ$k_Rate[i]
  
  # Calculate the time (days) required to reach 228 kg and 375 kg body weights
  t_228 <- find_t(228, A_est, B_est, k_est)
  t_375 <- find_t(375, A_est, B_est, k_est)
  
  # Create a new row combining all statistics
  row_data <- data.frame(
    ID = h,
    Asymptote_A_kg = A_est,
    HPD_Lower_95 = df_champ$HPD_Lower_95[i],
    HPD_Upper_95 = df_champ$HPD_Upper_95[i],
    B_Param = B_est,
    k_Rate = k_est,
    RMSE = df_champ$RMSE[i],
    Bayes_R2 = df_champ$Bayes_R2[i],
    LOOIC = df_champ$LOOIC[i],
    elpd_diff = df_champ$elpd_diff[i],
    se_diff = df_champ$se_diff[i],
    t_228_Days = t_228,
    PTI_Days = t_375
  )
  results_temp <- bind_rows(results_temp, row_data)
}

# 4. Calculation of stAUGC (Synchronized Area Under the Growth Curve)
# Calculate the duration (window) for animals to grow from 228 kg to 375 kg
results_temp$Duration_228_to_375 <- results_temp$PTI_Days - results_temp$t_228_Days

# Find the duration of the fastest animal (The fixed biological time window)
min_duration <- min(results_temp$Duration_228_to_375, na.rm = TRUE)

cat("\nFixed Biological Time Window (Fastest animal):", round(min_duration, 1), "Days\n")

# Determine the Profitability Area (stAUGC) by computing the definite integral
final_results <- results_temp %>% 
  rowwise() %>% 
  mutate(stAUGC = {
    if(!is.na(t_228_Days) && t_228_Days > -999) {
      tryCatch({
        integral <- integrate(logistic_func, 
                              lower = t_228_Days, 
                              upper = t_228_Days + min_duration, 
                              A = Asymptote_A_kg, 
                              B = B_Param, 
                              k = k_Rate)
        integral$value
      }, error = function(e) NA)
    } else { 
      NA 
    }
  }) %>% 
  arrange(desc(stAUGC))

# 5. Rounding and Saving
final_results <- final_results %>% 
  mutate(across(c(t_228_Days, PTI_Days, Duration_228_to_375, stAUGC), ~round(., 1)))

write.csv(final_results, "03_Bayesian_Logistic_stAUGC_Results.csv", row.names = FALSE) 

cat("\nProcess Successful! The stAUGC data has been saved to the '03_Bayesian_Logistic_stAUGC_Results.csv' file.\n")
cat("\n--- TOP 6 MOST PROFITABLE ANIMALS ---\n")
print(head(final_results %>% select(ID, Asymptote_A_kg, PTI_Days, Duration_228_to_375, stAUGC)))
