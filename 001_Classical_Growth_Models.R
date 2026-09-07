### ==============================================================================
### SCRIPT 1: Comparison of Classical Growth Curve Models
### Aim: Fitting 7 different non-linear growth models to individual animal data
### ==============================================================================

if(!require(minpack.lm)) install.packages("minpack.lm")
if(!require(dplyr)) install.packages("dplyr")
if(!require(tidyr)) install.packages("tidyr")
if(!require(readxl)) install.packages("readxl")

library(minpack.lm)
library(dplyr)
library(tidyr)
library(readxl)

# 1. Load and Clean Data
df <- read_excel("heifers.xlsx")
colnames(df) <- c("ID", "Day", "BW")
df <- df %>% drop_na()
animals <- unique(df$ID)

results <- data.frame()

# 2. Mathematical Formulas of Growth Models
models <- list(
  Brody = BW ~ A * (1 - B * exp(-k * Day)),
  Gompertz = BW ~ A * exp(-B * exp(-k * Day)),
  Logistic = BW ~ A / (1 + B * exp(-k * Day)),
  VonBertalanffy = BW ~ A * (1 - B * exp(-k * Day))^3,
  Richards = BW ~ A * (1 - B * exp(-k * Day))^M,
  NegExponential = BW ~ C + (A - C) * (1 - exp(-k * Day)),
  Weibull = BW ~ A - B * exp(-k * Day^M)
)

# 3. Biological Starting Values (Priors for NLS)
start_vals <- list(
  Brody = c(A=500, B=1, k=0.01),
  Gompertz = c(A=500, B=1.5, k=0.01),
  Logistic = c(A=500, B=5, k=0.02),
  VonBertalanffy = c(A=500, B=0.5, k=0.01),
  Richards = c(A=500, B=0.5, k=0.01, M=1.5),
  NegExponential = c(A=500, C=150, k=0.01),
  Weibull = c(A=500, B=400, k=0.01, M=0.8)
)

cat("Starting classical model fitting...\n")
pb <- txtProgressBar(min = 0, max = length(animals), style = 3)

# 4. Iteration Over Each Animal and Model
for(i in seq_along(animals)) {
  animal <- animals[i]
  df_sub <- filter(df, ID == animal)
  
  if(nrow(df_sub) < 4) next
  
  for(model_name in names(models)) {
    fit <- tryCatch({
      nlsLM(formula = models[[model_name]], data = df_sub, start = start_vals[[model_name]], control = nls.lm.control(maxiter = 1000))
    }, error = function(e) return(NULL))
    
    if(!is.null(fit)) {
      coefs <- coef(fit)
      
      r2 <- 1 - (sum(residuals(fit)^2) / sum((df_sub$BW - mean(df_sub$BW))^2))
      
      row_data <- data.frame(
        ID = animal, Model = model_name, Status = "Successful",
        Asymptote_A = ifelse("A" %in% names(coefs), coefs["A"], NA),
        B_Param = ifelse("B" %in% names(coefs), coefs["B"], NA),
        k_Rate = ifelse("k" %in% names(coefs), coefs["k"], NA),
        M_Param = ifelse("M" %in% names(coefs), coefs["M"], NA),
        C_Param = ifelse("C" %in% names(coefs), coefs["C"], NA),
        R2 = r2, AIC = AIC(fit), BIC = BIC(fit)
      )
      results <- bind_rows(results, row_data)
    } else {
      row_data <- data.frame(ID = animal, Model = model_name, Status = "Failed", Asymptote_A=NA, B_Param=NA, k_Rate=NA, M_Param=NA, C_Param=NA, R2=NA, AIC=NA, BIC=NA)
      results <- bind_rows(results, row_data)
    }
  }
  setTxtProgressBar(pb, i)
}
close(pb)

# 5. Overall Summary for the Herd
summary_table <- results %>%
  filter(Status == "Successful") %>%
  group_by(Model) %>%
  summarise(
    Successful_Animals = n(),
    Mean_R2 = mean(R2, na.rm=TRUE), SD_R2 = sd(R2, na.rm=TRUE),
    Mean_AIC = mean(AIC, na.rm=TRUE), SD_AIC = sd(AIC, na.rm=TRUE),
    Mean_BIC = mean(BIC, na.rm=TRUE), SD_BIC = sd(BIC, na.rm=TRUE),
    Mean_Asymptote = mean(Asymptote_A, na.rm=TRUE), SD_Asymptote = sd(Asymptote_A, na.rm=TRUE),
    Mean_B = mean(B_Param, na.rm=TRUE), SD_B = sd(B_Param, na.rm=TRUE),
    Mean_k = mean(k_Rate, na.rm=TRUE), SD_k = sd(k_Rate, na.rm=TRUE)
  ) %>%
  arrange(Mean_AIC)

write.csv(results, "01_Classical_Individual_Results.csv", row.names = FALSE)
write.csv(summary_table, "01_Classical_Summary_Table.csv", row.names = FALSE)

cat("\nAnalysis complete! Results are saved as '01_Classical_Individual_Results.csv' and '01_Classical_Summary_Table.csv'.\n")
