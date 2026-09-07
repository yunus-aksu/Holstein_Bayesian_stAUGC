### ==============================================================================
### SCRIPT 2: Bayesian Comparison of 7 Models (All Parameters and Statistics)
### Aim: Estimating parameters with HPD intervals, RMSE, Bayes R2, and elpd_diff
### ==============================================================================

library(brms)
library(dplyr)
library(tidyr)
library(readxl)

# 1. Load Data
df <- read_excel("heifers.xlsx")
colnames(df) <- c("ID", "Day", "BW")
df <- df %>% drop_na()
animals <- unique(df$ID)


# 2. Literature-Supported Priors #Huang (2023): 581.76 kg, Perotto et all. (1992): 613.18 kg, Koenen and Groen (1996): 667.0 kg, Bayram et all. (2004): 496.0 kg, Campos-Barreiro et all. (2014): 751.0 kg
mean_A <- 621.8; sd_A <- 95.2
prior_A <- prior_string(paste0("normal(", mean_A, ",", sd_A, ")"), nlpar = "A", lb=450, ub=850)

prior_brody <- c(prior_A, prior(normal(1, 0.5), nlpar = "B", lb=0), prior(normal(0.01, 0.01), nlpar = "k", lb=0))
prior_gomp  <- c(prior_A, prior(normal(3, 2), nlpar = "B", lb=0), prior(normal(0.02, 0.01), nlpar = "k", lb=0))
prior_log   <- c(prior_A, prior(normal(6, 3), nlpar = "B", lb=0), prior(normal(0.02, 0.01), nlpar = "k", lb=0))
prior_vonb  <- c(prior_A, prior(normal(0.6, 0.2), nlpar = "B", lb=0), prior(normal(0.02, 0.01), nlpar = "k", lb=0))
prior_rich  <- c(prior_A, prior(normal(0.5, 0.5), nlpar = "B", lb=0), prior(normal(0.01, 0.01), nlpar = "k", lb=0), prior(normal(1.5, 1), nlpar = "M", lb=0))
prior_neg   <- c(prior_A, prior(normal(150, 50), nlpar = "C", lb=0), prior(normal(0.01, 0.01), nlpar = "k", lb=0))
prior_weib  <- c(prior_A, prior(normal(400, 100), nlpar = "B", lb=0), prior(normal(0.01, 0.01), nlpar = "k", lb=0), prior(normal(0.8, 0.2), nlpar = "M", lb=0))

ctrl <- list(adapt_delta = 0.99, max_treedepth = 15)
iter_cnt <- 2000

results_all <- data.frame()

# Function to extract metrics (Parameters, HPD intervals, RMSE, R2, LOOIC)
get_metrics <- function(fit, model_name, df_sub) {
  if(is.null(fit)) return(NULL)
  f <- fixef(fit)
  
  # Parameter A and 95% HPD (Credible) Intervals
  A_est <- ifelse("A_Intercept" %in% rownames(f), f["A_Intercept", "Estimate"], NA)
  A_low <- ifelse("A_Intercept" %in% rownames(f), f["A_Intercept", "Q2.5"], NA)
  A_up  <- ifelse("A_Intercept" %in% rownames(f), f["A_Intercept", "Q97.5"], NA)
  
  # Other Parameters
  B <- ifelse("B_Intercept" %in% rownames(f), f["B_Intercept", "Estimate"], NA)
  k <- ifelse("k_Intercept" %in% rownames(f), f["k_Intercept", "Estimate"], NA)
  M <- ifelse("M_Intercept" %in% rownames(f), f["M_Intercept", "Estimate"], NA)
  C <- ifelse("C_Intercept" %in% rownames(f), f["C_Intercept", "Estimate"], NA)
  
  # Goodness-of-Fit Metrics
  pred <- predict(fit)[, "Estimate"]
  rmse <- sqrt(mean((df_sub$BW - pred)^2))
  r2 <- bayes_R2(fit)[1, "Estimate"]
  
  loo_obj <- tryCatch({ loo(fit) }, error = function(e) NULL)
  looic <- ifelse(!is.null(loo_obj), loo_obj$estimates["looic", "Estimate"], NA)
  
  list(Model = model_name, Asymptote_A_kg = A_est, HPD_Lower_95 = A_low, HPD_Upper_95 = A_up, 
       B_Param = B, k_Rate = k, M_param = M, C_param = C, 
       RMSE = rmse, Bayes_R2 = r2, LOOIC = looic, loo_obj = loo_obj)
}

cat("Extracting ALL Parameters and MCMC Statistics for the 7 Models...\n")
pb <- txtProgressBar(min = 0, max = length(animals), style = 3)

# 3. Iteration Loop
for(i in seq_along(animals)) {
  h <- animals[i]
  df_sub <- filter(df, ID == h)
  if(nrow(df_sub) < 4) next
  
  # brm Models (seed=12345)
  f_brody <- tryCatch({ brm(bf(BW ~ A * (1 - B * exp(-k * Day)), A + B + k ~ 1, nl=TRUE), data=df_sub, prior=prior_brody, family=gaussian(), chains=2, iter=iter_cnt, cores=2, control=ctrl, silent=2, refresh=0, seed=12345) }, error=function(e) NULL)
  f_gomp  <- tryCatch({ brm(bf(BW ~ A * exp(-B * exp(-k * Day)), A + B + k ~ 1, nl=TRUE), data=df_sub, prior=prior_gomp, family=gaussian(), chains=2, iter=iter_cnt, cores=2, control=ctrl, silent=2, refresh=0, seed=12345) }, error=function(e) NULL)
  f_log   <- tryCatch({ brm(bf(BW ~ A / (1 + B * exp(-k * Day)), A + B + k ~ 1, nl=TRUE), data=df_sub, prior=prior_log, family=gaussian(), chains=2, iter=iter_cnt, cores=2, control=ctrl, silent=2, refresh=0, seed=12345) }, error=function(e) NULL)
  f_vonb  <- tryCatch({ brm(bf(BW ~ A * (1 - B * exp(-k * Day))^3, A + B + k ~ 1, nl=TRUE), data=df_sub, prior=prior_vonb, family=gaussian(), chains=2, iter=iter_cnt, cores=2, control=ctrl, silent=2, refresh=0, seed=12345) }, error=function(e) NULL)
  f_rich  <- tryCatch({ brm(bf(BW ~ A * (1 - B * exp(-k * Day))^M, A + B + k + M ~ 1, nl=TRUE), data=df_sub, prior=prior_rich, family=gaussian(), chains=2, iter=iter_cnt, cores=2, control=ctrl, silent=2, refresh=0, seed=12345) }, error=function(e) NULL)
  f_neg   <- tryCatch({ brm(bf(BW ~ C + (A - C) * (1 - exp(-k * Day)), A + C + k ~ 1, nl=TRUE), data=df_sub, prior=prior_neg, family=gaussian(), chains=2, iter=iter_cnt, cores=2, control=ctrl, silent=2, refresh=0, seed=12345) }, error=function(e) NULL)
  f_weib  <- tryCatch({ brm(bf(BW ~ A - B * exp(-k * Day^M), A + B + k + M ~ 1, nl=TRUE), data=df_sub, prior=prior_weib, family=gaussian(), chains=2, iter=iter_cnt, cores=2, control=ctrl, silent=2, refresh=0, seed=12345) }, error=function(e) NULL)
  
  # Extract Results (Including HPD intervals)
  m_list <- list(
    get_metrics(f_brody, "Brody", df_sub), get_metrics(f_gomp, "Gompertz", df_sub),
    get_metrics(f_log, "Logistic", df_sub), get_metrics(f_vonb, "VonBertalanffy", df_sub),
    get_metrics(f_rich, "Richards", df_sub), get_metrics(f_neg, "NegExponential", df_sub),
    get_metrics(f_weib, "Weibull", df_sub)
  )
  
  # Collect LOO objects of successfully executed models
  m_list <- m_list[!sapply(m_list, is.null)]
  loo_objects <- lapply(m_list, function(x) x$loo_obj)
  names(loo_objects) <- sapply(m_list, function(x) x$Model)
  loo_objects <- loo_objects[!sapply(loo_objects, is.null)]
  
  # Calculate elpd_diff and se_diff (using loo_compare)
  diff_df <- data.frame(Model = character(), elpd_diff = numeric(), se_diff = numeric())
  if(length(loo_objects) > 1) {
    comp <- loo_compare(loo_objects)
    diff_df <- data.frame(Model = rownames(comp), elpd_diff = comp[,"elpd_diff"], se_diff = comp[,"se_diff"])
  } else if (length(loo_objects) == 1) {
    diff_df <- data.frame(Model = names(loo_objects), elpd_diff = 0, se_diff = 0)
  }
  
  # Append to Main Table
  for(m in m_list) {
    diff_info <- diff_df %>% filter(Model == m$Model)
    e_diff <- ifelse(nrow(diff_info) > 0, diff_info$elpd_diff, NA)
    s_diff <- ifelse(nrow(diff_info) > 0, diff_info$se_diff, NA)
    
    row_data <- data.frame(
      ID = h, Model = m$Model, 
      Asymptote_A_kg = m$Asymptote_A_kg, HPD_Lower_95 = m$HPD_Lower_95, HPD_Upper_95 = m$HPD_Upper_95,
      B_Param = m$B_Param, k_Rate = m$k_Rate, M_param = m$M_param, C_param = m$C_param,
      RMSE = m$RMSE, Bayes_R2 = m$Bayes_R2, LOOIC = m$LOOIC, elpd_diff = e_diff, se_diff = s_diff
    )
    results_all <- bind_rows(results_all, row_data)
  }
  setTxtProgressBar(pb, i)
}
close(pb)

# Rounding and Saving
results_all <- results_all %>%
  mutate(across(c(Asymptote_A_kg, HPD_Lower_95, HPD_Upper_95, M_param, C_param, RMSE), ~round(., 2)),
         across(c(B_Param, Bayes_R2, LOOIC, elpd_diff, se_diff), ~round(., 4)),
         across(c(k_Rate), ~round(., 5)))

write.csv(results_all, "02_Bayesian_7_Models_Full_Stats.csv", row.names = FALSE)
cat("\n\nAnalysis complete! All parameters are saved in the '02_Bayesian_7_Models_Full_Stats.csv' file.\n")
