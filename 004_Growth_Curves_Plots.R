### ==============================================================================
### SCRIPT 4: Advanced Publication Plots
### Aim: Bayesian vs Classical comparison, Spaghetti Plot, PTI, and stAUGC concept
### ==============================================================================

if(!require(ggplot2)) install.packages("ggplot2")
if(!require(dplyr)) install.packages("dplyr")
if(!require(minpack.lm)) install.packages("minpack.lm")
if(!require(readxl)) install.packages("readxl")

library(ggplot2)
library(dplyr)
library(minpack.lm)
library(readxl)

# 1. Load Raw Data and Bayesian Results
df <- read_excel("heifers.xlsx")
colnames(df) <- c("ID", "Day", "BW")

# Load the final results from Script 3
final_results <- read.csv("03_Bayesian_Logistic_stAUGC_Results.csv")

cat("Generating publication plots...\n")

### ==============================================================================
### FIGURE 1: CLASSICAL vs BAYESIAN COMPARISON (The Failure of Classical Models)
### Animal ID: 3636
### ==============================================================================
df_3636 <- df %>% filter(ID == 3636)

# Classical Logistic Model (Fails / Overestimates)
fit_class_log <- tryCatch({
  nlsLM(BW ~ A / (1 + B * exp(-k * Day)), data = df_3636, start = list(A=500, B=5, k=0.02), control = nls.lm.control(maxiter = 1000))
}, error=function(e) NULL)
coef_log <- coef(fit_class_log)
classical_log_curve <- function(t) { coef_log["A"] / (1 + coef_log["B"] * exp(-coef_log["k"] * t)) }

# Classical Gompertz Model (Fails / Overestimates)
fit_class_gomp <- tryCatch({
  nlsLM(BW ~ A * exp(-B * exp(-k * Day)), data = df_3636, start = list(A=500, B=1.5, k=0.01), control = nls.lm.control(maxiter = 1000))
}, error=function(e) NULL)
coef_gomp <- coef(fit_class_gomp)
classical_gomp_curve <- function(t) { coef_gomp["A"] * exp(-coef_gomp["B"] * exp(-coef_gomp["k"] * t)) }

# Bayesian Logistic Model (Best)
bayes_3636 <- final_results %>% filter(ID == 3636)
bayesian_log <- function(t) { bayes_3636$Asymptote_A_kg / (1 + bayes_3636$B_Param * exp(-bayes_3636$k_Rate * t)) }
bayesian_upper <- function(t) { bayes_3636$HPD_Upper_95 / (1 + bayes_3636$B_Param * exp(-bayes_3636$k_Rate * t)) }
bayesian_lower <- function(t) { bayes_3636$HPD_Lower_95 / (1 + bayes_3636$B_Param * exp(-bayes_3636$k_Rate * t)) }

t_seq1 <- seq(0, 1500, length.out = 200)
df_fig1 <- data.frame(
  Day = t_seq1, 
  Classical_Log = classical_log_curve(t_seq1), 
  Classical_Gomp = classical_gomp_curve(t_seq1), 
  Bayesian = bayesian_log(t_seq1), 
  Upper = bayesian_upper(t_seq1), 
  Lower = bayesian_lower(t_seq1)
)

fig1 <- ggplot(df_fig1, aes(x = Day)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = Bayesian, color = "Bayesian Logistic (95% HPD)"), linewidth = 1.2) +
  geom_line(aes(y = Classical_Log, color = "Classical Logistic"), linewidth = 1.2, linetype = "dashed") +
  geom_line(aes(y = Classical_Gomp, color = "Classical Gompertz"), linewidth = 1.2, linetype = "dotdash") +
  geom_point(data = df_3636, aes(x = Day, y = BW), color = "black", size = 2.5) +
  scale_color_manual(name = "Model Fit", values = c("Bayesian Logistic (95% HPD)" = "blue", "Classical Logistic" = "red", "Classical Gompertz" = "green")) +
  coord_cartesian(xlim = c(0, 1500), ylim = c(0, 3500)) +
  labs(title = "Bayesian 95% HPD vs. Classical Models (ID: 3636)", x = "Relative Age (Days)", y = "Body Weight (kg)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("Figure_1_Classical_vs_Bayesian.png", plot = fig1, width = 8, height = 6, dpi = 300)

### ==============================================================================
### FIGURE 2: SPAGHETTI PLOT (Individual Growth Trajectories vs Population)
### ==============================================================================
t_seq2 <- seq(0, 1500, length.out = 100)
df_spaghetti <- data.frame()

for(i in 1:nrow(final_results)) {
  w_vals <- final_results$Asymptote_A_kg[i] / (1 + final_results$B_Param[i] * exp(-final_results$k_Rate[i] * t_seq2))
  df_spaghetti <- rbind(df_spaghetti, data.frame(Day = t_seq2, BW = w_vals, ID = as.factor(final_results$ID[i])))
}

df_mean <- df_spaghetti %>%
  group_by(Day) %>%
  summarise(BW = mean(BW, na.rm = TRUE), ID = "Population Mean")

fig2 <- ggplot() +
  geom_line(data = df_spaghetti, aes(x = Day, y = BW, group = ID), color = "gray70", alpha = 0.6, linewidth = 0.8) +
  geom_line(data = df_mean, aes(x = Day, y = BW, color = "Population Mean"), linewidth = 2) +
  scale_color_manual(name = "", values = c("Population Mean" = "red")) +
  labs(title = "Individual Growth Trajectories vs. Population Mean", x = "Relative Age (Days)", y = "Body Weight (kg)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("Figure_2_Spaghetti_Plot.png", plot = fig2, width = 8, height = 6, dpi = 300)

### ==============================================================================
### FIGURE 3: PTI COMPARISON (Fast vs Slow Maturing)
### Animal IDs: 3798 (Fast) vs 4998 (Slow)
### ==============================================================================
fast_pti <- final_results %>% filter(ID == 3798)
slow_pti <- final_results %>% filter(ID == 4998)

fast_curve <- function(t) { fast_pti$Asymptote_A_kg / (1 + fast_pti$B_Param * exp(-fast_pti$k_Rate * t)) }
slow_curve <- function(t) { slow_pti$Asymptote_A_kg / (1 + slow_pti$B_Param * exp(-slow_pti$k_Rate * t)) }

fig3 <- ggplot(data.frame(Day = c(0, 450)), aes(x = Day)) +
  stat_function(fun = fast_curve, aes(color = "Fast Maturing (ID: 3798)"), linewidth = 1.2) +
  stat_function(fun = slow_curve, aes(color = "Slow Maturing (ID: 4998)"), linewidth = 1.2) +
  geom_hline(yintercept = 375, linetype = "dashed", color = "blue", linewidth = 1) +
  annotate("text", x = 10, y = 390, label = "Target Insemination Weight (375 kg)", color = "blue", fontface = "bold", hjust = 0) +
  annotate("segment", x = fast_pti$PTI_Days, y = 0, xend = fast_pti$PTI_Days, yend = 375, linetype = "dotted", color = "forestgreen", linewidth = 1) +
  annotate("segment", x = slow_pti$PTI_Days, y = 0, xend = slow_pti$PTI_Days, yend = 375, linetype = "dotted", color = "firebrick", linewidth = 1) +
  scale_color_manual(name = "Animal Profiles", values = c("Fast Maturing (ID: 3798)" = "forestgreen", "Slow Maturing (ID: 4998)" = "firebrick")) +
  annotate("text", x = fast_pti$PTI_Days + 15, y = 100, label = paste0("PTI: ", round(fast_pti$PTI_Days,1), " days"), color = "forestgreen", fontface="bold") +
  annotate("text", x = slow_pti$PTI_Days + 15, y = 100, label = paste0("PTI: ", round(slow_pti$PTI_Days,1), " days"), color = "firebrick", fontface="bold") +
  labs(title = "Predicted Time to Insemination (PTI)", x = "Relative Age (Days)", y = "Body Weight (kg)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("Figure_3_PTI_Comparison.png", plot = fig3, width = 8, height = 6, dpi = 300)

### ==============================================================================
### FIGURE 4: stAUGC COMPARISON (Synchronized Profitability Area)
### Fixed biological time window = min_duration (Fastest animal)
### Animal IDs: 3651 (High stAUGC) vs 3638 (Low stAUGC)
### ==============================================================================
high_staugc <- final_results %>% filter(ID == 3651)
low_staugc <- final_results %>% filter(ID == 3638)
min_duration <- min(final_results$Duration_228_to_375, na.rm = TRUE)

t_sync <- seq(0, min_duration, length.out = 100)

st_fast_BW <- high_staugc$Asymptote_A_kg / (1 + high_staugc$B_Param * exp(-high_staugc$k_Rate * (t_sync + high_staugc$t_228_Days)))
st_slow_BW <- low_staugc$Asymptote_A_kg / (1 + low_staugc$B_Param * exp(-low_staugc$k_Rate * (t_sync + low_staugc$t_228_Days)))

df_stAUGC <- rbind(
  data.frame(SyncTime = t_sync, BW = st_fast_BW, Animal = "High stAUGC (ID: 3651)"),
  data.frame(SyncTime = t_sync, BW = st_slow_BW, Animal = "Low stAUGC (ID: 3638)")
)

fig4 <- ggplot(df_stAUGC, aes(x = SyncTime, y = BW, fill = Animal, color = Animal)) +
  geom_area(alpha = 0.4, position = "identity") +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 375, linetype = "dashed", color = "blue", linewidth = 1) +
  geom_hline(yintercept = 228, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_vline(xintercept = min_duration, linetype = "dotted", color = "black", linewidth = 1) +
  scale_x_continuous(breaks = c(0, 50, 100, round(min_duration, 1))) +
  scale_fill_manual(name = "Profitability Profile", values = c("High stAUGC (ID: 3651)" = "forestgreen", "Low stAUGC (ID: 3638)" = "firebrick")) +
  scale_color_manual(name = "Profitability Profile", values = c("High stAUGC (ID: 3651)" = "darkgreen", "Low stAUGC (ID: 3638)" = "darkred")) +
  annotate("text", x = 2, y = 385, label = "Target Insemination Weight (375 kg)", color = "blue", fontface = "bold", hjust = 0) +
  annotate("text", x = 2, y = 233, label = "Synchronization Baseline (228 kg)", size=4, hjust = 0) +
  labs(title = "Synchronized Truncated AUGC (stAUGC)", x = "Synchronized Biological Time Window (Days)", y = "Body Weight (kg)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("Figure_4_stAUGC_Comparison.png", plot = fig4, width = 8, height = 6, dpi = 300)

cat("\nAll 4 plots have been successfully saved\n")
