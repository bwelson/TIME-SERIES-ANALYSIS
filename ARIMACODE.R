# ============================================================================
# COMPREHENSIVE CRUDE OIL PRICE ANALYSIS USING ARIMA
# Structured Workflow with Proper Diagnostics
# Author: [BENTUM WELSON]
# Date: November 2023 
# ============================================================================

# ============================================================================
# SECTION 1: SETUP & LIBRARIES
# ============================================================================

# Clear environment
rm(list = ls())
gc()

# Load required libraries
library(forecast)      # ARIMA modeling
library(tseries)       # Stationarity tests (ADF, KPSS, PP)
library(lmtest)        # Ljung-Box test
library(FinTS)         # ARCH test
library(tidyverse)     # Data manipulation
library(ggplot2)       # Visualization
library(gridExtra)     # Multiple plots
library(knitr)         # Reporting
library(seastests)     # Seasonality tests

 Set working directory (adjust to your path)
setwd("~/your/project/folder")

# Create output directories
dir.create("output", showWarnings = FALSE)
dir.create("output/figures", showWarnings = FALSE)
dir.create("output/tables", showWarnings = FALSE)
dir.create("output/models", showWarnings = FALSE)

# ============================================================================
# SECTION 2: DATA LOADING & PREPROCESSING
# ============================================================================

cat("\n========== SECTION 2: DATA PREPROCESSING ==========\n")
# Load data
crude <- read.csv("crude.csv")
cat("✓ Data loaded:", nrow(crude), "rows\n")

# Convert to long format (If Needed)
crude_long <- crude %>%
  pivot_longer(cols = Jan:Dec, names_to = "Month", values_to = "Price") %>%
  mutate(
    Month_Num = match(Month, month.abb),
    Date = as.Date(paste(Year, Month_Num, "01", sep = "-")),
    Price = ifelse(Price == 0, NA, Price)  # Replace 0 with NA
  ) %>%
  arrange(Date) %>%
  filter(!is.na(Price))  # Remove missing values

at("✓ Data cleaned:", nrow(crude_long), "observations retained\n")
cat("  Date range:", as.character(min(crude_long$Date)), "to", 
    as.character(max(crude_long$Date)), "\n")

# Create time series object
crude_ts <- ts(crude_long$Price, 
               start = c(crude_long$Year[1], crude_long$Month_Num[1]), 
               frequency = 12)

# Summary statistics
cat("\n--- Summary Statistics ---\n")
print(summary(crude_ts))
cat("SD:", round(sd(crude_ts), 2), "\n")
cat("CV:", round(sd(crude_ts)/mean(crude_ts)*100, 2), "%\n")


# ============================================================================
# SECTION 3: EXPLORATORY DATA ANALYSIS
# ============================================================================

cat("\n========== SECTION 3: EXPLORATORY ANALYSIS ==========\n")

# Time series plot
png("output/figures/01_timeseries_original.png", width = 1200, height = 600)
plot(crude_ts, 
     main = "Monthly Average Brent Crude Oil Prices (1999-2023)",
     ylab = "Price (USD/barrel)",
     xlab = "Year",
     col = "darkblue",
     lwd = 2)
grid()
dev.off()
cat("✓ Time series plot saved\n")

# Histogram and Q-Q plot
png("output/figures/02_distribution.png", width = 1200, height = 600)
par(mfrow = c(1, 2))
hist(crude_ts, breaks = 30, col = "steelblue", 
     main = "Distribution of Prices",
     xlab = "Price (USD/barrel)")
qqnorm(crude_ts, main = "Q-Q Plot")
qqline(crude_ts, col = "red", lwd = 2)
par(mfrow = c(1, 1))
dev.off()
cat("✓ Distribution plots saved\n")


# ============================================================================
# SECTION 4: COMPREHENSIVE STATIONARITY TESTING
# ============================================================================

cat("\n========== SECTION 4: STATIONARITY TESTING ==========\n")

# Function for comprehensive stationarity tests
test_stationarity <- function(ts_data, name, alpha = 0.05) {
  cat("\n--- Testing:", name, "---\n")
  
  # 1. ADF Test (H0: Unit root/non-stationary)
  adf <- adf.test(ts_data)
  
  # 2. PP Test (Phillips-Perron - robust to heteroskedasticity)
  pp <- PP.test(ts_data)
  
  # 3. KPSS Test (H0: Stationary)
  kpss <- kpss.test(ts_data)
  
  results <- data.frame(
    Test = c("ADF", "PP", "KPSS"),
    Statistic = c(adf$statistic, pp$statistic, kpss$statistic),
    P_Value = c(adf$p.value, pp$p.value, kpss$p.value),
    Conclusion = c(
      ifelse(adf$p.value < alpha, "Stationary", "Non-stationary"),
      ifelse(pp$p.value < alpha, "Stationary", "Non-stationary"),
      ifelse(kpss$p.value > alpha, "Stationary", "Non-stationary")
    )
  )
  
  print(results)
  
  # Consensus decision
  stationary_count <- sum(results$Conclusion == "Stationary")
  
  if(stationary_count >= 2) {
    cat("\n✓ CONSENSUS: Series is STATIONARY (", stationary_count, "/3 tests)\n")
    return(list(stationary = TRUE, data = ts_data, d = 0))
  } else {
    cat("\n✗ CONSENSUS: Series is NON-STATIONARY (", 3-stationary_count, "/3 tests)\n")
    cat("   → Applying first differencing...\n")
    
    ts_diff <- diff(ts_data)
    
    # Re-test differenced series
    adf_diff <- adf.test(ts_diff)
    pp_diff <- PP.test(ts_diff)
    kpss_diff <- kpss.test(ts_diff)
    
    results_diff <- data.frame(
      Test = c("ADF", "PP", "KPSS"),
      Statistic = c(adf_diff$statistic, pp_diff$statistic, kpss_diff$statistic),
      P_Value = c(adf_diff$p.value, pp_diff$p.value, kpss_diff$p.value),
      Conclusion = c(
        ifelse(adf_diff$p.value < alpha, "Stationary", "Non-stationary"),
        ifelse(pp_diff$p.value < alpha, "Stationary", "Non-stationary"),
        ifelse(kpss_diff$p.value > alpha, "Stationary", "Non-stationary")
      )
    )
    
    cat("\n   Differenced Series Results:\n")
    print(results_diff)
    
    stationary_count_diff <- sum(results_diff$Conclusion == "Stationary")
    
    if(stationary_count_diff >= 2) {
      cat("\n✓ Differenced series is STATIONARY (", stationary_count_diff, "/3 tests)\n")
      return(list(stationary = TRUE, data = ts_diff, d = 1, 
                  original_results = results, diff_results = results_diff))
    } else {
      cat("\n⚠ Warning: Series may require d=2 or have structural breaks\n")
      return(list(stationary = FALSE, data = ts_diff, d = 1,
                  original_results = results, diff_results = results_diff))
    }
  }
}

# Test original series
stationarity_result <- test_stationarity(crude_ts, "Original Series")

# Save results
stationarity_df <- if(stationarity_result$d == 0) {
  stationarity_result$original_results
} else {
  rbind(
    cbind(Series = "Original", stationarity_result$original_results),
    cbind(Series = "Differenced", stationarity_result$diff_results)
  )
}

write.csv(stationarity_df, "output/tables/stationarity_tests.csv", row.names = FALSE)
cat("\n✓ Stationarity test results saved\n")

# ============================================================================
# SECTION 5: ACF/PACF ANALYSIS
# ============================================================================

cat("\n========== SECTION 5: ACF/PACF ANALYSIS ==========\n")

png("output/figures/03_acf_pacf_comparison.png", width = 1200, height = 800)
par(mfrow = c(2, 2))

# Original series
acf(crude_ts, main = "ACF - Original Series", lag.max = 48)
pacf(crude_ts, main = "PACF - Original Series", lag.max = 48)

# Differenced series (if applicable)
if(stationarity_result$d > 0) {
  acf(stationarity_result$data, main = "ACF - Differenced Series", lag.max = 48)
  pacf(stationarity_result$data, main = "PACF - Differenced Series", lag.max = 48)
}

par(mfrow = c(1, 1))
dev.off()
cat("✓ ACF/PACF plots saved\n")

# ============================================================================
# SECTION 6: SEASONALITY TESTING
# ============================================================================

cat("\n========== SECTION 6: SEASONALITY TESTING ==========\n")

# Visual seasonality check
png("output/figures/04_seasonality_check.png", width = 1200, height = 800)
par(mfrow = c(2, 1))
monthplot(crude_ts, main = "Seasonal Subseries Plot", ylab = "Price", col = "darkblue")
boxplot(Price ~ Month, data = crude_long, 
        main = "Price Distribution by Month",
        xlab = "Month", ylab = "Price (USD/barrel)",
        col = rainbow(12), las = 2)
par(mfrow = c(1, 1))
dev.off()

# --- Formal Statistical Seasonality Tests ---
cat("\n--- Formal Seasonality Tests ---\n")

qs_result <- qs(crude_ts)
combined_result <- combined_test(crude_ts)

# Extract single p-values
qs_p <- qs_result$Pval
combined_p <- min(combined_result$Pval, na.rm = TRUE)  # handle multiple values

cat("QS Test p-value:", round(qs_p, 4),
    ifelse(qs_p < 0.05, "→ Seasonal ✓", "→ Not seasonal ✗"), "\n")

cat("Combined Test p-value(s):", paste(round(combined_result$Pval, 4), collapse = ", "),
    "\nMinimum p-value:", round(combined_p, 4),
    ifelse(combined_p < 0.05, "→ Seasonal ✓", "→ Not seasonal ✗"), "\n")

# Logical OR for final decision
seasonality_present <- (qs_p < 0.05) || (combined_p < 0.05)

cat("\n✓ CONCLUSION:",
    ifelse(seasonality_present, "Seasonality detected", "No significant seasonality"), "\n")


# ============================================================================
# SECTION 7: AUTOMATIC MODEL SELECTION
# ============================================================================

cat("\n========== SECTION 7: MODEL SELECTION ==========\n")

# Fit auto.arima (no constraints)
cat("\n--- Running auto.arima (exhaustive search) ---\n")
auto_model <- auto.arima(crude_ts,
                         seasonal = seasonality_present,
                         stepwise = FALSE,
                         approximation = FALSE,
                         trace = TRUE)
cat("\n✓ Selected Model:\n")
print(summary(auto_model))

# Extract model specification
model_order <- arimaorder(auto_model)
cat("\nModel Specification: ARIMA(", model_order[1], ",", 
    model_order[2], ",", model_order[3], ")", sep = "")
if(seasonality_present) {
  cat("(", model_order[4], ",", model_order[5], ",", 
      model_order[6], ")[", model_order[7], "]", sep = "")
}
cat("\n")

# ============================================================================
# SECTION 8: COMPREHENSIVE RESIDUAL DIAGNOSTICS
# ============================================================================

cat("\n========== SECTION 8: RESIDUAL DIAGNOSTICS ==========\n")

# Extract residuals
residuals_model <- residuals(auto_model)

# Create diagnostic results table
diagnostic_results <- data.frame(
  Test = character(),
  Statistic = numeric(),
  P_Value = numeric(),
  Conclusion = character(),
  stringsAsFactors = FALSE
)

# 1. Ljung-Box Test (Autocorrelation)
lb_test <- Box.test(residuals_model, lag = 20, type = "Ljung-Box", fitdf = sum(model_order[1:3]))
diagnostic_results <- rbind(diagnostic_results, data.frame(
  Test = "Ljung-Box (Autocorrelation)",
  Statistic = as.numeric(lb_test$statistic),
  P_Value = lb_test$p.value,
  Conclusion = ifelse(lb_test$p.value > 0.05, "✓ No autocorrelation", "✗ Autocorrelation present")
))

# 2. Shapiro-Wilk Test (Normality)
sw_test <- shapiro.test(residuals_model)
diagnostic_results <- rbind(diagnostic_results, data.frame(
  Test = "Shapiro-Wilk (Normality)",
  Statistic = as.numeric(sw_test$statistic),
  P_Value = sw_test$p.value,
  Conclusion = ifelse(sw_test$p.value > 0.05, "✓ Normal residuals", "✗ Non-normal residuals")
))

# 3. ARCH Test (Heteroskedasticity)
arch_test <- ArchTest(residuals_model, lags = 20)
diagnostic_results <- rbind(diagnostic_results, data.frame(
  Test = "ARCH (Heteroskedasticity)",
  Statistic = as.numeric(arch_test$statistic),
  P_Value = arch_test$p.value,
  Conclusion = ifelse(arch_test$p.value > 0.05, "✓ No ARCH effects", "✗ ARCH effects present")
))

# 4. Jarque-Bera Test (Alternative normality test)
jb_test <- jarque.bera.test(residuals_model)
diagnostic_results <- rbind(diagnostic_results, data.frame(
  Test = "Jarque-Bera (Normality)",
  Statistic = as.numeric(jb_test$statistic),
  P_Value = jb_test$p.value,
  Conclusion = ifelse(jb_test$p.value > 0.05, "✓ Normal residuals", "✗ Non-normal residuals")
))

cat("\n--- Residual Diagnostic Results ---\n")
print(diagnostic_results)

# Save diagnostic results
write.csv(diagnostic_results, "output/tables/residual_diagnostics.csv", row.names = FALSE)

# Visualize diagnostics
png("output/figures/05_residual_diagnostics.png", width = 1400, height = 1000)
par(mfrow = c(2, 2))

# Time plot of residuals
plot(residuals_model, main = "Residuals Over Time", 
     ylab = "Residuals", type = "l", col = "darkblue")
abline(h = 0, col = "red", lty = 2)

# ACF of residuals
acf(residuals_model, main = "ACF of Residuals", lag.max = 40)

# Histogram
hist(residuals_model, breaks = 30, col = "steelblue",
     main = "Distribution of Residuals", xlab = "Residuals", probability = TRUE)
lines(density(residuals_model), col = "red", lwd = 2)
curve(dnorm(x, mean = mean(residuals_model), sd = sd(residuals_model)), 
      add = TRUE, col = "blue", lwd = 2, lty = 2)
legend("topright", c("Actual", "Normal"), col = c("red", "blue"), 
       lty = c(1, 2), lwd = 2)

# Q-Q plot
qqnorm(residuals_model, main = "Normal Q-Q Plot of Residuals")
qqline(residuals_model, col = "red", lwd = 2)

par(mfrow = c(1, 1))
dev.off()
cat("\n✓ Diagnostic plots saved\n")

# ============================================================================
# SECTION 9: INTERPRETATION & RECOMMENDATIONS
# ============================================================================

cat("\n========== SECTION 9: MODEL EVALUATION ==========\n")

# Critical test: Ljung-Box
if(diagnostic_results$P_Value[1] > 0.05) {
  cat("\n✓ MODEL IS ADEQUATE:\n")
  cat("  - No autocorrelation in residuals (Ljung-Box p =", 
      round(diagnostic_results$P_Value[1], 4), ")\n")
  cat("  - Model captures temporal structure\n")
  model_adequate <- TRUE
} else {
  cat("\n✗ MODEL IS INADEQUATE:\n")
  cat("  - Autocorrelation remains in residuals (Ljung-Box p =", 
      round(diagnostic_results$P_Value[1], 4), ")\n")
  cat("  - Consider: Higher order ARIMA, seasonal terms, or structural breaks\n")
  model_adequate <- FALSE
}

# Non-normality handling
if(diagnostic_results$P_Value[2] < 0.05 | diagnostic_results$P_Value[4] < 0.05) {
  cat("\n⚠ NON-NORMAL RESIDUALS DETECTED:\n")
  cat("  Impact Assessment:\n")
  cat("  • Point forecasts: ✓ Still valid (ARIMA is robust)\n")
  cat("  • Confidence intervals: ⚠ May be unreliable\n")
  cat("  • Hypothesis tests: ⚠ Use with caution\n")
  cat("\n  Recommendations:\n")
  cat("  1. Use bootstrapped confidence intervals\n")
  cat("  2. Report prediction intervals with caveat\n")
  cat("  3. Consider: GARCH for volatility, or robust methods\n")
  cat("  4. Note: Non-normality is COMMON in commodity prices\n")
}

# ARCH effects
if(diagnostic_results$P_Value[3] < 0.05) {
  cat("\n⚠ ARCH EFFECTS DETECTED:\n")
  cat("  • Volatility clustering present\n")
  cat("  • Recommendation: Use ARIMA-GARCH model for better forecasts\n")
}

# ============================================================================
# SECTION 10: TRAIN-TEST VALIDATION
# ============================================================================

cat("\n========== SECTION 10: OUT-OF-SAMPLE VALIDATION ==========\n")


# Create train-test split (80-20)
train_size <- floor(0.958 * length(crude_ts))
test_size <- length(crude_ts) - train_size

train_ts <- window(crude_ts, end = time(crude_ts)[train_size])
test_ts <- window(crude_ts, start = time(crude_ts)[train_size + 1])

cat("Training set:", train_size, "observations\n")
cat("Test set:", test_size, "observations\n")
cat("Split date:", time(crude_ts)[train_size], "\n")

# Fit model on training data
model_train <- auto.arima(train_ts, seasonal = seasonality_present)

# Generate forecasts
png("output/figures/11_forecast_original.png", width = 1400, height = 800)
forecast_test <- forecast(model_train, h = test_size)
print(plot(forecast_test))

grid()
dev.off()
cat("\n✓ Forecast plot saved\n")

# Calculate accuracy
accuracy_metrics <- accuracy(forecast_test, test_ts)
cat("\n--- Forecast Accuracy ---\n")
print(accuracy_metrics)

# Save accuracy results
write.csv(accuracy_metrics, "output/tables/forecast_accuracy.csv")

# Visualize forecasts
png("output/figures/06_forecast_validation.png", width = 1400, height = 800)
plot(crude_ts, 
     xlim = c(start(train_ts)[1], end(test_ts)[1]),
     ylim = c(min(crude_ts) * 0.9, max(crude_ts) * 1.1),
     main = "Out-of-Sample Forecast Validation",
     ylab = "Price (USD/barrel)",
     xlab = "Year",
     lwd = 2)

abline(v = time(crude_ts)[train_size], col = "gray30", lty = 2, lwd = 2)
text(time(crude_ts)[train_size], max(crude_ts) * 1.05, 
     "Train | Test", pos = 4, col = "gray30")

lines(forecast_test$mean, col = "blue", lwd = 2)
polygon(c(time(forecast_test$mean), rev(time(forecast_test$mean))),
        c(forecast_test$lower[,2], rev(forecast_test$upper[,2])),
        col = rgb(0, 0, 1, 0.2), border = NA)

legend("topleft",
       legend = c("Actual", "Forecast", "95% CI"),
       col = c("black", "blue", "blue"),
       lty = c(1, 1, 1),
       lwd = c(2, 2, 10))
grid()
dev.off()
cat("\n✓ Forecast2 plot saved\n")

#MODEL COMPARISON
data_plot <- data.frame(
  Time = time(crude_ts),
  Actual = as.numeric(crude_ts),
  Fitted = as.numeric(fitted(auto_model))
)

png("output/figures/07_actual_vs_fitted.png", width = 1400, height = 800)
ggplot(data_plot, aes(x = Time)) +
  geom_line(aes(y = Actual, color = "Actual")) +
  geom_line(aes(y = Fitted, color = "Fitted")) +
  scale_color_manual(values = c("Actual" = "black", "Fitted" = "blue")) +
  labs(title = "Actual vs Fitted ARIMA Model",
       y = "Output", x = "Time", color = "Series") +
  theme_minimal()
dev.off()
cat("\n✓ actual vs fitted\n")


# ============================================================================
# SECTION 11: FINAL MODEL & SAVE
# ============================================================================

cat("\n========== SECTION 11: FINAL MODEL ==========\n")

# Save final model
saveRDS(auto_model, "output/models/final_arima_model.rds")
cat("✓ Model saved to: output/models/final_arima_model.rds\n")

# Create summary report
summary_report <- list(
  Model_Specification = paste0("ARIMA(", paste(model_order[1:3], collapse = ","), ")"),
  Data_Range = paste(start(crude_ts), "to", end(crude_ts)),
  N_Observations = length(crude_ts),
  Stationarity = ifelse(stationarity_result$d > 0, 
                        paste("Differencing required (d =", stationarity_result$d, ")"),
                        "Stationary"),
  Seasonality = ifelse(seasonality_present, "Present", "Not detected"),
  Model_Adequate = model_adequate,
  Test_RMSE = accuracy_metrics[2, "RMSE"],
  Test_MAPE = accuracy_metrics[2, "MAPE"],
  Residual_Issues = paste(
    ifelse(diagnostic_results$P_Value[1] < 0.05, "Autocorrelation", ""),
    ifelse(diagnostic_results$P_Value[2] < 0.05, "Non-normality", ""),
    ifelse(diagnostic_results$P_Value[3] < 0.05, "ARCH effects", "")
  )
)

# Save summary
capture.output(print(summary_report), file = "output/tables/model_summary.txt")

cat("\n========== ANALYSIS COMPLETE ==========\n")
cat("\nAll outputs saved to 'output/' directory:\n")
cat("  • Figures: output/figures/\n")
cat("  • Tables: output/tables/\n")
cat("  • Models: output/models/\n")

cat("\n========== KEY FINDINGS ==========\n")
cat("1. Stationarity:", ifelse(stationarity_result$d > 0,
                               paste("Achieved after", stationarity_result$d, "difference(s)"),
                               "Series is stationary"), "\n")
cat("2. Seasonality:", ifelse(seasonality_present, "Detected", "Not significant"), "\n")
cat("3. Selected Model:", summary_report$Model_Specification, "\n")
cat("4. Model Adequacy:", ifelse(model_adequate, "✓ Adequate", "✗ Needs improvement"), "\n")
cat("5. Test RMSE:", round(summary_report$Test_RMSE, 2), "USD/barrel\n")
cat("6. Test MAPE:", round(summary_report$Test_MAPE, 2), "%\n")

if(!model_adequate) {
  cat("\n⚠ NEXT STEPS:\n")
  cat("  1. Check for outliers and structural breaks\n")
  cat("  2. Consider seasonal differencing if patterns persist\n")
  cat("  3. Try intervention analysis for known events (2008 crisis, COVID)\n")
}

if(diagnostic_results$P_Value[2] < 0.05) {
  cat("\n📝 NOTE ON NON-NORMAL RESIDUALS:\n")
  cat("  This is EXPECTED and ACCEPTABLE for crude oil prices.\n")
  cat("  Literature: Commodity prices typically show:\n")
  cat("  • Fat tails (extreme events)\n")
  cat("  • Volatility clustering\n")
  cat("  • Non-Gaussian distributions\n")
  cat("  ➜ Point forecasts remain valid; use bootstrapped CIs\n")
}

cat("\n✓ Script execution complete!\n")



#ONLY USE IF VARIANCE ISN'T STABLE OR NORMAL

# ============================================================================
# 12. BOX-COX TRANSFORMATION (ALTERNATIVE)
# ============================================================================
cat("\n========== BOX-COX TRANSFORMATION ==========\n")
# Box-Cox finds optimal lambda for transformation
lambda <- BoxCox.lambda(crude_ts)
cat("\nOptimal Box-Cox lambda:", lambda, "\n")

cat("\nInterpretation:\n")
cat("λ = 1   → No transformation needed\n")
cat("λ = 0   → Log transformation (what we used)\n")
cat("λ = 0.5 → Square root transformation\n")
cat("λ = -1  → Inverse transformation\n")

if(abs(lambda) < 0.1) {
  cat("\n✓ Box-Cox suggests log transformation (λ ≈ 0)\n")
} else {
  cat("\n⚠ Box-Cox suggests different transformation (λ =", lambda, ")\n")
}

# Apply Box-Cox transformation
crude_ts_boxcox <- BoxCox(crude_ts, lambda = lambda)

# Plot comparison
png("output/figures/08_original_vs_boxcox.png", width = 1400, height = 800)
par(mfrow = c(2, 1))
plot(crude_ts, main = "Original", ylab = "Price", col = "darkblue", lwd = 2)
plot(crude_ts_boxcox, main = paste("Box-Cox Transform (λ =", round(lambda, 2), ")"), 
     ylab = "Transformed", col = "darkgreen", lwd = 2)
par(mfrow = c(1, 1))
dev.off()
cat("\n✓ original vs boxcox\n")


cat("\n========== SECTION 4: STATIONARITY TESTING ==========\n")
# Function for comprehensive stationarity tests
test_stationarity <- function(ts_data, name, alpha = 0.05) {
  cat("\n--- Testing:", name, "---\n")
  
  # 1. ADF Test (H0: Unit root/non-stationary)
  adf <- adf.test(ts_data)
  
  # 2. PP Test (Phillips-Perron - robust to heteroskedasticity)
  pp <- PP.test(ts_data)
  
  # 3. KPSS Test (H0: Stationary)
  kpss <- kpss.test(ts_data)
  
  results <- data.frame(
    Test = c("ADF", "PP", "KPSS"),
    Statistic = c(adf$statistic, pp$statistic, kpss$statistic),
    P_Value = c(adf$p.value, pp$p.value, kpss$p.value),
    Conclusion = c(
      ifelse(adf$p.value < alpha, "Stationary", "Non-stationary"),
      ifelse(pp$p.value < alpha, "Stationary", "Non-stationary"),
      ifelse(kpss$p.value > alpha, "Stationary", "Non-stationary")
    )
  )
  
  print(results)
  
  # Consensus decision
  stationary_count <- sum(results$Conclusion == "Stationary")
  
  if(stationary_count >= 2) {
    cat("\n✓ CONSENSUS: Series is STATIONARY (", stationary_count, "/3 tests)\n")
    return(list(stationary = TRUE, data = ts_data, d = 0))
  } else {
    cat("\n✗ CONSENSUS: Series is NON-STATIONARY (", 3-stationary_count, "/3 tests)\n")
    cat("   → Applying first differencing...\n")
    
    ts_diff <- diff(ts_data)
    
    # Re-test differenced series
    adf_diff <- adf.test(ts_diff)
    pp_diff <- PP.test(ts_diff)
    kpss_diff <- kpss.test(ts_diff)
    
    results_diff <- data.frame(
      Test = c("ADF", "PP", "KPSS"),
      Statistic = c(adf_diff$statistic, pp_diff$statistic, kpss_diff$statistic),
      P_Value = c(adf_diff$p.value, pp_diff$p.value, kpss_diff$p.value),
      Conclusion = c(
        ifelse(adf_diff$p.value < alpha, "Stationary", "Non-stationary"),
        ifelse(pp_diff$p.value < alpha, "Stationary", "Non-stationary"),
        ifelse(kpss_diff$p.value > alpha, "Stationary", "Non-stationary")
      )
    )
    
    cat("\n   Differenced boxcox Series Results:\n")
    print(results_diff)
    
    stationary_count_diff <- sum(results_diff$Conclusion == "Stationary")
    
    if(stationary_count_diff >= 2) {
      cat("\n✓ Differenced series is STATIONARY (", stationary_count_diff, "/3 tests)\n")
      return(list(stationary = TRUE, data = ts_diff, d = 1, 
                  original_results = results, diff_results = results_diff))
    } else {
      cat("\n⚠ Warning: Series may require d=2 or have structural breaks\n")
      return(list(stationary = FALSE, data = ts_diff, d = 1,
                  original_results = results, diff_results = results_diff))
    }
  }
}

# Test boxcox series
stationarity_result <- test_stationarity(crude_ts_boxcox, "Boxcox Series")

# Save results
stationarity_df <- if(stationarity_result$d == 0) {
  stationarity_result$original_results
} else {
  rbind(
    cbind(Series = "boxcox Original", stationarity_result$original_results),
    cbind(Series = "boxcox Differenced", stationarity_result$diff_results)
  )
}

write.csv(stationarity_df, "output/tables/boxcoxstationarity_tests.csv", row.names = FALSE)
cat("\n✓ boxcox Stationarity test results saved\n")

#MODEL SELECTION
auto_model_box <- auto.arima(crude_ts_boxcox,
                             seasonal = FALSE,
                             stationary = FALSE,
                             stepwise = FALSE,
                             approximation = FALSE,
                             trace = TRUE)



cat("\n========== METHOD 1: BACK-TRANSFORMED PLOT (RECOMMENDED) ==========\n")

# Step 1: Get fitted values in Box-Cox scale
fitted_boxcox <- fitted(auto_model_box)

# Step 2: Back-transform to original scale
fitted_original <- InvBoxCox(fitted_boxcox, lambda)

# EXPLANATION:
# - InvBoxCox() reverses the Box-Cox transformation
# - If Y = BoxCox(X, λ), then X = InvBoxCox(Y, λ)
# - This gives fitted values in original USD scale

# Step 3: Create data frame for plotting
data_plot_backtransformed <- data.frame(
  Time = time(crude_ts),
  Actual = as.numeric(crude_ts),           # Original scale (USD)
  Fitted = as.numeric(fitted_original)     # Back-transformed to USD
)

# Step 4: Plot (both on same scale)
png("output/figures/09_actual_vs_fitted_boxcox_BACKTRANSFORMED.png", 
    width = 1400, height = 800)

ggplot(data_plot_backtransformed, aes(x = Time)) +
  geom_line(aes(y = Actual, color = "Actual"), linewidth = 1) +
  geom_line(aes(y = Fitted, color = "Fitted"), linewidth = 1) +
  scale_color_manual(values = c("Actual" = "black", "Fitted" = "darkgreen")) +
  labs(
    title = "Box-Cox Model: Actual vs Fitted (Back-Transformed to Original Scale)",
    subtitle = paste0("λ = ", round(lambda, 3)),
    y = "Price (USD/barrel)", 
    x = "Time", 
    color = "Series"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "top"
  )

dev.off()

cat("✓ Back-transformed plot saved\n")
cat("  Location: output/figures/09_actual_vs_fitted_boxcox_BACKTRANSFORMED.png\n\n")

# Print sample values to show transformation
cat("--- Sample Comparison ---\n")
cat("Time        | Actual (USD) | Fitted BoxCox | Fitted USD (Back-transformed)\n")
cat("------------|--------------|---------------|---------------------------------\n")
for(i in c(1, 50, 100, 150, 200, length(crude_ts))) {
  cat(sprintf("%s | $%-11.2f | %-13.3f | $%-11.2f\n",
              format(time(crude_ts)[i], nsmall = 2),
              crude_ts[i],
              fitted_boxcox[i],
              fitted_original[i]))
}

cat("\nNOTICE:\n")
cat("  • Actual and Back-transformed Fitted are in SAME UNITS (USD)\n")
cat("  • Box-Cox fitted values are in DIFFERENT UNITS (transformed space)\n")
cat("  • This is why we MUST back-transform for comparison!\n\n")


# ============================================================================
# METHOD 2: WITHOUT BACK-TRANSFORMATION (NOT RECOMMENDED)
# ============================================================================

cat("\n========== METHOD 2: WITHOUT BACK-TRANSFORMATION (NOT RECOMMENDED) ==========\n")

# Plot in Box-Cox space
data_plot_transformed <- data.frame(
  Time = time(crude_ts_boxcox),
  Actual_BoxCox = as.numeric(crude_ts_boxcox),
  Fitted_BoxCox = as.numeric(fitted_boxcox)
)

png("output/figures/10_actual_vs_fitted_boxcox_NOT_BACKTRANSFORMED.png", 
    width = 1400, height = 800)

ggplot(data_plot_transformed, aes(x = Time)) +
  geom_line(aes(y = Actual_BoxCox, color = "Actual (Box-Cox)"), linewidth = 1) +
  geom_line(aes(y = Fitted_BoxCox, color = "Fitted (Box-Cox)"), linewidth = 1) +
  scale_color_manual(values = c("Actual (Box-Cox)" = "black", 
                                "Fitted (Box-Cox)" = "orange")) +
  labs(
    title = "Box-Cox Model: Actual vs Fitted (NOT Back-Transformed)",
    subtitle = "⚠️ WARNING: Both series in transformed space (hard to interpret!)",
    y = "Box-Cox Transformed Values", 
    x = "Time", 
    color = "Series"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "red"),
    legend.position = "top"
  )

dev.off()

cat("✓ Non-back-transformed plot saved (for comparison only)\n")
cat("  Location: output/figures/10_actual_vs_fitted_boxcox_NOT_BACKTRANSFORMED.png\n\n")

cat("⚠️  PROBLEM WITH THIS PLOT:\n")
cat("  • Y-axis values are meaningless (Box-Cox units)\n")
cat("  • Can't tell if price predictions are good in USD terms\n")
cat("  • Hard to communicate results to stakeholders\n")
cat("  • Use ONLY for technical diagnostics, not reporting!\n\n")



# ============================================================================
# CALCULATE FIT QUALITY METRICS (ON ORIGINAL SCALE)
# ============================================================================

cat("\n========== FIT QUALITY ASSESSMENT ==========\n")

# Calculate errors (must be on same scale!)
errors <- crude_ts - fitted_original  # Both in USD

# Metrics
mae <- mean(abs(errors))
rmse <- sqrt(mean(errors^2))
mape <- mean(abs(errors / crude_ts)) * 100
mse <- mean(errors^2)

cat("\nBox-Cox Model Performance (Original Scale):\n")
cat("  MAE:  $", round(mae, 2), "\n", sep = "")
cat("  RMSE: $", round(rmse, 2), "\n", sep = "")
cat("  MAPE: ", round(mape, 2), "%\n", sep = "")
cat("  MSE:  ", round(mse, 2), "\n\n", sep = "")

cat("INTERPRETATION:\n")
cat("  • MAPE = ", round(mape, 2), "% means fitted values are off by ~", 
    round(mape, 1), "% on average\n", sep = "")
cat("  • RMSE = $", round(rmse, 2), " means typical error is about $", 
    round(rmse, 0), "\n", sep = "")

if(mape < 10) {
  cat("  ✓ EXCELLENT in-sample fit!\n")
} else if(mape < 15) {
  cat("  ✓ GOOD in-sample fit\n")
} else if(mape < 25) {
  cat("  ⚠ MODERATE fit - consider improvements\n")
} else {
  cat("  ✗ POOR fit - model may not be capturing patterns well\n")
}


# ============================================================================
# RESIDUAL PLOT (ALSO BACK-TRANSFORMED)
# ============================================================================

cat("\n========== RESIDUAL ANALYSIS ==========\n")

png("output/figures/12_boxcox_residuals_backtransformed.png", 
    width = 1400, height = 1000)

par(mfrow = c(2, 2))

# Residuals over time (in USD)
plot(time(crude_ts), errors, type = "l", 
     main = "Residuals Over Time (Original Scale)",
     ylab = "Residuals (USD)", xlab = "Time",
     col = "darkblue", lwd = 1.5)
abline(h = 0, col = "red", lty = 2, lwd = 2)
grid()

# Residuals vs Fitted
plot(fitted_original, errors,
     main = "Residuals vs Fitted Values",
     xlab = "Fitted Values (USD)", ylab = "Residuals (USD)",
     pch = 20, col = rgb(0, 0, 0, 0.5))
abline(h = 0, col = "red", lty = 2, lwd = 2)
grid()

# Histogram of residuals
hist(errors, breaks = 30, col = "steelblue",
     main = "Distribution of Residuals",
     xlab = "Residuals (USD)", probability = TRUE)
lines(density(errors), col = "red", lwd = 2)

# Q-Q plot
qqnorm(errors, main = "Normal Q-Q Plot of Residuals")
qqline(errors, col = "red", lwd = 2)

par(mfrow = c(1, 1))
dev.off()

cat("✓ Residual diagnostics saved\n\n")

str(fitted_original)


# ============================================================================
# COMPARISON WITH ORIGINAL MODEL (IF AVAILABLE)
# ============================================================================

if(exists("auto_model")) {
  cat("\n========== COMPARING BOX-COX VS ORIGINAL MODEL ==========\n")
  
  # Original model fitted values
  fitted_original_model <- fitted(auto_model)
  errors_original_model <- crude_ts - fitted_original_model
  
  # Metrics
  mae_orig <- mean(abs(errors_original_model))
  rmse_orig <- sqrt(mean(errors_original_model^2))
  mape_orig <- mean(abs(errors_original_model / crude_ts)) * 100
  
  comparison_table <- data.frame(
    Model = c("Original Scale", "Box-Cox Transformed"),
    MAE = c(mae_orig, mae),
    RMSE = c(rmse_orig, rmse),
    MAPE = c(mape_orig, mape)
  )
  
  cat("\n--- Model Comparison (In-Sample Fit) ---\n")
  print(comparison_table)
  
  # Determine winner
  cat("\nBest Model by Metric:\n")
  cat("  MAE:  ", comparison_table$Model[which.min(comparison_table$MAE)], "\n")
  cat("  RMSE: ", comparison_table$Model[which.min(comparison_table$RMSE)], "\n")
  cat("  MAPE: ", comparison_table$Model[which.min(comparison_table$MAPE)], "\n")
  
  # Save the comparison table
  write.csv(comparison_table, "output/tables/model_comparison_boxcox_vs_original.csv", row.names = FALSE)
  cat("\n✓ comparison table results saved\n")
  
  # Close any lingering devices
  graphics.off()
  
  # Reopen a clean one
  png("output/figures/13_original_vs_boxcox_fitted.png", width = 1400, height = 800)
  
  plot(as.numeric(time(crude_ts)), as.numeric(crude_ts), type = "l",
       main = "Comparing Original & Box-Cox Model Fits with Actual",
       ylab = "Price (USD/barrel)", xlab = "Time",
       col = "black", lwd = 2, ylim = range(crude_ts))
  
  lines(as.numeric(time(crude_ts)), as.numeric(fitted_original_model),
        col = "blue", lwd = 1.5, lty = 2)
  
  lines(as.numeric(time(crude_ts)), as.numeric(fitted_original),
        col = "red", lwd = 1.5, lty = 3)
  
  legend("topleft",
         c("Actual", "Original Model", "Box-Cox Model"),
         col = c("black", "blue", "red"),
         lty = c(1, 2, 3),
         lwd = c(2, 1.5, 1.5))
  grid()
  
  dev.off()
  cat("\n✓ Comparison plot saved successfully\n")
}


# ============================================================================
# KEY TAKEAWAYS & RECOMMENDATIONS
# ============================================================================


cat("1. ALWAYS BACK-TRANSFORM FOR REPORTING\n")
cat("   ✓ Use InvBoxCox(fitted_values, lambda)\n")
cat("   ✓ This converts Box-Cox scale → Original scale (USD)\n")
cat("   ✓ Allows direct comparison with actual prices\n\n")

cat("2. WHY BACK-TRANSFORM?\n")
cat("   • Interpretability: Stakeholders understand USD, not Box-Cox units\n")
cat("   • Fair comparison: Can't compare different scales\n")
cat("   • Accurate metrics: MAE, RMSE, MAPE must be on same scale\n\n")

cat("3. WHEN TO USE EACH PLOT:\n")
cat("   📊 Back-transformed plot:\n")
cat("      → For reports, presentations, decision-making\n")
cat("      → To assess practical forecast accuracy\n")
cat("      → To communicate with non-technical audience\n\n")
cat("   📊 Non-back-transformed plot:\n")
cat("      → For technical diagnostics only\n")
cat("      → To check model fit in transformed space\n")
cat("      → Internal analysis purposes\n\n")

cat("4. FORECASTING:\n")
cat("   When you generate forecasts:\n")
cat("   forecast_boxcox <- forecast(auto_model_box, h = 12)\n")
cat("   forecast_original <- InvBoxCox(forecast_boxcox$mean, lambda)\n")
cat("   ↑ ALWAYS back-transform before reporting!\n\n")

cat("5. BOX-COX LAMBDA INTERPRETATION:\n")
cat("   Your λ = ", round(lambda, 3), "\n", sep = "")
if(abs(lambda) < 0.1) {
  cat("   → Close to 0: Similar to log transformation\n")
  cat("   → Stabilizes variance when it increases with level\n")
} else if(abs(lambda - 0.5) < 0.1) {
  cat("   → Close to 0.5: Square root transformation\n")
  cat("   → Moderately stabilizes variance\n")
} else if(abs(lambda - 1) < 0.1) {
  cat("   → Close to 1: No transformation needed\n")
  cat("   → Original scale is already appropriate\n")
} else {
  cat("   → Custom transformation for optimal variance stabilization\n")
}



cat("\n # ============================================================================ \n")
cat(" # SAVE RESULTS TO TABLE \n")


# Create results summary
boxcox_results <- data.frame(
  Lambda = lambda,
  MAE = mae,
  RMSE = rmse,
  MAPE = mape,
  Model_Specification = paste0(auto_model_box$arma[1], ",", 
                               auto_model_box$arma[6], ",", 
                               auto_model_box$arma[2])
)

write.csv(boxcox_results, "output/tables/boxcox_model_results.csv", row.names = FALSE)
cat("✓ Results saved to: output/tables/boxcox_model_results.csv\n")




cat("\n==========  RESIDUAL DIAGNOSTICS (BOX-COX) ==========\n")

# Extract residuals from Box-Cox model
residuals_box <- residuals(auto_model_box)

# Initialize diagnostics table
diagnostic_results_box <- data.frame(
  Test = character(),
  Statistic = numeric(),
  P_Value = numeric(),
  Conclusion = character(),
  stringsAsFactors = FALSE
)

# 1. Ljung-Box Test (Autocorrelation)
lb_test_box <- Box.test(residuals_box, lag = 20, type = "Ljung-Box", 
                        fitdf = sum(auto_model_box$arma[1:3]))
diagnostic_results_box <- rbind(diagnostic_results_box, data.frame(
  Test = "Ljung-Box (Autocorrelation)",
  Statistic = as.numeric(lb_test_box$statistic),
  P_Value = lb_test_box$p.value,
  Conclusion = ifelse(lb_test_box$p.value > 0.05, "✓ No autocorrelation", "✗ Autocorrelation present")
))

# 2. Shapiro-Wilk Test (Normality)
sw_test_box <- shapiro.test(residuals_box)
diagnostic_results_box <- rbind(diagnostic_results_box, data.frame(
  Test = "Shapiro-Wilk (Normality)",
  Statistic = as.numeric(sw_test_box$statistic),
  P_Value = sw_test_box$p.value,
  Conclusion = ifelse(sw_test_box$p.value > 0.05, "✓ Normal residuals", "✗ Non-normal residuals")
))

# 3. ARCH Test (Heteroskedasticity)
arch_test_box <- ArchTest(residuals_box, lags = 20)
diagnostic_results_box <- rbind(diagnostic_results_box, data.frame(
  Test = "ARCH (Heteroskedasticity)",
  Statistic = as.numeric(arch_test_box$statistic),
  P_Value = arch_test_box$p.value,
  Conclusion = ifelse(arch_test_box$p.value > 0.05, "✓ No ARCH effects", "✗ ARCH effects present")
))

# 4. Jarque-Bera Test (Alternative normality test)
jb_test_box <- jarque.bera.test(residuals_box)
diagnostic_results_box <- rbind(diagnostic_results_box, data.frame(
  Test = "Jarque-Bera (Normality)",
  Statistic = as.numeric(jb_test_box$statistic),
  P_Value = jb_test_box$p.value,
  Conclusion = ifelse(jb_test_box$p.value > 0.05, "✓ Normal residuals", "✗ Non-normal residuals")
))

# Print results
cat("\n--- Box-Cox Model Residual Diagnostic Results ---\n")
print(diagnostic_results_box)

# Save to CSV
write.csv(diagnostic_results_box, "output/tables/residual_diagnostics_boxcox.csv", row.names = FALSE)
cat("\n✓ Box-Cox residual diagnostics saved: output/tables/residual_diagnostics_boxcox.csv\n")

cat("\n BOXCOX TRANSFORMATION DIDN'T NECESSARILY IMPROVE OF MODEL OR PREVENT THE LAG PREDICTION, MAYBE I SHOULD TRY GARCH?? \n")



#EXPORTING ALL WORK AND OUTPUT (LATEX)
# knitr::stitch('BESTARIMA.r')
# browseURL('BESTARIMA.pdf')
