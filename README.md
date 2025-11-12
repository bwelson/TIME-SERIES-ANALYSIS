# Time Series Analysis of Brent Crude Oil Prices

[![R](https://img.shields.io/badge/R-4.4.2-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

Comprehensive time series analysis comparing ARIMA and Box-Cox transformed models for Brent crude oil price forecasting (1999-2023). Implements rigorous stationarity testing, residual diagnostics, and out-of-sample validation.

## Key Findings

- **Model**: ARIMA(1,1,0) optimal per AIC criteria
- **In-sample fit**: MAPE 6.71%, RMSE $5.29/barrel
- **Diagnostics**: Model adequate (no autocorrelation), but exhibits ARCH effects
- **Box-Cox transformation**: Minimal improvement (λ = -0.672)

## Repository Structure
```
├── data/
│   └── crude.csv                    # Monthly Brent crude prices
├── scripts/
│   └── BESTARIMA.R                  # Complete analysis pipeline
├── output/
│   ├── figures/                     # All diagnostic plots
│   ├── tables/                      # Results CSVs
│   └── models/                      # Saved ARIMA objects
├── manuscript/
│   └── manuscript.tex               # Publication draft
└── README.md
```

## Methodology

1. **Stationarity Testing**: ADF, PP, KPSS (consensus approach)
2. **Model Selection**: Exhaustive auto.arima search
3. **Diagnostics**: Ljung-Box, Shapiro-Wilk, ARCH, Jarque-Bera tests
4. **Validation**: 80-20 train-test split, RMSE/MAE/MAPE metrics
5. **Transformation**: Box-Cox with proper back-transformation

## Results Summary

| Model | MAE | RMSE | MAPE | ARCH Effects |
|-------|-----|------|------|--------------|
| ARIMA(1,1,0) | 4.03 | 5.29 | 6.71% | Yes (p=0.001) |
| Box-Cox ARIMA(0,1,4) | 4.06 | 5.42 | 6.70% | No (p=0.121) |

## Usage
```r
# Install required packages
install.packages(c("forecast", "tseries", "lmtest", "FinTS", "tidyverse"))

# Run analysis
source("scripts/BESTARIMA.R")

# Load saved model
model <- readRDS("output/models/final_arima_model.rds")
```

## Citation

If you use this code, please cite:
```
[Your Name]. (2025). Time Series Modeling of Brent Crude Oil Prices: 
A Comparative Analysis of ARIMA and Box-Cox Transformed Models. 
GitHub repository: https://github.com/[username]/crude-oil-arima-analysis
```

## License

MIT License - see LICENSE file for details

## Contact

- **Author**: [BENTUM WELSON]
- **Email**: bwelson523.email@gmail.com
- **LinkedIn**: [www.linkedin.com/in/bentumwelson523/]
