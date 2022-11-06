library(DescTools)
library(latex2exp)


param_options <- expand.grid(ifr_level = c("L", "M", "H"), # index over using lower, middle, uppper end of range
                             infectious_period = c(1,1.5,2.0),
                             lengthscale = c(8,12,16)
)

NV <- nrow(param_options)

## --- LOAD MODEL FITS: EFFECTS, FORECASTS from 74 WEEKS, FORECASTS from 77 WEEKS ---
effect_summary <- matrix(nrow = NV, ncol = 3, dimnames = list(NULL, c("Mean", "Low", "High")))
for(version in 1:NV){
  sfit <- readRDS(paste0("output/us_analysis/sfit_90_", version, ".rds"))
  effect_summary[version,] <- sfit$summary$Global[2,c(1,3,4)]

}

median_forecast_74 <- matrix(nrow = 90, ncol = NV)
for(version in 1:NV){
  sfit <- readRDS(paste0("output/us_analysis/sfit_74_", version, ".rds"))
  sfor <- smesir_forecast(90-74, sfit,
                          new_x = list(delta_sq = delta_prevalence[75:90,]^2), 
                          new_vaccinations = covid_vaccinations[75:90,])
  median_forecast_74[,version] <- rowSums(apply(sfor$expected_samples, c(1,3), quantile, 0.6))
}

median_forecast_77 <- matrix(nrow = 90, ncol = NV)
for(version in 1:NV){
  sfit <- readRDS(paste0("output/us_analysis/sfit_77_", version, ".rds"))
  sfor <- smesir_forecast(90-77, sfit,
                          new_x = list(delta_sq = delta_prevalence[78:90,]^2), 
                          new_vaccinations = covid_vaccinations[78:90,])
  median_forecast_77[,version] <- rowSums(apply(sfor$expected_samples, c(1,3), quantile, 0.6))
}


## --- PREPARE PLOT LABELS ---
# row labels for the 27 effect estimates
plot_labels <- c()
for(v in 1:NV){
  plot_labels[v] <- TeX(paste0(
    sprintf(r'(IFR $= %s$, $\gamma^{-1} = %1.1f$, )', param_options[v,1], param_options[v,2]), 
    '\u2113', 
    sprintf(r'($= %1.0f$)',param_options[v,3])
  ))
}

# monthly date vector for the independent variable in the forecasting plot
weeks <- ymd(rownames(covid_deaths))
yii <- year(min(weeks))
mii <- month(min(weeks))
monthly <- c()
while(yii < year(max(weeks)) || mii <= month(max(weeks))){
  monthly <- c(monthly, paste0(yii,"-",mii,"-01"))
  
  mii <- mii + 1
  if(mii > 12){
    yii <- yii + 1
    mii <- 1
  }
}
monthly <- ymd(monthly)


## --- CREATE THE PLOTS ---
png("output/us_analysis/sensitivity/effect_sensitivity.png", width = 6, height = 12, units = 'in', res = 300)
PlotDot(effect_summary[,1],
        pch = rep(c(1,2,4), times = 9),
        args.errbars = list(from = effect_summary[,2], to = effect_summary[,3], col = rep(1:3, each = 3)),
        labels = plot_labels,
        xlim = c(0,2), groups = rep(1:3, each = 9), color = rep(1:3, each = 3),
        main = "Posterior Distribution of\nDelta Prevalence Coefficient"
        )
dev.off()

PLOTSTART <- 50 # lower bound on time
png("output/us_analysis/sensitivity/forecast_sensitivity.png", width = 8, height = 8, units = 'in', res = 300)
plot(weeks[PLOTSTART:90], rowSums(covid_deaths[PLOTSTART:90,]), xaxt = 'n', xlab = NA, ylab = NA, main = "US COVID-19 Deaths\nForecast Comparison", type = "l",lwd=2)
axis(1, monthly,format(monthly, "%b %y"), cex.axis = 1)
for(version in 1:NV){
  lines(weeks[75:90], median_forecast_74[75:90,version], 
        col = floor((version - 1)/ 9) + 1)
  lines(weeks[78:90], median_forecast_77[78:90,version], 
        col = floor((version - 1)/ 9) + 1, lty = 3)
}
legend_labels <- c(TeX('\u2113 $= 8$'), TeX('\u2113 $= 12$'), TeX('\u2113 $= 16$'))
legend("topleft",legend = legend_labels, col = 1:3, lty = 1)
dev.off()

