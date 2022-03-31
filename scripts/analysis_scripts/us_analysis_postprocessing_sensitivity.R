library(DescTools)
library(latex2exp)

# sensitivity analysis, scanning over different versions of the model fit
endpoint <- 90 # so far we've only focused on the 90 week conclusions for our sensitivity analysis
NV <- 27 # we conducted the analysis with 27 different epi parameter specifications

# load previously fitted models
sfit <- lapply(1:NV, function(v){
  readRDS(paste0("output/us_analysis/sfit_",endpoint, "_", v, ".rds"))
})

# Currently, the below just checks:
# - sensitivity of the MCMC algorithm (NOT STORED)
# - sensitivity of IIP inferences for NY and VT
# - sensitivity of remaining vulnerable for all states (UGLY)

# In the future, we should add
# Sensitivity of our 90 week conclusions:
# - coefficient of delta squared
# - long term deaths upper bound
global_delta_effect <- matrix(nrow = NV, ncol = 3) 
long_run_total_deaths <- matrix(nrow = NV, ncol = 3)

# Sensitivity of our 77 week foreast conclusions:
# - 0.85 quantile forecasts
forecast_examples <- matrix(nrow = NV, ncol = (90 - 77))

# Section 0 - Worst-case Effective Sample Size for each simulation
for(v in 1:NV){
  print(min(unlist(lapply(sfit[[v]]$mcmc_diagnostics,function(x) min(x[,2])))))
}

## Section 1 - Initial Infectious Population (NY and VT)
# Compute posterior CIs
new_york_iip <- matrix(nrow = NV, ncol = 3)
vermont_iip <- matrix(nrow = NV, ncol = 3)
for(v in 1:NV){
  new_york_iip[v,] <- sfit[[v]]$summary[["new york"]][3,c(3,1,4)]
  vermont_iip[v,] <- sfit[[v]]$summary[["vermont"]][3,c(3,1,4)]
}

# Compute prior CIs
iip_prior <- rexp(1e5,rate = 1/10000)
iip_prior_summary <- c(quantile(iip_prior, 0.025), 
                       mean(iip_prior), 
                       quantile(iip_prior, 0.975))

new_york_iip <- rbind(iip_prior_summary, new_york_iip)
vermont_iip <- rbind(iip_prior_summary, vermont_iip)
rownames(new_york_iip) <- NULL
rownames(vermont_iip) <- NULL

# Create plot labels
plot_labels <- c()
for(v in 1:NV){
  plot_labels[v] <- TeX(paste0(
    sprintf(r'(IFR $= %1.0f$, $\gamma^{-1} = %1.1f$, )', param_options[v,1],param_options[v,2]), 
    '\u2113', 
    sprintf(r'($= %1.0f$)',param_options[v,3])
  ))
}
plot_labels <- c("Prior Distribution", plot_labels)

# Create Plot (NY)
quartz(type = "pdf", 
       file = "figures/us_analysis/sensitivity_iip_ny.pdf",
       height = 9, width = 5)
PlotDot(new_york_iip[,2],
        labels = plot_labels,
        #groups = c(1,rep(2,27)),
        args.errbars = list(from = new_york_iip[,1], to = new_york_iip[,3]),
        xlab = "Initial Infectious Population",
        main = "Sensitivity Analysis\nNew York")
dev.off()

# Create Plot (VT)
quartz(type = "pdf", 
       file = "figures/us_analysis/sensitivity_iip_vt.pdf", 
       height = 9, width = 5)
PlotDot(vermont_iip[,2],
        labels = plot_labels,
        args.errbars = list(from = vermont_iip[,1], to = vermont_iip[,3]),
        xlab = "Initial Infectious Population",
        main = "Sensitivity Analysis\nVermont")
dev.off()

# Create plot (both)
quartz(type = "pdf", 
       file = "figures/us_analysis/sensitivity_iip_both.pdf",
       height = 9, width = 9)
par(mfrow = c(1,2))
PlotDot(new_york_iip[,2],
        labels = plot_labels,
        args.errbars = list(from = new_york_iip[,1], to = new_york_iip[,3]),
        xlab = "Initial Infectious Population",
        main = "New York")

PlotDot(vermont_iip[,2],
        args.errbars = list(from = vermont_iip[,1], to = vermont_iip[,3]),
        xlab = "Initial Infectious Population",
        main = "Vermont")
dev.off()




## Section 2 - Remaining Vulnerable (States as rows, replot for each version)
remaining_vulnerable <- list()

for(v in 1:NV){
  region_names <- sfit[[v]]$region_names
  K <- length(region_names)
  
  vulnerable <- matrix(nrow = K, ncol = 3)
  rownames(vulnerable) <- region_names
  for(k in 1:K){
    vulnerable[k,] <- compute_min_vulnerable(sfit[[v]],k)
  }
  remaining_vulnerable[[v]] <- vulnerable
}

quartz(type = "pdf", 
       file = "figures/us_analysis/remaining_vulnerable.pdf", 
       height = 9, width = 5)
# plot for first version
PlotDot(remaining_vulnerable[[1]][,1],
        labels = region_names,
        args.errbars = list(
          from = remaining_vulnerable[[1]][,2], 
          to = remaining_vulnerable[[1]][,3])
        )
# iterate over the remaining versions
lapply(remaining_vulnerable[2:NV],
       function(vnb) {
         PlotDot(vnb[,1], 
                 args.errbars = list(
                   from = vnb[,2], 
                   to = vnb[,3]), 
                 add = TRUE)
       })
dev.off()




