library(DescTools)
library(latex2exp)

# load previously fitted models
param_options <- expand.grid(ifr_level = c("L","M","H"), 
                             infectious_period = c("L","M","H"),
                             lengthscale = c("L","M","H")
)
NV <- nrow(param_options)
effect_summary <- matrix(nrow = NV, ncol = 3, dimnames = list(NULL, c("Mean", "Low", "High")))
for(version in 1:NV){
  sfit <- readRDS(paste0("output/us_analysis/sfit_90_", version, ".rds"))
  effect_summary[version,] <- sfit$summary$Global[2,c(1,3,4)]
  print(sfit$mcmc_diagnostics[,1])
}


## It looks like not all of these mixed well - we should check and then
## re-run the sensitivity analysis with the fixed large delta!!!
PlotDot(effect_summary[,1], 
        args.errbars = list(from = effect_summary[,2], to = effect_summary[,3]),
        labels = paste(param_options[,1], param_options[,2], param_options[,3]),
        xlim = c(0,2)
        )
