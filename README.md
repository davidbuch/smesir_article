# smesir_article
Code to reproduce the results and figures in the "smesir" article

Reproduce the results from the article "Explaining transmission rate variations and forecasting epidemic spread in multiple regions with a semiparametric mixed effects SIR model", by Buch, Johndrow, and Dunson. The scripts in this repository rely on the ```smesir``` software package, which can be installed from within an R session using the command ```devtools::install_github("davidbuch/smesir")```. In the ```data``` directory, we include some data files from the supplement of "Regression Models for Understanding COVID-19 Epidemic Dynamics With Incomplete Data" by Quick, Dey, and Lin (2021).

To reproduce results from the article, run ```./scripts/master_script.sh``` from the command line. Results will populate the ``output`` directory. *Notice that the script begins with two "for loops" which iterate over a range of hyperparameters considered in a sensitivity analysis and a simulation study. If this is run sequentially, executing the complete set of sensitivity analyses and simulation studies might take anywhere from a few days to two months, based on our rough estimates.*

To reproduce specific results, the `R` scripts can be run separately, for example ```scripts/analysis_scripts/sim_study_case.R``` can be run to generate all results from our simulation case study. 

To run the 81 different versions of the US COVID-19 analysis in parallel, as well as to parallelize the simulation study runs, we have provided SLURM directives, so that ```scripts/master_script_slurm.sh``` can be executed on a computing cluster. This will generate data files in the ```output``` directory which can then be acted upon locally by the R scripts called in the latter portion of ```scripts/master_script.sh```.
