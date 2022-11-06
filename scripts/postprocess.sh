#!/bin/bash
#SBATCH -o ./output/us_analysis/logs/postprocess_results.txt                                                              
#SBATCH -e ./output/us_analysis/logs/postprocess_errors.err                                                                     
#SBATCH -p dunsonlab --account=dunsonlab                                                                                
#SBATCH --mem=10G

module load R
Rscript ./scripts/analysis_scripts/sim_study_meta_postprocess.R
Rscript ./scripts/analysis_scripts/rolling_qdl_comparison.R
Rscript ./scripts/analysis_scripts/us_analysis_nonparametric.R
Rscript ./scripts/analysis_scripts/us_analysis_postprocessing_endpoints.R
Rscript ./scripts/analysis_scripts/us_analysis_postprocessing_sensitivity.R
module unload R