#!/bin/bash
#SBATCH -o ./output/us_analysis/logs/results_%a.txt                                                              
#SBATCH -e ./output/us_analysis/logs/errors_%a.err                                                                     
#SBATCH -p dunsonlab --account=dunsonlab                                                                                
#SBATCH --mem=10G
#SBATCH -a 1-81

module load R
Rscript ./scripts/analysis_scripts/us_analysis_master_script.R
module unload R