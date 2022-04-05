#!/bin/bash
#SBATCH -o ./output/sim_study_case/logs/results.txt                                                              
#SBATCH -e ./output/sim_study_case/logs/errors.err                                                                     
#SBATCH -p dunsonlab --account=dunsonlab                                                                                
#SBATCH --mem=10G                                                                                                       

module load R
Rscript ./scripts/analysis_scripts/sim_study_case.R
module unload R
