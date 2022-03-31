#!/bin/bash
#SBATCH -o ./output/sim_study_meta/logs/results_%a.txt                                                              
#SBATCH -e ./output/sim_study_meta/logs/errors_%a.err                                                                     
#SBATCH -p dunsonlab --account=dunsonlab                                                                                
#SBATCH --mem=10G                                                                                                       
#SBATCH -a 1-10       

module load R
Rscript ./scripts/analysis_scripts/sim_study_meta_master_script.R
module unload R