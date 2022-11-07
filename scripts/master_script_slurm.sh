#!/bin/bash

sbatch scripts/analysis_scripts/sim_study_case_slurm.sh
sbatch scripts/analysis_scripts/sim_study_meta_slurm.sh
sbatch scripts/analysis_scripts/us_analysis_slurm.sh
