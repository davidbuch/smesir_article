#!/bin/bash

for i in {1 .. 10};
do
Rscript ./scripts/analysis_scripts/sim_study_meta_master_script.R  $i
done 

for i in {1..81};
do
Rscript ./scripts/analysis_scripts/us_analysis_master_script.R  $i
done 

Rscript ./scripts/analysis_scripts/sim_study_case.R
Rscript ./scripts/analysis_scripts/sim_study_meta_postprocess.R
Rscript ./scripts/analysis_scripts/rolling_qdl_comparison.R
Rscript ./scripts/analysis_scripts/us_analysis_nonparametric.R
Rscript ./scripts/analysis_scripts/us_analysis_postprocessing_endpoints.R
Rscript ./scripts/analysis_scripts/us_analysis_postprocessing_sensitivity.R