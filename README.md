bayesFloodStage

README file last updated by Christopher Piecuch, cpiecuch-at-whoi-dot-edu, Thu Jul 17 2025

Basic description

Citation

This code was generated to produce Bayesian modeling results for LIA' VLM presented in the main text of:

Lewright, L., Austerman, J., Piecuch, C. G., Adhikari, S., Davis, J. L., Milne, G. A., Paxman, G. J. G., “Projections of 21st Century Sea-Level Fall Along Coastal Greenland”, Nature Communications, xxx.

Please cite this reference when using this code.

Contents

Text files
* Berg_Uplift.txt: data file from Berg et al. (2024) https://doi.org/10.1029/2023GL104851
* Copyright: copyright statement
* License: license statement

MATLAB .mat files
* RF_correction_BERG.mat: Reference-frame correction at GNSS locations
* VLM_st_Berg.mat: 14,040 model predictions of LIA' VLM at GNSS locations
* deglacialRrate.mat: Deglacial GIA VLM estimates with uncertainties at GNSS locations

Main MATLAB .m files
* bayes_main_code.m: Bayesian model code. Execute to produce posterior solutions.

Supporting MATLAB .m files
* EarthDistances.m: Compute great-circle distance between latitudes and longitudes (called in bayes_main_code)
* delete_burn_in.m: Delete burn in (called in bayes_main_code)
* drchrnd.m: Draw dirichlet random variable (called in bayes_main_code)
* initialize_output.m: Initialize posterior model output (called in bayes_main_code)
* prepare_gps_data.m: Load and format GPS VLM rate data (called in bayes_main_code)
* set_hyperparameters.m: Set hyperparameters in Bayesian model code (called in bayes_main_code)
* set_initial_values.m: Set initial conditions in Bayesian model code (called in bayes_main_code)
* update_all_arrays.m: Update model output for each iteration of the Gibbs/Metropolis sampler in Bayesian model code (called in bayes_main_code)

Note on usage of code and MATLAB release This code was written and executed (to produce Bayesian modeling results for LIA' VLM in Lewright et al. 2025) using MATLAB release version MATLAB_R2024a.
