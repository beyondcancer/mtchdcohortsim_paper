//ITERATION THIS IS FOR
global nsim = 2000
global type = "phase1_ss"
global true = 0


// PUT ALL ESTIMATE FILES TOGETHER AND SAVE TO THE MAIN RESULTS
global workdir /// Directory for simulation repetitions
global doDir /// Dofile directory for preformance measures
global graphDir /// Dofile directory for graphs
global resultsDir /// Directory where to store results
global rawDir ///Directory of simulation results


//ONLY RUN IF YOU NEED THE RAW ACCESS FILES//
///

cd $workdir

do test_hpchr_test.do /// designed to run on hpc

///

clear 

cd $doDir
do sim_combine_estimates.do

//ESTIMATES DATASET:
use "$resultsDir\sim_estimates_${nsim}_${type}", clear

//characteristics of the data
preserve
bysort confounding ratio method getallpossible: gen count_repid = _N
by confounding ratio method getallpossible: keep if _n == 1   // Keep only one row per combination
keep if count_repid != 2000             // Keep only those combinations where count_repid is not equal to 2000
noi disp "there was an error in:"
list confounding ratio method getallpossible count_repid    // Display the results
restore

// CALCULATE SUMMARY ESTIMATES

cd $doDir
do sim_summary_simsum

cd $doDir
do sim_comp_performance

// RUN GRAPHS

//Estimates

//Graphs
cd $graphDir
//do sim_graph_simsum_phase1_ss
do sim_graph_bias_model_phase1_ss
do sim_graph_cov_model_phase1_ss
do sim_graph_empSE_ratio_phase1_ss
do sim_graph_rel_err_ratio_phase1_ss

/*
graph combine empSE_all_sconf rel_err_all_sconf , row(2) name(precision_all_smallconf, replace)
graph export "$resultsDir\sim_graph_precision_ratio_all_smallconf_${nsim}_${type}.emf", replace

graph combine empSE_all_nconf rel_err_all_nconf , row(2) name(precision_all_nconf, replace)
graph export "$resultsDir\sim_graph_precision_ratio_all_noconf_${nsim}_${type}.emf", replace

graph combine empSE_all_nconf rel_err_all_bconf , row(2) name(precision_all_bconf, replace)
graph export "$resultsDir\sim_graph_precision_ratio_all_bconf_${nsim}_${type}.emf", replace
*/

cd $doDir
do sim_graph_barcode_control_phase1_ss.do
do sim_graph_barcode_confounding_phase1_ss.do
do sim_graph_slope_phase1_ss.do

graph combine slope_emp_sconf slope_rel_err_sconf, col(1) name(comb_sconf, replace) imargin(0 0 0 0)
graph combine barcode_control_avail comb_sconf, col(2)  name(dgm_control_avail, replace)
graph export "$resultsDir\dgm_control_avail_${nsim}_${type}.emf", replace

graph combine slope_emp_10 slope_rel_err_10, col(1) name(comb_10, replace) imargin(0 0 0 0)
graph combine barcode_confounding comb_10, col(2) name(dgm_confounding, replace) 
graph export "$resultsDir\dgm_confounding_${nsim}_${type}.emf", replace

cd $doDir
do tool_emf_to_word.do


use "$resultsDir\sim_summary_simsumbc1_${nsim}_${type}", clear
export delimited "$resultsDir\chatacteristics_bc_${nsim}_${type}.csv", replace


use "$resultsDir\sim_estimates_${nsim}_${type}", clear
drop if ratio == 0
collapse (mean) nobs (median) nomatch, by(confounding getallpossible ratio)
list

export delimited "$resultsDir\nomatch_${nsim}_${type}.csv", replace

use "$resultsDir\sim_summary_simsummtch1_${nsim}_${type}", clear




graph dir
local graphs `r(list)'
foreach g of local graphs {
    if strpos("`g'", "_all_") == 0 {
        graph drop `g'
    }

}
