//ITERATION THIS IS FOR
global nsim = 2000
global type = "phase1_ss"
global true = 0

global conf_type "1_smallcon" "1_noconf" "10_smallco" "10_noconf"


// PUT ALL ESTIMATE FILES TOGETHER AND SAVE TO THE MAIN RESULTS
global doDir "H:\mtchdcohortsim"
global resultsDir "H:\mtchdcohortsim\Main_RESULTS"
global rawDir "H:\mtchdcohortsim\hpc\Output_2000_12NOV2024"


//ONLY RUN IF YOU NEED THE RAW ACCESS FILES//

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

// RUN GRAPHS

//Estimates

//Graphs
cd $doDir
do sim_graph_simsum_phase1_ss
do sim_graph_empSE_ratio_phase1_ss
do sim_graph_rel_err_ratio_phase1_ss
do sim_graph_bias_model_phase1_ss
do sim_graph_cov_model_phase1_ss

graph combine empSE_all_smallconf rel_err_all_smallconf , row(2) name(precision_all_smallconf, replace)
graph export "$resultsDir\sim_graph_precision_ratio_all_smallconf_${nsim}_${type}.emf", replace

graph combine empSE_all_noconf rel_err_all_noconf , row(2) name(precision_all_noconf, replace)
graph export "$resultsDir\sim_graph_precision_ratio_all_noconf_${nsim}_${type}.emf", replace

use "$resultsDir\sim_summary_simsumbc1_${nsim}_${type}", clear
list if test == "Bias" & confounding == "10_smallco"
list if test == "Bias" & confounding == "1_smallcon"


use "$resultsDir\sim_estimates_${nsim}_${type}", clear
drop if ratio == 0
collapse (mean) nobs (median) nomatch, by(confounding getallpossible ratio)
list


use "$resultsDir\sim_summary_simsummtch1_${nsim}_${type}", clear


graph dir
local graphs `r(list)'
foreach g of local graphs {
    if strpos("`g'", "_all_") == 0 {
        graph drop `g'
    }

}
