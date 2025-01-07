use "${resultsDir}\rep_cohortdeets_2000_phase1_ss.dta", clear

drop if time_mtch == .
duplicates drop

collapse  (mean) avg_time = time_mtch (median) median_time = time_mtch (max) max_time = time_mtch (count) nobs = rep_id (min) min_time = time_mtch, by(confounding getallpossible ratio)

list

export excel using "${resultsDir}\computation_performance.xlsx", firstrow(variables) replace


//////////////

use "${resultsDir}\rep_cohortdeets_2000_phase1_ss.dta", clear

keep if time_mtch ==.
duplicates drop

collapse (mean) median_set set_25 set_75 full_sets (count) nobs = rep_id, by(confounding getallpossible ratio)
 
list

export excel using "${resultsDir}\matched_sets.xlsx", firstrow(variables) replace