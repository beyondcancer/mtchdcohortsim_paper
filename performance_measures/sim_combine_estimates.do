cd $rawDir

local bc_est: dir . files "bc_mysimestimates*.dta"
clear
append using `bc_est'
save "bc_estimates_${nsim}_${type}.dta", replace
clear


local rep_est: dir . files "rep_mysimestimates*.dta"
clear
append using `rep_est'
save "rep_estimates_${nsim}_${type}.dta", replace
clear

use "bc_estimates_${nsim}_${type}.dta", clear
list in 1/10

append using "rep_estimates_${nsim}_${type}.dta"

save "${resultsDir}\sim_estimates_${nsim}_${type}", replace

local rep_deets: dir . files "rep_cohortdeets*.dta"
clear
append using `rep_deets'
save "${resultsDir}\rep_cohortdeets_${nsim}_${type}.dta", replace
clear


/*
local bc_state: dir . files "bc_mysimstate*.dta"
clear
append using `bc_state'
save "${resultsDir}\sim_bc_states_${nsim}_${type}.dta", replace
clear

local rep_state: dir . files "rep_mysimstate*.dta"
clear
append using `rep_state'
save "${resultsDir}\sim_rep_states_${nsim}_${type}.dta", replace
*/
clear