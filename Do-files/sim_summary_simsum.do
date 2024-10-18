

use "$resultsDir\sim_estimates_${nsim}_${type}", clear

//drop if confounding != "no_conf" | ratio != 1 | getallpossible != "noreplacement"

sort rep_id

//For each simulation:
replace getallpossible = "Base cohort" if getallpossible == ""

preserve

drop if getallpossible == "Base cohort"

simsum theta, se(se) methodvar(method) id(rep_id) true($true) mcse format(%6.3f %6.1f %6.0f) bias empse modelse relerror cover listsep by(confounding ratio getallpossible) saving($resultsDir\sim_summary_simsummtch_${nsim}_${type}, replace)

restore


preserve

drop if getallpossible != "Base cohort"

simsum theta, se(se) methodvar(method) id(rep_id) true($true) mcse format(%6.3f %6.1f %6.0f) listsep bias empse modelse relerror cover by(confounding) saving($resultsDir\sim_summary_simsumbc_${nsim}_${type}, replace)





use "$resultsDir\sim_summary_simsummtch_${nsim}_${type}", clear


gen type = _n


forvalues i = 3/8 {
gen MCSE_low`i'= theta`i' - (1.96 * theta`i'_mcse)
gen MCSE_hi`i'= theta`i' + (1.96 * theta`i'_mcse)

drop theta`i'_mcse
}


reshape long theta MCSE_low MCSE_hi, i(type) j(method)

label define methodlab 1 "base_uadj" 2 "base_adj"  3 "cox_uadj" 4 "cox_adj" 5 "cox_match" 6 "p_uadj" 7 "p_adj" 8"p_match"
label values method methodlab

rename theta avg
rename statcode test

replace test = "Bias" if test == "bias"
replace test = "Coverage of 95% CIs" if test == "cover"
replace test = "Empirical SE" if test == "empse"
replace test = "Relative % error in ModSE" if test == "relerror"

drop type
drop statnum

save "$resultsDir\sim_summary_simsummtch1_${nsim}_${type}", replace





use "$resultsDir\sim_summary_simsumbc_${nsim}_${type}", clear


gen type = _n


forvalues i = 1/2 {
gen MCSE_low`i'= theta`i' - (1.96 * theta`i'_mcse)
gen MCSE_hi`i'= theta`i' + (1.96 * theta`i'_mcse)

drop theta`i'_mcse
}


reshape long theta MCSE_low MCSE_hi, i(type) j(method)

label define methodlab 1 "base_uadj" 2 "base_adj"  3 "cox_uadj" 4 "cox_adj" 5 "cox_match" 6 "p_uadj" 7 "p_adj" 8"p_match"
label values method methodlab

rename theta avg
rename statcode test

replace test = "Bias" if test == "bias"
replace test = "Model Based SE" if test == "modelse"
replace test = "Coverage of 95% CIs" if test == "cover"
replace test = "Empirical SE" if test == "empse"

drop type


save "$resultsDir\sim_summary_simsumbc1_${nsim}_${type}", replace
