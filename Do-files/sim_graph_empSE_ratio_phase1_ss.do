
set scheme gg_tableau 



foreach conf in "10_noconf" "10_smallco" "7_noconf" "7_smallcon" "5_noconf" "5_smallcon" {
	
use "$resultsDir\sim_summary_simsummtch1_${nsim}_${type}", clear

gen prev = "10% prevalence" if strpos(confounding, "10") > 0
replace prev = "5% prevalence" if strpos(confounding, "5") > 0
replace prev = "7% prevalence" if strpos(confounding, "7") > 0

keep if test == "Empirical SE" 

keep if confounding == "`conf'" 

local prev = prev[1]

cap gen ratio_offset_1 = ratio + 0.08
cap gen ratio_offset_2 = ratio + 0.08

//empSE

twoway connected avg ratio_offset_1 if method ==6 & getallpossible == "noreplacement", mcol("121 153 157") lcol("121 153 157")  ///
|| rcap MCSE_lo MCSE_hi ratio_offset_1 if method == 6 & getallpossible == "noreplacement", mcol("121 153 157") lcol("121 153 157") ///
|| connected avg ratio_offset_2 if method ==7 & getallpossible == "noreplacement", mcol("39 182 122") lcol("39 182 122") ///
|| rcap MCSE_lo MCSE_hi ratio_offset_2 if method == 7 & getallpossible == "noreplacement", mcol("39 182 122") lcol("39 182 122") ///
|| connected avg ratio if method ==8 & getallpossible == "noreplacement", mcol("0 70 90") lcol("0 70 90")   ///
|| rcap MCSE_lo MCSE_hi ratio if method == 8 & getallpossible == "noreplacement",  mcol("0 70 90") lcol("0 70 90") ///
///
|| connected avg ratio_offset_2 if method ==6 & getallpossible == "replacement",  mcol("121 153 157") lcol("121 153 157") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio_offset_2 if method == 6 & getallpossible == "replacement",  mcol("121 153 157") lcol("121 153 157") lpattern(dash) ///
|| connected avg ratio if method == 7 & getallpossible == "replacement", mcol("39 182 122") lcol("39 182 122") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio if method == 7 & getallpossible == "replacement", mcol("39 182 122") lcol("39 182 122") lpattern(dash) ///
|| connected avg ratio if method ==8 & getallpossible == "replacement",mcol("0 70 90") lcol("0 70 90") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio if method == 8 & getallpossible == "replacement", mcol("0 70 90") lcol("0 70 90") lpattern(dash)  ///
	xtitle("") title("Poisson models", size(small)) ///
	xlab(1 "1:1" 3 "3:1" 5 "5:1" 10 "10:1", labsize(tiny)) ylab(0.02(0.01)0.1, labsize(tiny) ) name(empSE_p_norep, replace) ///
	legend(off)
	
	//(order( 1 "Unadjusted" 3 "Adjusted" 5 "Matched") subtitle("No replacement") ring(0) cols(3) pos(12)) ///


//graph save "empSE_p_${getallpossible}_${nreps}.gph",   replace
//graph export "empSE_p_${getallpossible}_${nreps}.emf", replace


///
twoway connected avg ratio_offset_1 if method ==3 & getallpossible == "noreplacement", mcol("121 153 157") lcol("121 153 157") ///
|| rcap MCSE_lo MCSE_hi ratio_offset_1 if method == 3 & getallpossible == "noreplacement",  mcol("121 153 157") lcol("121 153 157") ///
|| connected avg ratio if method == 4 & getallpossible == "noreplacement", mcol("39 182 122") lcol("39 182 122") ///
|| rcap MCSE_lo MCSE_hi ratio if method == 4 & getallpossible == "noreplacement",  mcol("39 182 122") lcol("39 182 122") ///
|| connected avg ratio if method ==5 & getallpossible == "noreplacement", mcol("0 70 90") lcol("0 70 90") ///
|| rcap MCSE_lo MCSE_hi ratio if method == 5 & getallpossible == "noreplacement", mcol("0 70 90") lcol("0 70 90") ///
///
|| connected avg ratio_offset_2 if method == 3 & getallpossible == "replacement", mcol("121 153 157") lcol("121 153 157") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio_offset_2 if method == 3 & getallpossible == "replacement",  mcol("121 153 157") lcol("121 153 157") lpattern(dash) ///
|| connected avg ratio if method == 4 & getallpossible == "replacement", mcol("39 182 122") lcol("39 182 122") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio if method == 4 & getallpossible == "replacement", mcol("39 182 122") lcol("39 182 122") lpattern(dash) ///
|| connected avg ratio if method ==5 & getallpossible == "replacement", mcol("0 70 90") lcol("0 70 90") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio if method == 5 & getallpossible == "replacement", mcol("0 70 90") lcol("0 70 90") lpattern(dash)  ///
	xtitle("") title("Cox models", size(small)) ///
	xlab(1 "1:1" 3 "3:1" 5 "5:1" 10 "10:1", labsize(tiny)) ylab(0.02(0.01) 0.1, noticks nolabels) name(empSE_c_norep, replace) ///
	legend(off)
	
	//(order( 7 "Unadjusted" 9 "Adjusted" 11 "Matched") subtitle("Replacement") ring(0) cols(3) pos(12)) ///

//graph save "empSE_c_${getallpossible}_${nreps}.gph",   replace
//graph export "empSE_c_${getallpossible}_${nreps}.emf", replace

graph combine empSE_p_norep empSE_c_norep, imargin(0 0 0 0) title("`prev'", size(small)) xcommon ycommon name(empSE_`conf'_ratio, replace) //legendfrom(empSE_p_${getallpossible})


graph export "$resultsDir\sim_graph_empSE_`conf'_ratio_${nsim}_${type}.emf", replace
}


graph combine empSE_10_smallco_ratio empSE_7_smallcon_ratio empSE_5_smallcon_ratio, xcommon ycommon imargin(0 0 0 0) c(3) name(empSE_all_smallconf, replace) title("Empirical SE", size(small))
 
graph export "$resultsDir\sim_graph_empSE_ratio_all_smallconf_${nsim}_${type}.emf", replace


graph combine empSE_10_noconf_ratio empSE_7_noconf_ratio empSE_5_noconf_ratio, xcommon ycommon  imargin(0 0 0 0) c(3) name(empSE_all_snoconf, replace) title("Empirical SE", size(small))
 
graph export "$resultsDir\sim_graph_empSE_ratio_all_noconf_${nsim}_${type}.emf", replace

