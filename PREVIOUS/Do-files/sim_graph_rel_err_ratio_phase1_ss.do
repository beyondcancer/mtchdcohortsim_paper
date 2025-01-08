
set scheme gg_tableau 



foreach conf in "1_smallcon" "1_noconf" "10_smallco" "10_noconf" {	
	
use "$resultsDir\sim_summary_simsummtch1_${nsim}_${type}", clear

gen prev = "10% prevalence" if strpos(confounding, "10") > 0
replace prev = "1% prevalence" if strpos(confounding, "1") > 0 & strpos(confounding, "0") == 0
//replace prev = "7% prevalence" if strpos(confounding, "7") > 0

keep if test == "Relative % error in ModSE" 

keep if confounding == "`conf'" 

local prev = prev[1]

gen ratio_offset_1 = ratio + 0.08
gen ratio_offset_2 = ratio + 0.08

//rel_err

twoway connected avg ratio_offset_1 if method ==5 & getallpossible == "noreplacement", mcol("121 153 157") lcol("121 153 157")  ///
|| rcap MCSE_lo MCSE_hi ratio_offset_1 if method == 5 & getallpossible == "noreplacement", mcol("121 153 157") lcol("121 153 157") ///
|| connected avg ratio_offset_2 if method ==6 & getallpossible == "noreplacement", mcol("39 182 122") lcol("39 182 122") ///
|| rcap MCSE_lo MCSE_hi ratio_offset_2 if method == 6 & getallpossible == "noreplacement", mcol("39 182 122") lcol("39 182 122") ///
|| connected avg ratio if method ==7 & getallpossible == "noreplacement", mcol("0 70 90") lcol("0 70 90")   ///
|| rcap MCSE_lo MCSE_hi ratio if method == 7 & getallpossible == "noreplacement",  mcol("0 70 90") lcol("0 70 90") ///
///
|| connected avg ratio_offset_2 if method ==5 & getallpossible == "replacement",  mcol("121 153 157") lcol("121 153 157") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio_offset_2 if method == 5 & getallpossible == "replacement",  mcol("121 153 157") lcol("121 153 157") lpattern(dash) ///
|| connected avg ratio if method == 6 & getallpossible == "replacement", mcol("39 182 122") lcol("39 182 122") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio if method == 6 & getallpossible == "replacement", mcol("39 182 122") lcol("39 182 122") lpattern(dash) ///
|| connected avg ratio if method ==7 & getallpossible == "replacement",mcol("0 70 90") lcol("0 70 90") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio if method == 7 & getallpossible == "replacement", mcol("0 70 90") lcol("0 70 90") lpattern(dash) yline(0, lcol(black))  ///
	xtitle("") title("Cox models", size(small)) ///
	xlab(1 "1:1" 3 "3:1" 5 "5:1" 10 "10:1") ylab(, nolabels) name(rel_err_ca_norep, replace) ///
	legend(off)
	
	//(order( 1 "Unadjusted" 3 "Adjusted" 5 "Matched") subtitle("No replacement") ring(0) cols(3) pos(12)) ///


//graph save "rel_err_p_${getallpossible}_${nreps}.gph",   replace
//graph export "rel_err_p_${getallpossible}_${nreps}.emf", replace


///
twoway connected avg ratio_offset_1 if method ==8 & getallpossible == "noreplacement", mcol("121 153 157") lcol("121 153 157") ///
|| rcap MCSE_lo MCSE_hi ratio_offset_1 if method == 8 & getallpossible == "noreplacement",  mcol("121 153 157") lcol("121 153 157") ///
|| connected avg ratio if method == 9 & getallpossible == "noreplacement", mcol("39 182 122") lcol("39 182 122") ///
|| rcap MCSE_lo MCSE_hi ratio if method == 9 & getallpossible == "noreplacement",  mcol("39 182 122") lcol("39 182 122") ///
|| connected avg ratio if method ==10 & getallpossible == "noreplacement", mcol("0 70 90") lcol("0 70 90") ///
|| rcap MCSE_lo MCSE_hi ratio if method == 10 & getallpossible == "noreplacement", mcol("0 70 90") lcol("0 70 90") ///
///
|| connected avg ratio_offset_2 if method == 8 & getallpossible == "replacement", mcol("121 153 157") lcol("121 153 157") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio_offset_2 if method == 8 & getallpossible == "replacement",  mcol("121 153 157") lcol("121 153 157") lpattern(dash) ///
|| connected avg ratio if method == 9 & getallpossible == "replacement", mcol("39 182 122") lcol("39 182 122") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio if method == 9 & getallpossible == "replacement", mcol("39 182 122") lcol("39 182 122") lpattern(dash) ///
|| connected avg ratio if method ==10 & getallpossible == "replacement", mcol("0 70 90") lcol("0 70 90") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio if method ==10 & getallpossible == "replacement", mcol("0 70 90") lcol("0 70 90") lpattern(dash) yline(0, lcol(black)) ///
	xtitle("") title("Poisson models", size(small)) ///
	xlab(1 "1:1" 3 "3:1" 5 "5:1" 10 "10:1") ylab(, ) 	name(rel_err_pa_norep, replace) ///
	legend(off)
	
	//(order( 7 "Unadjusted" 9 "Adjusted" 11 "Matched") subtitle("Replacement") ring(0) cols(3) pos(12)) ///



twoway connected avg ratio_offset_1 if method ==11 & getallpossible == "noreplacement", mcol("121 153 157") lcol("121 153 157")  ///
|| rcap MCSE_lo MCSE_hi ratio_offset_1 if method == 11 & getallpossible == "noreplacement", mcol("121 153 157") lcol("121 153 157") ///
|| connected avg ratio_offset_2 if method ==12 & getallpossible == "noreplacement", mcol("39 182 122") lcol("39 182 122") ///
|| rcap MCSE_lo MCSE_hi ratio_offset_2 if method == 12 & getallpossible == "noreplacement", mcol("39 182 122") lcol("39 182 122") ///
|| connected avg ratio if method ==13 & getallpossible == "noreplacement", mcol("0 70 90") lcol("0 70 90")   ///
|| rcap MCSE_lo MCSE_hi ratio if method == 13 & getallpossible == "noreplacement",  mcol("0 70 90") lcol("0 70 90") ///
///
|| connected avg ratio_offset_2 if method ==11 & getallpossible == "replacement",  mcol("121 153 157") lcol("121 153 157") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio_offset_2 if method ==11 & getallpossible == "replacement",  mcol("121 153 157") lcol("121 153 157") lpattern(dash) ///
|| connected avg ratio if method == 12 & getallpossible == "replacement", mcol("39 182 122") lcol("39 182 122") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio if method == 12 & getallpossible == "replacement", mcol("39 182 122") lcol("39 182 122") lpattern(dash) ///
|| connected avg ratio if method ==13 & getallpossible == "replacement",mcol("0 70 90") lcol("0 70 90") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio if method == 13 & getallpossible == "replacement", mcol("0 70 90") lcol("0 70 90") lpattern(dash) yline(0, lcol(black))  ///
	xtitle("") title("Cox models", size(small)) ///
	xlab(1 "1:1" 3 "3:1" 5 "5:1" 10 "10:1") ylab(, nolabels) name(rel_err_ct_norep, replace) ///
	legend(off)
	
	//(order( 1 "Unadjusted" 3 "Adjusted" 5 "Matched") subtitle("No replacement") ring(0) cols(3) pos(12)) ///


//graph save "rel_err_p_${getallpossible}_${nreps}.gph",   replace
//graph export "rel_err_p_${getallpossible}_${nreps}.emf", replace


///
twoway connected avg ratio_offset_1 if method ==14 & getallpossible == "noreplacement", mcol("121 153 157") lcol("121 153 157") ///
|| rcap MCSE_lo MCSE_hi ratio_offset_1 if method == 14 & getallpossible == "noreplacement",  mcol("121 153 157") lcol("121 153 157") ///
|| connected avg ratio if method == 15 & getallpossible == "noreplacement", mcol("39 182 122") lcol("39 182 122") ///
|| rcap MCSE_lo MCSE_hi ratio if method == 15 & getallpossible == "noreplacement",  mcol("39 182 122") lcol("39 182 122") ///
|| connected avg ratio if method ==16 & getallpossible == "noreplacement", mcol("0 70 90") lcol("0 70 90") ///
|| rcap MCSE_lo MCSE_hi ratio if method == 16 & getallpossible == "noreplacement", mcol("0 70 90") lcol("0 70 90") ///
///
|| connected avg ratio_offset_2 if method == 14 & getallpossible == "replacement", mcol("121 153 157") lcol("121 153 157") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio_offset_2 if method == 14 & getallpossible == "replacement",  mcol("121 153 157") lcol("121 153 157") lpattern(dash) ///
|| connected avg ratio if method == 15 & getallpossible == "replacement", mcol("39 182 122") lcol("39 182 122") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio if method == 15 & getallpossible == "replacement", mcol("39 182 122") lcol("39 182 122") lpattern(dash) ///
|| connected avg ratio if method ==16 & getallpossible == "replacement", mcol("0 70 90") lcol("0 70 90") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio if method ==16 & getallpossible == "replacement", mcol("0 70 90") lcol("0 70 90") lpattern(dash) yline(0, lcol(black)) ///
	xtitle("") title("Poisson models", size(small)) ///
	xlab(1 "1:1" 3 "3:1" 5 "5:1" 10 "10:1") ylab(, ) 	name(rel_err_pt_norep, replace) ///
	legend(off)
		
	
//graph save "rel_err_c_${getallpossible}_${nreps}.gph",   replace
//graph export "rel_err_c_${getallpossible}_${nreps}.emf", replace

graph combine rel_err_pt_norep rel_err_ct_norep, imargin(0 0 0 0) xcommon ycommon title("`prev' - time", size(small)) name(rel_err_t_`conf'_ratio, replace) //legendfrom(rel_err_p_${getallpossible})

graph combine rel_err_pa_norep rel_err_ca_norep, imargin(0 0 0 0) xcommon ycommon title("`prev' - age", size(small)) name(rel_err_a_`conf'_ratio, replace) //legendfrom(rel_err_p_${getallpossible})


graph export "$resultsDir\sim_graph_rel_err_`conf'_ratio_${nsim}_${type}.emf", replace
}

//ADD IN GRAPHS IF YOU NEED MULTIPLE GRAPHS BY CONFOUNDING

//graph combine rel_err_10_smallco_ratio rel_err_7_smallcon_ratio rel_err_5_smallcon_ratio, xcommon ycommon imargin(0 0 0 0) c(3) name(rel_err_all_smallconf, replace) title("Relative % error of the Model SE", size(small))

 graph combine rel_err_t_10_smallco_ratio rel_err_a_10_smallco_ratio rel_err_t_1_smallcon_ratio rel_err_a_1_smallcon_ratio, xcommon ycommon imargin(0 0 0 0) c(2) r(2) name(rel_err_all_smallconf, replace) title("Relative % error of the Model SE", size(small))
 
graph export "$resultsDir\sim_graph_rel_err_ratio_all_smallconf_${nsim}_${type}.emf", replace


graph combine rel_err_t_10_noconf_ratio rel_err_a_10_noconf_ratio rel_err_t_1_noconf_ratio rel_err_a_1_noconf_ratio, xcommon ycommon imargin(0 0 0 0) c(2) r(2) name(rel_err_all_noconf, replace) title("Relative % error of the Model SE", size(small))

graph export "$resultsDir\sim_graph_rel_err_ratio_all_noconf_${nsim}_${type}.emf", replace