
set scheme white_tableau 



foreach conf in  "10_nconf" "10_sconf" "10_bconf" "1_nconf" "1_sconf" "1_bconf" {	
	
use "$resultsDir\sim_summary_simsummtch1_${nsim}_${type}", clear

gen prev = "10% prevalence" if strpos(confounding, "10") > 0
replace prev = "1% prevalence" if strpos(confounding, "1") > 0 & strpos(confounding, "0") == 0
//replace prev = "7% prevalence" if strpos(confounding, "7") > 0

keep if test == "Relative % error in ModSE" 

keep if confounding == "`conf'" 

local prev = prev[1]

gen ratio_offset_1 = ratio + 0.1
gen ratio_offset_2 = ratio + 0.2
gen ratio_offset_3 = ratio + 0.3

gen ratio_offset_4 = ratio - 0.1
gen ratio_offset_5 = ratio - 0.2
gen ratio_offset_6 = ratio - 0.3

//rel_err

twoway connected avg ratio_offset_1 if method ==5 & getallpossible == "noreplacement", mcol("254 80 0") lcol("254 80 0")  ///
|| rcap MCSE_lo MCSE_hi ratio_offset_1 if method == 5 & getallpossible == "noreplacement", mcol("254 80 0") lcol("254 80 0") ///
|| connected avg ratio_offset_2 if method ==6 & getallpossible == "noreplacement", mcol("0 174 199") lcol("0 174 199") ///
|| rcap MCSE_lo MCSE_hi ratio_offset_2 if method == 6 & getallpossible == "noreplacement", mcol("0 174 199") lcol("0 174 199") ///
|| connected avg ratio_offset_3 if method ==7 & getallpossible == "noreplacement", mcol("0 191 111") lcol("0 191 111")   ///
|| rcap MCSE_lo MCSE_hi ratio_offset_3 if method == 7 & getallpossible == "noreplacement",  mcol("0 191 111") lcol("0 191 111") ///
///
|| connected avg ratio_offset_4 if method ==5 & getallpossible == "replacement",  mcol("254 80 0") lcol("254 80 0") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio_offset_4 if method == 5 & getallpossible == "replacement",  mcol("254 80 0") lcol("254 80 0") lpattern(dash) ///
|| connected avg ratio_offset_5 if method == 6 & getallpossible == "replacement", mcol("0 174 199") lcol("0 174 199") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio_offset_5 if method == 6 & getallpossible == "replacement", mcol("0 174 199") lcol("0 174 199") lpattern(dash) ///
|| connected avg ratio_offset_6 if method ==7 & getallpossible == "replacement",mcol("0 191 111") lcol("0 191 111") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio_offset_6 if method == 7 & getallpossible == "replacement", mcol("0 191 111") lcol("0 191 111") lpattern(dash) yline(0, lcol(black))  ///
	xtitle("") title("Cox models", size(small)) ///
	xlab(1 "1:1" 3 "3:1" 5 "5:1" 10 "10:1", labsize(medium)) ylab(, format(%9.0g) nolabels) name(rel_err_ca_norep, replace) ///
	legend(off)

twoway connected avg ratio_offset_1 if method ==11 & getallpossible == "noreplacement", mcol("254 80 0") lcol("254 80 0")  ///
|| rcap MCSE_lo MCSE_hi ratio_offset_1 if method == 11 & getallpossible == "noreplacement", mcol("254 80 0") lcol("254 80 0") ///
|| connected avg ratio_offset_2 if method ==12 & getallpossible == "noreplacement", mcol("0 174 199") lcol("0 174 199") ///
|| rcap MCSE_lo MCSE_hi ratio_offset_2 if method == 12 & getallpossible == "noreplacement", mcol("0 174 199") lcol("0 174 199") ///
|| connected avg ratio_offset_3 if method ==13 & getallpossible == "noreplacement", mcol("0 191 111") lcol("0 191 111")   ///
|| rcap MCSE_lo MCSE_hi ratio_offset_3 if method == 13 & getallpossible == "noreplacement",  mcol("0 191 111") lcol("0 191 111") ///
///
|| connected avg ratio_offset_4 if method ==11 & getallpossible == "replacement",  mcol("254 80 0") lcol("254 80 0") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio_offset_4 if method ==11 & getallpossible == "replacement",  mcol("254 80 0") lcol("254 80 0") lpattern(dash) ///
|| connected avg ratio_offset_5 if method == 12 & getallpossible == "replacement", mcol("0 174 199") lcol("0 174 199") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio_offset_5 if method == 12 & getallpossible == "replacement", mcol("0 174 199") lcol("0 174 199") lpattern(dash) ///
|| connected avg ratio_offset_6 if method ==13 & getallpossible == "replacement",mcol("0 191 111") lcol("0 191 111") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio_offset_6 if method == 13 & getallpossible == "replacement", mcol("0 191 111") lcol("0 191 111") lpattern(dash) yline(0, lcol(black))  ///
	xtitle("") title("Cox models", size(small)) ///
	xlab(1 "1:1" 3 "3:1" 5 "5:1" 10 "10:1", labsize(medium)) ylab(, format(%9.0g) nolabels) name(rel_err_ct_norep, replace) ///
	legend(off)
	

///
twoway connected avg ratio_offset_1 if method ==8 & getallpossible == "noreplacement", mcol("254 80 0") lcol("254 80 0") ///
|| rcap MCSE_lo MCSE_hi ratio_offset_1 if method == 8 & getallpossible == "noreplacement",  mcol("254 80 0") lcol("254 80 0") ///
|| connected avg ratio_offset_2 if method == 9 & getallpossible == "noreplacement", mcol("0 174 199") lcol("0 174 199") ///
|| rcap MCSE_lo MCSE_hi ratio_offset_2 if method == 9 & getallpossible == "noreplacement",  mcol("0 174 199") lcol("0 174 199") ///
|| connected avg ratio_offset_3 if method ==10 & getallpossible == "noreplacement", mcol("0 191 111") lcol("0 191 111") ///
|| rcap MCSE_lo MCSE_hi ratio_offset_3 if method == 10 & getallpossible == "noreplacement", mcol("0 191 111") lcol("0 191 111") ///
///
|| connected avg ratio_offset_4 if method == 8 & getallpossible == "replacement", mcol("254 80 0") lcol("254 80 0") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio_offset_4 if method == 8 & getallpossible == "replacement",  mcol("254 80 0") lcol("254 80 0") lpattern(dash) ///
|| connected avg ratio_offset_5 if method == 9 & getallpossible == "replacement", mcol("0 174 199") lcol("0 174 199") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio_offset_5 if method == 9 & getallpossible == "replacement", mcol("0 174 199") lcol("0 174 199") lpattern(dash) ///
|| connected avg ratio_offset_6 if method == 10 & getallpossible == "replacement", mcol("0 191 111") lcol("0 191 111") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio_offset_6 if method == 10 & getallpossible == "replacement", mcol("0 191 111") lcol("0 191 111") lpattern(dash) yline(0, lcol(black)) ///
	xtitle("") title("Poisson models", size(small)) ///
	xlab(1 "1:1" 3 "3:1" 5 "5:1" 10 "10:1", labsize(medium)) ylab(, format(%9.0g) labsize(medium)) 	name(rel_err_pa_norep, replace) ///
	legend(off)

///
twoway connected avg ratio_offset_1 if method == 14 & getallpossible == "noreplacement", mcol("254 80 0") lcol("254 80 0") ///
|| rcap MCSE_lo MCSE_hi ratio_offset_1 if method == 14 & getallpossible == "noreplacement",  mcol("254 80 0") lcol("254 80 0") ///
|| connected avg ratio_offset_2 if method == 15 & getallpossible == "noreplacement", mcol("0 174 199") lcol("0 174 199") ///
|| rcap MCSE_lo MCSE_hi ratio_offset_2 if method == 15 & getallpossible == "noreplacement",  mcol("0 174 199") lcol("0 174 199") ///
|| connected avg ratio_offset_3 if method ==16 & getallpossible == "noreplacement", mcol("0 191 111") lcol("0 191 111") ///
|| rcap MCSE_lo MCSE_hi ratio_offset_3 if method == 16 & getallpossible == "noreplacement", mcol("0 191 111") lcol("0 191 111") ///
///
|| connected avg ratio_offset_4 if method == 14 & getallpossible == "replacement", mcol("254 80 0") lcol("254 80 0") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio_offset_4 if method == 14 & getallpossible == "replacement",  mcol("254 80 0") lcol("254 80 0") lpattern(dash) ///
|| connected avg ratio_offset_5 if method == 15 & getallpossible == "replacement", mcol("0 174 199") lcol("0 174 199") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio_offset_5 if method == 15 & getallpossible == "replacement", mcol("0 174 199") lcol("0 174 199") lpattern(dash) ///
|| connected avg ratio_offset_6 if method == 16 & getallpossible == "replacement", mcol("0 191 111") lcol("0 191 111") lpattern(dash) msymbol(triangle) ///
|| rcap MCSE_lo MCSE_hi ratio_offset_6 if method == 16 & getallpossible == "replacement", mcol("0 191 111") lcol("0 191 111") lpattern(dash) yline(0, lcol(black)) ///
	xtitle("") title("Poisson models", size(small)) ///
	xlab(1 "1:1" 3 "3:1" 5 "5:1" 10 "10:1", labsize(medium)) ylab(, format(%9.0g) labsize(medium)) 	name(rel_err_pt_norep, replace) ///
	legend(off)
		
	
//graph save "rel_err_c_${getallpossible}_${nreps}.gph",   replace
//graph export "rel_err_c_${getallpossible}_${nreps}.emf", replace

graph combine rel_err_pt_norep rel_err_ct_norep, imargin(0 0 0 0) xcommon ycommon title("Time in study", size(small)) name(rel_err_t_`conf'_ratio, replace) //legendfrom(rel_err_p_${getallpossible})

graph combine rel_err_pa_norep rel_err_ca_norep, imargin(0 0 0 0) xcommon ycommon title("Age", size(small)) name(rel_err_a_`conf'_ratio, replace) //legendfrom(rel_err_p_${getallpossible})

graph combine rel_err_t_`conf'_ratio rel_err_a_`conf'_ratio, name(rel_err_`conf'_ratio, replace) title("Relative % error in Model SE") row(2) cols(2)

graph export "$resultsDir\sim_graph_rel_err_`conf'_ratio_${nsim}_${type}.emf", replace
}


//ADD IN GRAPHS IF YOU NEED MULTIPLE GRAPHS BY CONFOUNDING


 graph combine rel_err_t_10_sconf_ratio rel_err_a_10_sconf_ratio rel_err_t_1_sconf_ratio rel_err_a_1_sconf_ratio, xcommon ycommon imargin(0 0 0 0) c(2) r(2) name(rel_err_all_sconf, replace) title("Relative % error of the Model SE", size(small))
 
graph export "$resultsDir\sim_graph_rel_err_ratio_all_sconf_${nsim}_${type}.emf", replace


graph combine rel_err_t_10_nconf_ratio rel_err_a_10_nconf_ratio rel_err_t_1_nconf_ratio rel_err_a_1_nconf_ratio, xcommon ycommon imargin(0 0 0 0) c(2) r(2) name(rel_err_all_nconf, replace) title("Relative % error of the Model SE", size(small))
 
 
graph export "$resultsDir\sim_graph_rel_err_ratio_all_nconf_${nsim}_${type}.emf", replace


graph combine rel_err_t_10_bconf_ratio rel_err_a_10_bconf_ratio rel_err_t_1_bconf_ratio rel_err_a_1_bconf_ratio, xcommon ycommon imargin(0 0 0 0) c(2) r(2) name(rel_err_all_bconf, replace) title("Relative % error of the Model SE", size(small))
 
 
graph export "$resultsDir\sim_graph_rel_err_ratio_all_bconf_${nsim}_${type}.emf", replace