set scheme white_tableau

use "$resultsDir\sim_summary_simsummtch1_${nsim}_${type}", clear

drop if test != "Bias"
label define ratio_lab 1 "ratio 1:1" 3 "ratio 3:1" 5 "ratio 5:1" 10 "ratio 10:1"
label val ratio ratio_lab

gen prev = "10% prevalence" if strpos(confounding, "10") > 0
replace prev = "1% prevalence" if strpos(confounding, "1") > 0 & strpos(confounding, "0") == 0



gen method_reordered = .
replace method_reordered = 1 if method == 5  // cox_uadj_a
replace method_reordered = 2 if method == 11 // cox_uadj_t
replace method_reordered = 3 if method == 8  // p_uadj_a
replace method_reordered = 4 if method == 14 // p_uadj_t
replace method_reordered = 5 if method == 6  // cox_adj_a
replace method_reordered = 6 if method == 12 // cox_adj_t
replace method_reordered = 7 if method == 9  // p_adj_a
replace method_reordered = 8 if method == 15 // p_adj_t
replace method_reordered = 9 if method == 7  // cox_match_a
replace method_reordered = 10 if method == 13 // cox_match_t
replace method_reordered = 11 if method == 10 // p_match_a
replace method_reordered = 12 if method == 16 // p_match_t

gen method_offset_1 = method_reordered + 0.08
gen method_offset_2 = method_reordered - 0.08


//replace prev = "7% prevalence" if strpos(confounding, "7") > 0


foreach conf in  "10_nconf" "10_sconf" "10_bconf" "1_nconf" "1_sconf" "1_bconf"  {	

	
preserve

drop if confounding!= `"`conf'"'

local prev = prev[1]

    graph twoway scatter avg method_offset_1 if getallpossible == "replacement", mcol("255 184 28") msize(medium) ///
|| rcap MCSE_low MCSE_hi method_offset_1 if getallpossible == "replacement", lcol("255 184 28") lwidth(medium) ///
|| scatter avg method_offset_2 if getallpossible == "noreplacement", mcol("75 0 45") msize(medium) /// 
|| rcap MCSE_low MCSE_hi method_offset_2 if getallpossible == "noreplacement", lcol("75 0 45") lwidth(medium) ///
xlabel(	1 "Cox-UA (Age)" 2 "Cox-UA (Time)" 3 "Pois-UA (Age)" 4 "Pois-UA (Time)" ///
       5 "Cox-A (Age)" 6 "Cox-A (Time)" 7 "Pois-A (Age)" 8 "Pois-A (Time)" ///
       9 "Cox-Strat (Age)" 10 "Cox-Strat (Time)" 11 "Pois-Cond (Age)" 12 "Pois-Cond (Time)", angle(45) labsize(small)) yline(0, lcol(black)) xtitle("") ///
 by(ratio, noyrescale imargin(0.6 0 0 0) note("") legend(off) title("Bias of the estimate", size(medium))) ylab(-0.04(0.02)0.1, format(%9.2f) labsize(small))  name(bias_`conf'_model, replace)

restore 

graph export "$resultsDir\sim_graph_bias_model_`conf'_${nsim}_${type}.emf", replace

}



graph combine bias_10_sconf_model bias_1_sconf_model, xcommon ycommon imargin(0 0 0 0) title("Bias of the estimate",size(medium)) name(bias_all_smallconf, replace)
 
graph export "$resultsDir\sim_graph_bias_model_all_smallconf_${nsim}_${type}.emf", replace


graph combine bias_10_nconf_model bias_1_nconf_model, xcommon ycommon imargin(0 0 0 0) name(bias_all_noconf, replace)
 
graph export "$resultsDir\sim_graph_bias_model_all_noconf_${nsim}_${type}.emf", replace

graph combine bias_10_bconf_model bias_1_bconf_model, xcommon ycommon imargin(0 0 0 0) name(bias_all_bconf, replace)
 
graph export "$resultsDir\sim_graph_bias_model_all_bconf_${nsim}_${type}.emf", replace

