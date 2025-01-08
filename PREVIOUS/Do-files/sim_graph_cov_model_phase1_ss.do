set scheme gg_tableau

use "$resultsDir\sim_summary_simsummtch1_${nsim}_${type}", clear

drop if test != "Coverage of 95% CIs"
label define ratio_lab 1 "ratio 1:1" 3 "ratio 3:1" 5 "ratio 5:1" 10 "ratio 10:1"
label val ratio ratio_lab

gen prev = "10% prevalence" if strpos(confounding, "10") > 0
replace prev = "1% prevalence" if strpos(confounding, "1") > 0 & strpos(confounding, "0") == 0
//replace prev = "7% prevalence" if strpos(confounding, "7") > 0

//foreach conf in "no_conf" "small_conf" "mid_conf" "big_conf"{
	
//foreach conf in "10_noconf" "10_smallco" "7_noconf" "7_smallcon" "5_noconf" "5_smallcon" {
foreach conf in "1_smallcon" "1_noconf" "10_smallco" "10_noconf" {	

	
preserve

drop if confounding!= `"`conf'"'

local prev = prev[1]

    graph twoway scatter avg method if getallpossible == "replacement", mcol("121 153 157") ///
|| rcap MCSE_low MCSE_hi method if getallpossible == "replacement", lcol("121 153 157")  ///
|| scatter avg method if getallpossible == "noreplacement", mcol("0 70 90") /// 
|| rcap MCSE_low MCSE_hi method if getallpossible == "noreplacement", lcol("0 70 90") ///
xlabel(	5 "cox_uadj_a" 6 "cox_adj_a" 7 "cox_match_a" 8 "p_uadj_a" 9 "p_adj_a" 10 "p_match_a" ///
							11 "cox_uadj_t" 12 "cox_adj_t" 13 "cox_match_t" 14 "p_uadj_t" 15 "p_adj_t" 16 "p_match_t", labsize(2) angle(90)) yline(95, lcol(black)) xtitle("") ///
 by(ratio, noyrescale note("") legend(off) title(`prev', size(vsmall))) ylab(, format(%9.2f) labsize(tiny)) name(cov_`conf'_model, replace)

restore 

graph export "$resultsDir\sim_graph_cov_model_`conf'_${nsim}_${type}.emf", replace

}


graph combine cov_10_smallco_model cov_1_smallcon_model, xcommon ycommon imargin(0 0 0 0) name(cov_all_smallconf, replace)
 
 
graph export "$resultsDir\sim_graph_cov_model_all_smallconf_${nsim}_${type}.emf", replace


graph combine cov_10_noconf_model cov_1_noconf_model, xcommon ycommon imargin(0 0 0 0) name(cov_all_noconf, replace)
 
graph export "$resultsDir\sim_graph_cov_model_all_noconf_${nsim}_${type}.emf", replace