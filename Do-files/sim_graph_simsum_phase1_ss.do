set scheme white_tableau

use "$resultsDir\sim_summary_simsummtch1_${nsim}_${type}", clear

drop if test == "modelse"

foreach conf in "1_smallcon" "1_noconf" "10_smallco" "10_noconf" {

	
preserve

drop if confounding!= `"`conf'"'

    graph twoway scatter avg method if getallpossible == "replacement", mcol("255 127 14") ///
|| rcap MCSE_low MCSE_hi method if getallpossible == "replacement", lcol("255 127 14") ///
|| scatter avg method if getallpossible == "noreplacement", mcol("44 160 44") /// 
|| rcap MCSE_low MCSE_hi method if getallpossible == "noreplacement", lcol("44 160 44") ///
xlabel(5 "cox_uadj_a" 6 "cox_adj_a" 7 "cox_match_a" 8 "p_uadj_a" 9 "p_adj_a" 10 "p_match_a" ///
	11 "cox_uadj_t" 12 "cox_adj_t" 13 "cox_match_t" 14 "p_uadj_t" 15 "p_adj_t" 16 "p_match_t", angle(vertical)) xtitle("") ///
 by(test ratio, yrescale note("") compact col(4) noedgelabel legend(off) title(`conf')) ylab(, format(%9.2f) labsize(vsmall)) name(test_`conf'_${type}, replace)

restore 


gr_edit .plotregion1.yaxis1[1].reset_rule -0.01 0.03 0.015 , tickset(major) ruletype(range) 
gr_edit .plotregion1.yaxis1[2].reset_rule -0.01 0.03 0.015 , tickset(major) ruletype(range) 
gr_edit .plotregion1.yaxis1[3].reset_rule -0.01 0.03 0.015 , tickset(major) ruletype(range) 
gr_edit .plotregion1.yaxis1[4].reset_rule -0.01 0.03 0.015 , tickset(major) ruletype(range) 

gr_edit .plotregion1.yaxis1[5].reset_rule 70 100 10 , tickset(major) ruletype(range) 
gr_edit .plotregion1.yaxis1[6].reset_rule 70 100 10 , tickset(major) ruletype(range) 
gr_edit .plotregion1.yaxis1[7].reset_rule 70 100 10 , tickset(major) ruletype(range) 
gr_edit .plotregion1.yaxis1[8].reset_rule 70 100 10 , tickset(major) ruletype(range)

gr_edit .plotregion1.yaxis1[9].reset_rule 0.04 0.08 0.015, tickset(major) ruletype(range) 
gr_edit .plotregion1.yaxis1[10].reset_rule 0.04 0.08 0.015  , tickset(major) ruletype(range) 
gr_edit .plotregion1.yaxis1[11].reset_rule 0.04 0.08 0.015  , tickset(major) ruletype(range) 
gr_edit .plotregion1.yaxis1[12].reset_rule 0.04 0.08 0.015  , tickset(major) ruletype(range) 

gr_edit .plotregion1.yaxis1[13].reset_rule -50 5 15  , tickset(major) ruletype(range) 
gr_edit .plotregion1.yaxis1[14].reset_rule -50 5 15  , tickset(major) ruletype(range) 
gr_edit .plotregion1.yaxis1[15].reset_rule -50 5 15  , tickset(major) ruletype(range) 
gr_edit .plotregion1.yaxis1[16].reset_rule -50 5 15  , tickset(major) ruletype(range) 


graph play "H:\mtchdcohortsim\remove_yaxis.grec"
graph play "H:\mtchdcohortsim\add_lines.grec"
 
graph export "$resultsDir\summary_`conf'_${nsim}_${type}.emf", replace

}



