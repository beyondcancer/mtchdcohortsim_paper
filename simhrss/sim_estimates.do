//quietly{
foreach l in $conf_test {
	local exp_ithbi = word("$conf_expbi", `l')
	local out_ithbi = word("$conf_outbi", `l')
	local exp_ithcont = word("$conf_expcont", `l')
	local out_ithcont = word("$conf_outcont", `l')
	local type_ith = word("$conf_type", `l')
	local lambda_ith = word("$lambda_exp", `l')
noi disp "Initiating the `type_ith' simulation loop, this is loop `l'"

 tempname bc_estimates rep_estimates	
 postfile `bc_estimates' float(nobs) str10 confounding int(rep_id) int(ratio) int(method) float(theta se lci uci pval) using bc_mysimestimates_`type_ith'_sim${nsim_start}_${nsim}.dta, replace	
postfile `rep_estimates' str10 confounding str20 getallpossible int(rep_id) int(ratio) int(method) float(nomatch theta se lci uci pval) using rep_mysimestimates_`type_ith'_sim${nsim_start}_${nsim}.dta, replace 

tempname bc_states rep_states	
postfile `bc_states' str10 confounding int(rep_id) str2000 bc1 str2000 bc2 str2000 bc3 using bc_mysimstates_`type_ith'_sim${nsim_start}_${nsim}, replace
postfile `rep_states' str10 confounding int(rep_id) str20 implementation int(ratio) str2000 state1 str2000 state2 str2000 state3 using rep_mysimstates_`type_ith'_sim${nsim_start}_${nsim}.dta, replace ///figure this out later

forval i= $nsim_start/$nsim {	
	
noi disp "this is rep `i' out of " ($nsim - $nsim_start + 1)

create_bc, nobs($n) loghr($loghr) i(`i') lambda_exp(`lambda_ith') age_exp(`exp_ithcont') age_out(`out_ithcont') gender_exp(`exp_ithbi') gender_out(`out_ithbi') conf_type("`type_ith'")

post `bc_states' ("`type_ith'") (`i') (substr(r(bc_state),1,2000)) (substr(r(bc_state),2001, 2000)) (substr(r(bc_state),4001,.)) 

post `bc_estimates' (r(total_n)) ("`type_ith'") (`i') (0) (1) (r(base_uadj_est)) (r(base_uadj_se)) (r(base_uadj_lci)) (r(base_uadj_uci)) (r(base_uadj_pval))
post `bc_estimates' (r(total_n)) ("`type_ith'") (`i') (0) (2) (r(base_adj_est)) (r(base_adj_se)) (r(base_adj_lci)) (r(base_adj_uci)) (r(base_adj_pval))

noi disp "base cohort created : rep `i' `type_ith' lambda_exp is `lambda_ith'"	

///Info on the cohort
		count if tag 
		qui local total_n = r(N)
		noi disp "Total patients: " `total_n' 
		tab exposed outcome, row
		tabulate exposed, matcell(freq_matrix) matrow(var_values)
		local exp_n = freq_matrix[2,1]
		noi disp "Total exposed: " `exp_n' 
		noi disp (`exp_n'/`total_n'*100) "% of exposed"
		
		
///Matching loops
foreach rep of global implementation {	
		 noi disp "implementation = `rep'"
		foreach j of global ratio { // iteration of matching ratio
			noi disp "ratio = `j'"	
			

		matched_sim, getallpossible(`rep') ratio(`j') i(`i') conf_type("`type_ith'")
		
		post `rep_states' ("`type_ith'") (`i') ("`rep'") (`j') (substr(r(state_`rep'_`j'),1,2000)) (substr(r(state_`rep'_`j'),2001, 2000)) (substr(r(state_`rep'_`j'),4001,.)) 

		post `rep_estimates' ("`type_ith'") ("`rep'") (`i') (`j') (3) (r(nomatch)) (r(cox_uadj_est)) (r(cox_uadj_se))  (r(cox_uadj_lci)) (r(cox_uadj_uci))  (r(cox_uadj_pval))
		post `rep_estimates' ("`type_ith'") ("`rep'") (`i') (`j') (4) (r(nomatch)) (r(cox_adj_est))  (r(cox_adj_se))  (r(cox_adj_lci)) (r(cox_adj_uci))  (r(cox_adj_pval))
		post `rep_estimates' ("`type_ith'") ("`rep'") (`i') (`j') (5) (r(nomatch)) (r(cox_mtch_est)) (r(cox_mtch_se))  (r(cox_mtch_lci)) (r(cox_mtch_uci))  (r(cox_mtch_pval))
		post `rep_estimates' ("`type_ith'") ("`rep'") (`i') (`j') (6) (r(nomatch)) (r(p_uadj_est))  (r(p_uadj_se))  (r(p_uadj_lci)) (r(p_uadj_uci))  (r(p_uadj_pval))
		post `rep_estimates' ("`type_ith'") ("`rep'") (`i') (`j') (7) (r(nomatch)) (r(p_adj_est))  (r(p_adj_se))  (r(p_adj_lci)) (r(p_adj_uci))  (r(p_adj_pval))
		post `rep_estimates' ("`type_ith'") ("`rep'") (`i') (`j') (8) (r(nomatch)) (r(p_mtch_est))  (r(p_mtch_se))  (r(p_mtch_lci)) (r(p_mtch_uci))  (r(p_mtch_pval))
		
return list
noi disp "erasing `rep' `type_ith' matched cohort data"
noi erase "getmatchedcohort`i'`type_ith'`rep'.dta"
}
}

noi disp "rep `i' for `type_ith' finished ... erasing base cohort "
noi erase "bc_`i'_`type_ith'.dta"


	}
	
postclose `bc_states'
postclose `rep_states'
postclose `bc_estimates'
postclose `rep_estimates'

use "bc_mysimestimates_`type_ith'_sim${nsim_start}_${nsim}", clear
	
	label variable confounding "confounding"
	label variable rep_id "Rep num"
    label variable method "Method"
    label variable theta "θᵢ"
    label variable se "SE(θᵢ)"
	label variable lci "Lower 95% CI"
	label variable uci "Upper 95% CI"
	label variable pval "P value"
    label define methodlab 1 "base_uadj" 2 "base_adj"  3 "cox_uadj" 4 "cox_adj" 5 "cox_match" 6 "p_uadj" 7 "p_adj" 8"p_match"
        label values method methodlab
    sort rep_id method

gen bias = theta - $loghr
gen se2 = se^2
gen cov= 1 if lci <= $loghr & uci >= $loghr
replace cov= 0 if cov==. & lci!=. & uci!=.

sort method rep_id
save "bc_mysimestimates_`type_ith'_sim${nsim_start}_${nsim}", replace
list, noobs sepby(confounding method)

use "rep_mysimestimates_`type_ith'_sim${nsim_start}_${nsim}", clear 
*labels
    label variable rep_id "Rep num"
	label variable confounding  "confounding"
	label variable getallpossible "implementation"
	label variable nomatch "N exposed unmatched"
    label variable method "Method"
    label variable theta "θᵢ"
    label variable se "SE(θᵢ)"
	label variable lci "Lower 95% CI"
	label variable uci "Upper 95% CI"
	label variable pval "P value"
    label define methodlab 1 "base_uadj" 2 "base_adj"  3 "cox_uadj" 4 "cox_adj" 5 "cox_match" 6 "p_uadj" 7 "p_adj" 8"p_match"
        label values method methodlab
    sort getallpossible rep_id ratio method

gen bias = theta - $loghr
gen se2 = se^2
gen cov= 1 if lci <= $loghr & uci >= $loghr
replace cov= 0 if cov==. & lci!=. & uci!=.

sort getallpossible ratio method rep_id

	
save "rep_mysimestimates_`type_ith'_sim${nsim_start}_${nsim}",replace

//list, noobs sepby(getallpossible ratio method)



}

//}
noi disp "SIMULATION COMPLETE"