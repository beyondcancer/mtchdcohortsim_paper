
//quietly{
foreach l in $conf_test {
	local exp_ithbi = word("$conf_expbi", `l')
	local out_ithbi = word("$conf_outbi", `l')
	local type_ith = word("$conf_type", `l')
	local lambda_ith = word("$lambda_exp", `l')
	local gamma_ith = word("$gamma_exp", `l')
	local lambda_o_ith = word("$lambda_out", `l')
	local gamma_o_ith = word("$gamma_out", `l')
	local exp_ithm = word("$conf_expmult", `l')
	local out_ithm = word("$conf_outmult", `l')
	
noi disp "Initiating the `type_ith' simulation loop, this is loop `l'"

 tempname bc_estimates rep_estimates	
 postfile `bc_estimates' float(nobs) str20 confounding int(rep_id) int(ratio) int(method) float(theta se lci uci pval) using bc_mysimestimates_`type_ith'_sim${nsim_start}_${nsim}.dta, replace	
postfile `rep_estimates' str20 confounding str20 getallpossible int(rep_id) int(ratio) int(method) float(nomatch theta se lci uci pval) using rep_mysimestimates_`type_ith'_sim${nsim_start}_${nsim}.dta, replace 

tempname bc_states rep_states	
postfile `bc_states' str10 confounding int(rep_id) str2000 bc1 str2000 bc2 str1000 bc3 using bc_mysimstates_`type_ith'_sim${nsim_start}_${nsim}, replace
postfile `rep_states' str10 confounding int(rep_id) str20 implementation int(ratio) str2000 state1 str2000 state2 str1000 state3 using rep_mysimstates_`type_ith'_sim${nsim_start}_${nsim}.dta, replace ///figure this out later

tempname cohort_deets
postfile `cohort_deets' str20 confounding str20 getallpossible int(rep_id) int(ratio) int(n_set) int(median_set) int(set_25) int(set_75) float(full_sets) float(time_mtch) using rep_cohortdeets_`type_ith'_sim${nsim_start}_${nsim}.dta, replace

forval i= $nsim_start/$nsim {	

noi disp "this is rep `i' out of " ($nsim - $nsim_start + 1)

create_bc, nobs($n) loghr($loghr) i(`i') lambda_exp(`lambda_ith') gamma_exp(`gamma_ith') lambda_out(`lambda_o_ith') gamma_out(`gamma_o_ith')gender_exp(`exp_ithbi') gender_out(`out_ithbi') pracid_exp(`exp_ithm') pracid_out(`out_ithm') conf_type("`type_ith'")


post `bc_states' ("`type_ith'") (`i') (substr(r(bc_state),1,2000)) (substr(r(bc_state),2001, 2000)) (substr(r(bc_state),4001,.)) 

post `bc_estimates' (r(total_n)) ("`type_ith'") (`i') (0) (1) (r(tbase_uadj_est)) (r(tbase_uadj_se)) (r(tbase_uadj_lci)) (r(tbase_uadj_uci)) (r(tbase_uadj_pval))
post `bc_estimates' (r(total_n)) ("`type_ith'") (`i') (0) (2) (r(tbase_adj_est)) (r(tbase_adj_se)) (r(tbase_adj_lci)) (r(tbase_adj_uci)) (r(tbase_adj_pval))
post `bc_estimates' (r(total_n)) ("`type_ith'") (`i') (0) (3) (r(abase_uadj_est)) (r(abase_uadj_se)) (r(abase_uadj_lci)) (r(abase_uadj_uci)) (r(abase_uadj_pval))
post `bc_estimates' (r(total_n)) ("`type_ith'") (`i') (0) (4) (r(abase_adj_est)) (r(abase_adj_se)) (r(abase_adj_lci)) (r(abase_adj_uci)) (r(abase_adj_pval))

noi disp "base cohort created : rep `i' `type_ith' lambda_exp is `lambda_ith'"	

///Info on the cohort
		count if tag 
		qui local total_n = r(N)
		noi disp "Total patients: " `total_n' 
		tab exposed outcome, row
		qui tabulate exposed, matcell(freq_matrix) matrow(var_values)
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

		post `rep_estimates' ("`type_ith'") ("`rep'") (`i') (`j') (5) (r(nomatch)) (r(cox_uadj_a_est)) (r(cox_uadj_a_se))  (r(cox_uadj_a_lci)) (r(cox_uadj_a_uci))  (r(cox_uadj_a_pval))
		post `rep_estimates' ("`type_ith'") ("`rep'") (`i') (`j') (6) (r(nomatch)) (r(cox_adj_a_est))  (r(cox_adj_a_se))  (r(cox_adj_a_lci)) (r(cox_adj_a_uci))  (r(cox_adj_a_pval))
		post `rep_estimates' ("`type_ith'") ("`rep'") (`i') (`j') (7) (r(nomatch)) (r(cox_mtch_a_est)) (r(cox_mtch_a_se))  (r(cox_mtch_a_lci)) (r(cox_mtch_a_uci))  (r(cox_mtch_a_pval))
		post `rep_estimates' ("`type_ith'") ("`rep'") (`i') (`j') (8) (r(nomatch)) (r(p_uadj_a_est))  (r(p_uadj_a_se))  (r(p_uadj_a_lci)) (r(p_uadj_a_uci))  (r(p_uadj_a_pval))
		post `rep_estimates' ("`type_ith'") ("`rep'") (`i') (`j') (9) (r(nomatch)) (r(p_adj_a_est))  (r(p_adj_a_se))  (r(p_adj_a_lci)) (r(p_adj_a_uci))  (r(p_adj_a_pval))
		post `rep_estimates' ("`type_ith'") ("`rep'") (`i') (`j') (10) (r(nomatch)) (r(p_mtch_a_est))  (r(p_mtch_a_se))  (r(p_mtch_a_lci)) (r(p_mtch_a_uci))  (r(p_mtch_a_pval))
		post `rep_estimates' ("`type_ith'") ("`rep'") (`i') (`j') (11) (r(nomatch)) (r(cox_uadj_t_est)) (r(cox_uadj_t_se))  (r(cox_uadj_t_lci)) (r(cox_uadj_t_uci))  (r(cox_uadj_t_pval))
		post `rep_estimates' ("`type_ith'") ("`rep'") (`i') (`j') (12) (r(nomatch)) (r(cox_adj_t_est))  (r(cox_adj_t_se))  (r(cox_adj_t_lci)) (r(cox_adj_t_uci))  (r(cox_adj_t_pval))
		post `rep_estimates' ("`type_ith'") ("`rep'") (`i') (`j') (13) (r(nomatch)) (r(cox_mtch_t_est)) (r(cox_mtch_t_se))  (r(cox_mtch_t_lci)) (r(cox_mtch_t_uci))  (r(cox_mtch_t_pval))
		post `rep_estimates' ("`type_ith'") ("`rep'") (`i') (`j') (14) (r(nomatch)) (r(p_uadj_t_est))  (r(p_uadj_t_se))  (r(p_uadj_t_lci)) (r(p_uadj_t_uci))  (r(p_uadj_t_pval))
		post `rep_estimates' ("`type_ith'") ("`rep'") (`i') (`j') (15) (r(nomatch)) (r(p_adj_t_est))  (r(p_adj_t_se))  (r(p_adj_t_lci)) (r(p_adj_t_uci))  (r(p_adj_t_pval))
		post `rep_estimates' ("`type_ith'") ("`rep'") (`i') (`j') (16) (r(nomatch)) (r(p_mtch_t_est))  (r(p_mtch_t_se))  (r(p_mtch_t_lci)) (r(p_mtch_t_uci))  (r(p_mtch_t_pval))
	
noi return list


post `cohort_deets' ("`type_ith'") ("`rep'") (`i') (`j') (.) (.) (.) (.) (.) (r(mtch_time))

use "temp_`i'`type_ith'`rep'.dta", clear

//Number of controls per set
* Count the number of unexposed individuals in each setid
bysort setid: gen unexposed_count = sum(exposed == 0)
bysort setid: gen setid_count = _N


* Keep only the last observation per setid (representing the total count of unexposed in each setid)
bysort setid: keep if _n==_N
local total = _N
disp `total'

tab setid_count, matcell(freq) matrow(percent)  // Perform the tabulation and create matrices
local rows = rowsof(freq)  // Get the total number of rows in the matrix
// Save the total from the last row of the frequency matrix
local total_count = freq[`rows', 1]

local full_sets = (`total_count'/`total')*100 

display "Total: " `full_sets'

* Calculate the mean, 25th percentile, and 75th percentile of unexposed_count across all sets
qui summarize unexposed_count, detail

return list

post `cohort_deets' ("`type_ith'") ("`rep'") (`i') (`j') (r(N)) (r(p50)) (r(p25)) (r(p75)) (`full_sets') (r(mtch_time))

* Display results
di "Mean: " r(mean)
di "Median: " r(p50)
di "25th percentile (p25): " r(p25)
di "75th percentile (p75): " r(p75)

noi disp "erasing `rep' `type_ith' matched cohort data"
noi erase "getmatchedcohort`i'`type_ith'`rep'.dta"
noi erase "temp_`i'`type_ith'`rep'.dta"


}
}

noi disp "rep `i' for `type_ith' finished ... erasing base cohort "
noi erase "bc_`i'_`type_ith'.dta"


	}
	
postclose `bc_states'
postclose `rep_states'
postclose `bc_estimates'
postclose `rep_estimates'
postclose `cohort_deets'

use "bc_mysimestimates_`type_ith'_sim${nsim_start}_${nsim}", clear
	
	label variable confounding "confounding"
	label variable rep_id "Rep num"
    label variable method "Method"
    label variable theta "θᵢ"
    label variable se "SE(θᵢ)"
	label variable lci "Lower 95% CI"
	label variable uci "Upper 95% CI"
	label variable pval "P value"
    label define methodlab 1 "tbase_uadj" 2 "tbase_adj" 3 "abase_uadj" 4 "abase_adj" ///
							5 "cox_uadj_a" 6 "cox_adj_a" 7 "cox_match_a" 8 "p_uadj_a" 9 "p_adj_a" 10 "p_match_a" ///
							11 "cox_uadj_t" 12 "cox_adj_t" 13 "cox_match_t" 14 "p_uadj_t" 15 "p_adj_t" 16 "p_match_t"
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
    label define methodlab 1 "tbase_uadj" 2 "tbase_adj" 3 "abase_uadj" 4 "abase_adj" ///
							5 "cox_uadj_a" 6 "cox_adj_a" 7 "cox_match_a" 8 "p_uadj_a" 9 "p_adj_a" 10 "p_match_a" ///
							11 "cox_uadj_t" 12 "cox_adj_t" 13 "cox_match_t" 14 "p_uadj_t" 15 "p_adj_t" 16 "p_match_t"
    label values method methodlab
    sort getallpossible rep_id ratio method

gen bias = theta - $loghr
gen se2 = se^2
gen cov= 1 if lci <= $loghr & uci >= $loghr
replace cov= 0 if cov==. & lci!=. & uci!=.

sort getallpossible ratio method rep_id
list, noobs sepby(confounding method)

	
save "rep_mysimestimates_`type_ith'_sim${nsim_start}_${nsim}",replace


}

noi disp "SIMULATION COMPLETE"