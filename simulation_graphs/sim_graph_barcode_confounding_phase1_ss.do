set scheme white_tableau

// Load and prepare the data
use "$resultsDir//sim_summary_simsummtch1_2000_phase1_ss.dta", clear

keep if test == "Bias" | test == "Coverage of 95% CIs"
keep if strpos(confounding, "10_") > 0

// Create the bias indicator
gen indicator = (MCSE_hi >= 0.0001 & MCSE_low <= 0.0001) | (MCSE_hi >= 95 & MCSE_low <= 95)

// Assign numeric labels to ratio_lab with bigger gaps between ratios
gen ratio_lab = .
replace ratio_lab = 3 if ratio == 1
replace ratio_lab = 5 if ratio == 3  // Gap between 1:1 and 3:1
replace ratio_lab = 7 if ratio == 5 // Gap between 3:1 and 5:1
replace ratio_lab = 9 if ratio == 10 // Gap between 5:1 and 10:1

// Adjust offsets for better spacing and add gaps
gen ratio_offset = .
replace ratio_offset = ratio_lab - 0.5 if confounding == "10_nconf" & getallpossible == "noreplacement"
replace ratio_offset = ratio_lab - 0.3 if confounding == "10_nconf" & getallpossible == "replacement"
replace ratio_offset = ratio_lab - 0.1 if confounding == "10_sconf" & getallpossible == "noreplacement"
replace ratio_offset = ratio_lab + 0.1 if confounding == "10_sconf" & getallpossible == "replacement"
replace ratio_offset = ratio_lab + 0.3 if confounding == "10_bconf" & getallpossible == "noreplacement"
replace ratio_offset = ratio_lab + 0.5 if confounding == "10_bconf" & getallpossible == "replacement"

// Reorder methods for plotting
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

// Generate the single graph with confounding level labels and larger gaps between ratios
twoway ///
    scatter method_reordered ratio_offset if indicator == 1 & getallpossible == "noreplacement", ///
        msymbol(square) mcolor(maroon) msize(medium)  /// Without replacement
    || scatter method_reordered ratio_offset  if indicator == 1 & getallpossible == "replacement", ///
        msymbol(square) mcolor(gold) msize(medium)  /// With replacement
    , ylabel(1 "Cox-UA (Age)" 2 "Cox-UA (Time)" 3 "Pois-UA (Age)" 4 "Pois-UA (Time)" ///
             5 "Cox-A (Age)" 6 "Cox-A (Time)" 7 "Pois-A (Age)" 8 "Pois-A (Time)" ///
             9 "Cox-Strat (Age)" 10 "Cox-Strat (Time)" 11 "Pois-Cond (Age)" 12 "Pois-Cond (Time)", ///
             labsize(small) angle(-20)) ///
      xlabel(2.6 "nconf" 4.6 "nconf" 6.6 "nconf" 8.6 "nconf" ///
             3 "1:1 sconf" 5 "3:1 sconf" 7 "5:1 sconf" 9 "10:1 sconf" ///
             3.4 "bconf" 5.4 "bconf" 7.4 "bconf" 9.4 "bconf", ///
             labsize(small) angle(90)) ///
			 legend(off) ///
			 title("Bias = 0 Coverage = 95%", size(medium)) ///
			 xtitle("") ytitle("") name(barcode_confounding, replace)	 
   *   legend(order(1 "Without replacement" 2 "With replacement") col(1) size(small)) ///


// Save the graph
graph export "$resultsDir/barcode_confounding.emf", replace
