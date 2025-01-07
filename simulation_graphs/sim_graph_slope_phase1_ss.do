set scheme white_tableau
* Define file path
local resultsDir "$resultsDir//sim_summary_simsummtch1_2000_phase1_ss.dta"

* Load the data
use "`resultsDir'", clear

gen test1 = "emp" if test== "Empirical SE"
replace test1 = "rel_err" if test== "Relative % error in ModSE"

*local test_loop emp rel_err
*local conf_loop 10 sconf

local test_loop emp rel_err
local conf_loop 10 sconf



	* Define LSHTM colors

	local color3 = "30 34 170"     //(blue)
	local color4 = "0 174 199"   //  (light blue)

	local color5 = "13 82 87"   // (green)
	local color6 = "0 191 111"   // (light green)

	local color1 = "254 80 0"      // (orange)
	local color2 = "255 177 187"    // (yellow)

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
	
		// Offset based on getallpossible
	
	label define ratio_lab 1 "ratio 1:1" 3 "ratio 3:1" 5 "ratio 5:1" 10 "ratio 10:1"
	label val ratio ratio_lab

	gen timescale = ""
	replace timescale = "time in study" if inlist(method_reordered, 2, 4, 6, 8, 10, 12)
	replace timescale = "age" if inlist(method_reordered, 1, 3, 5, 7, 9, 11)

foreach t of local test_loop {
	foreach c of local conf_loop {
	
	preserve

	keep if test1 == "`t'"
	local title = test[1]

	if "`c'" == "10" {
		
		keep if strpos(confounding, "10") > 0
		* Encode confounding as a factor variable
		gen newconf = .
		replace newconf = 1 if confounding == "10_nconf"
		replace newconf = 2  if confounding == "10_sconf"
		replace newconf = 3  if confounding == "10_bconf"
		label define conf_lbl 1 "10_nconf" 2 "10_sconf" 3 "10_bconf"
		label values newconf conf_lbl
		tab newconf
		}

	if "`c'" == "sconf" {

		keep if strpos(confounding, "sconf") > 0
		
		gen newconf = .
		replace newconf = 1 if confounding == "1_sconf"
		replace newconf = 2 if confounding == "10_sconf"
		label define conf_lbl 1 "1_sconf" 2 "10_sconf" 
		label values newconf conf_lbl
		tab newconf
		}

	sort newconf ratio timescale 


	* slight offset if variable is time

	// Offset timescales
	gen conf_offset = . // Initialize conf_offset as a copy of newconf
	
	replace conf_offset = 1 if ratio == 1 
	replace conf_offset = 3 if ratio == 3
	replace conf_offset = 5 if ratio == 5
	replace conf_offset = 7 if ratio == 10
	
	replace conf_offset = conf_offset if newconf == 1
	replace conf_offset = conf_offset + 0.75 if newconf == 2
	replace conf_offset = conf_offset + 1.5 if newconf == 3
	

	replace conf_offset = conf_offset if method_reordered==1 | method_reordered ==2
	replace conf_offset = conf_offset + 0.2 if  method_reordered==3 | method_reordered ==4
	replace conf_offset = conf_offset + 0.4  if method_reordered==5 | method_reordered ==6
	replace conf_offset = conf_offset + 0.6  if method_reordered==7 | method_reordered == 8 
	replace conf_offset = conf_offset + 0.8 if method_reordered >= 9
	
	
	replace conf_offset = conf_offset if getallpossible == "replacement"

	label values conf_offset conf_lbl

	* Generate the plot for all methods
	twoway ///
		connected avg conf_offset if inlist(method_reordered, 1,2) & ratio== 1 & getallpossible == "noreplacement", ///
			lpattern(solid) msymbol(triangle) lcolor("`color1'") mcolor("`color1'") /// 
		|| connected avg conf_offset if inlist(method_reordered, 1,2) & ratio== 1 & getallpossible == "replacement", ///
			lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle) lcolor("`color1'") mcolor("`color1'") ///
		|| connected avg conf_offset if inlist(method_reordered, 3,4) & ratio== 1 & getallpossible == "noreplacement", ///
			lpattern(solid) msymbol(triangle) lcolor("`color2'") mcolor("`color2'") ///
		|| connected avg conf_offset if inlist(method_reordered, 3,4) & ratio== 1 & getallpossible == "replacement", ///
		   lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle) lcolor("`color2'") mcolor("`color2'") ///
		|| connected avg conf_offset if inlist(method_reordered, 5,6) & ratio== 1 & getallpossible == "noreplacement", ///
			lpattern(solid) msymbol(triangle) lcolor("`color3'") mcolor("`color3'") ///
		|| connected avg conf_offset if inlist(method_reordered, 5,6) & ratio== 1 & getallpossible == "replacement", ///
			lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle) lcolor("`color3'") mcolor("`color3'") ///
		|| connected avg conf_offset if inlist(method_reordered, 7,8) & ratio== 1 & getallpossible == "noreplacement", ///
			lpattern(solid) msymbol(triangle) lcolor("`color4'") mcolor("`color4'") ///
		|| connected avg conf_offset if inlist(method_reordered, 7,8) & ratio== 1 & getallpossible == "replacement", ///
			lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle) lcolor("`color4'") mcolor("`color4'") ///
		|| connected avg conf_offset if inlist(method_reordered, 9,10) & ratio== 1 & getallpossible == "noreplacement", ///
			lpattern(solid) msymbol(triangle) lcolor("`color5'") mcolor("`color5'") ///
		|| connected avg conf_offset if inlist(method_reordered, 9,10) & ratio== 1 & getallpossible == "replacement", ///
			lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle) lcolor("`color5'") mcolor("`color5'") ///	
		|| connected avg conf_offset if inlist(method_reordered, 11,12) & ratio== 1 & getallpossible == "noreplacement", ///
			lpattern(solid) msymbol(triangle) lcolor("`color6'") mcolor("`color6'") ///
		|| connected avg conf_offset if inlist(method_reordered, 11,12) & ratio== 1 & getallpossible == "replacement", ///
			lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle)  lcolor("`color6'") mcolor("`color6'") ///
			///
		|| connected avg conf_offset if inlist(method_reordered, 1,2) & ratio== 3 & getallpossible == "noreplacement", ///
			lpattern(solid) msymbol(triangle) lcolor("`color1'") mcolor("`color1'") /// 
		|| connected avg conf_offset if inlist(method_reordered, 1,2) & ratio== 3 & getallpossible == "replacement", ///
			lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle) lcolor("`color1'") mcolor("`color1'") ///
		|| connected avg conf_offset if inlist(method_reordered, 3,4) & ratio== 3 & getallpossible == "noreplacement", ///
			lpattern(solid) msymbol(triangle) lcolor("`color2'") mcolor("`color2'") ///
		|| connected avg conf_offset if inlist(method_reordered, 3,4) & ratio== 3 & getallpossible == "replacement", ///
		   lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle) lcolor("`color2'") mcolor("`color2'") ///
		|| connected avg conf_offset if inlist(method_reordered, 5,6) & ratio== 3 & getallpossible == "noreplacement", ///
			lpattern(solid) msymbol(triangle) lcolor("`color3'") mcolor("`color3'") ///
		|| connected avg conf_offset if inlist(method_reordered, 5,6) & ratio== 3 & getallpossible == "replacement", ///
			lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle) lcolor("`color3'") mcolor("`color3'") ///
		|| connected avg conf_offset if inlist(method_reordered, 7,8) & ratio== 3 & getallpossible == "noreplacement", ///
			lpattern(solid) msymbol(triangle) lcolor("`color4'") mcolor("`color4'") ///
		|| connected avg conf_offset if inlist(method_reordered, 7,8) & ratio== 3 & getallpossible == "replacement", ///
			lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle) lcolor("`color4'") mcolor("`color4'") ///
		|| connected avg conf_offset if inlist(method_reordered, 9,10) & ratio== 3 & getallpossible == "noreplacement", ///
			lpattern(solid) msymbol(triangle) lcolor("`color5'") mcolor("`color5'") ///
		|| connected avg conf_offset if inlist(method_reordered, 9,10) & ratio== 3 & getallpossible == "replacement", ///
			lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle) lcolor("`color5'") mcolor("`color5'") ///	
		|| connected avg conf_offset if inlist(method_reordered, 11,12) & ratio== 3 & getallpossible == "noreplacement", ///
			lpattern(solid) msymbol(triangle)  lcolor("`color6'") mcolor("`color6'") ///
		|| connected avg conf_offset if inlist(method_reordered, 11,12) & ratio== 3 & getallpossible == "replacement", ///
			lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle)  lcolor("`color6'") mcolor("`color6'") ///
			///
		|| connected avg conf_offset if inlist(method_reordered, 1,2) & ratio== 5 & getallpossible == "noreplacement", ///
			lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle) lcolor("`color1'") mcolor("`color1'") ///
		|| connected avg conf_offset if inlist(method_reordered, 1,2) & ratio== 5 & getallpossible == "replacement", ///
			lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle) lcolor("`color1'") mcolor("`color1'") ///
		|| connected avg conf_offset if inlist(method_reordered, 3,4) & ratio== 5 & getallpossible == "noreplacement", ///
			lpattern(solid) msymbol(triangle) lcolor("`color2'") mcolor("`color2'") ///
		|| connected avg conf_offset if inlist(method_reordered, 3,4) & ratio== 5 & getallpossible == "replacement", ///
		   lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle) lcolor("`color2'") mcolor("`color2'") ///
		|| connected avg conf_offset if inlist(method_reordered, 5,6) & ratio== 5 & getallpossible == "noreplacement", ///
			lpattern(solid) msymbol(triangle) lcolor("`color3'") mcolor("`color3'") ///
		|| connected avg conf_offset if inlist(method_reordered, 5,6) & ratio== 5 & getallpossible == "replacement", ///
			lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle) lcolor("`color3'") mcolor("`color3'") ///
		|| connected avg conf_offset if inlist(method_reordered, 7,8) & ratio== 5 & getallpossible == "noreplacement", ///
			lpattern(solid) msymbol(triangle) lcolor("`color4'") mcolor("`color4'") ///
		|| connected avg conf_offset if inlist(method_reordered, 7,8) & ratio== 5 & getallpossible == "replacement", ///
			lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle) lcolor("`color4'") mcolor("`color4'") ///
		|| connected avg conf_offset if inlist(method_reordered, 9,10) & ratio== 5 & getallpossible == "noreplacement", ///
			lpattern(solid) msymbol(triangle) lcolor("`color5'") mcolor("`color5'") ///
		|| connected avg conf_offset if inlist(method_reordered, 9,10) & ratio== 5 & getallpossible == "replacement", ///
			lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle) lcolor("`color5'") mcolor("`color5'") ///	
		|| connected avg conf_offset if inlist(method_reordered, 11,12) & ratio== 5 & getallpossible == "noreplacement", ///
			lpattern(solid) msymbol(triangle)  lcolor("`color6'") mcolor("`color6'") ///
		|| connected avg conf_offset if inlist(method_reordered, 11,12) & ratio== 5 & getallpossible == "replacement", ///
			lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle)  lcolor("`color6'") mcolor("`color6'") ///
			///
		|| connected avg conf_offset if inlist(method_reordered, 1,2) & ratio== 10 & getallpossible == "noreplacement", ///
			lpattern(solid) msymbol(triangle) lcolor("`color1'") mcolor("`color1'") /// 
		|| connected avg conf_offset if inlist(method_reordered, 1,2) & ratio== 10 & getallpossible == "replacement", ///
			lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle) lcolor("`color1'") mcolor("`color1'") ///
		|| connected avg conf_offset if inlist(method_reordered, 3,4) & ratio== 10 & getallpossible == "noreplacement", ///
			lpattern(solid) msymbol(triangle) lcolor("`color2'") mcolor("`color2'") ///
		|| connected avg conf_offset if inlist(method_reordered, 3,4) & ratio== 10 & getallpossible == "replacement", ///
		   lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle) lcolor("`color2'") mcolor("`color2'") ///
		|| connected avg conf_offset if inlist(method_reordered, 5,6) & ratio== 10 & getallpossible == "noreplacement", ///
			lpattern(solid) msymbol(triangle) lcolor("`color3'") mcolor("`color3'") ///
		|| connected avg conf_offset if inlist(method_reordered, 5,6) & ratio== 10 & getallpossible == "replacement", ///
			lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle) lcolor("`color3'") mcolor("`color3'") ///
		|| connected avg conf_offset if inlist(method_reordered, 7,8) & ratio== 10 & getallpossible == "noreplacement", ///
			lpattern(solid) msymbol(triangle) lcolor("`color4'") mcolor("`color4'") ///
		|| connected avg conf_offset if inlist(method_reordered, 7,8) & ratio== 10 & getallpossible == "replacement", ///
			lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle) lcolor("`color4'") mcolor("`color4'") ///
		|| connected avg conf_offset if inlist(method_reordered, 9,10) & ratio== 10 & getallpossible == "noreplacement", ///
			lpattern(solid) msymbol(triangle) lcolor("`color5'") mcolor("`color5'") ///
		|| connected avg conf_offset if inlist(method_reordered, 9,10) & ratio== 10 & getallpossible == "replacement", ///
			lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle) lcolor("`color5'") mcolor("`color5'") ///	
		|| connected avg conf_offset if inlist(method_reordered, 11,12) & ratio== 10 & getallpossible == "noreplacement", ///
			lpattern(solid) msymbol(triangle)  lcolor("`color6'") mcolor("`color6'") ///
		|| connected avg conf_offset if inlist(method_reordered, 11,12) & ratio== 10 & getallpossible == "replacement", ///
			lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle)  lcolor("`color6'") mcolor("`color6'") ///
			///
		|| connected avg conf_offset if inlist(method_reordered, 1,2) & ratio== 10 & getallpossible == "noreplacement", ///
			lpattern(solid) msymbol(triangle) lcolor("`color1'") mcolor("`color1'") /// 
		|| connected avg conf_offset if inlist(method_reordered, 1,2) & ratio== 10 & getallpossible == "replacement", ///
			lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle) lcolor("`color1'") mcolor("`color1'") ///
		|| connected avg conf_offset if inlist(method_reordered, 3,4) & ratio== 10 & getallpossible == "noreplacement", ///
			lpattern(solid) msymbol(triangle) lcolor("`color2'") mcolor("`color2'") ///
		|| connected avg conf_offset if inlist(method_reordered, 3,4) & ratio== 10 & getallpossible == "replacement", ///
		   lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle) lcolor("`color2'") mcolor("`color2'") ///
		|| connected avg conf_offset if inlist(method_reordered, 5,6) & ratio== 10 & getallpossible == "noreplacement", ///
			lpattern(solid) msymbol(triangle) lcolor("`color3'") mcolor("`color3'") ///
		|| connected avg conf_offset if inlist(method_reordered, 5,6) & ratio== 10 & getallpossible == "replacement", ///
			lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle) lcolor("`color3'") mcolor("`color3'") ///
		|| connected avg conf_offset if inlist(method_reordered, 7,8) & ratio== 10 & getallpossible == "noreplacement", ///
			lpattern(solid) msymbol(triangle) lcolor("`color4'") mcolor("`color4'") ///
		|| connected avg conf_offset if inlist(method_reordered, 7,8) & ratio== 10 & getallpossible == "replacement", ///
			lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle) lcolor("`color4'") mcolor("`color4'") ///
		|| connected avg conf_offset if inlist(method_reordered, 9,10) & ratio== 10 & getallpossible == "noreplacement", ///
			lpattern(solid) msymbol(triangle) lcolor("`color5'") mcolor("`color5'") ///
		|| connected avg conf_offset if inlist(method_reordered, 9,10) & ratio== 10 & getallpossible == "replacement", ///
			lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle) lcolor("`color5'") mcolor("`color5'") ///	
		|| connected avg conf_offset if inlist(method_reordered, 11,12) & ratio== 10 & getallpossible == "noreplacement", ///
			lpattern(solid) msymbol(triangle)  lcolor("`color6'") mcolor("`color6'") ///
		|| connected avg conf_offset if inlist(method_reordered, 11,12) & ratio== 10 & getallpossible == "replacement", ///
			lpattern(shortdash) mfcolor(none) mlwidth(medium) msymbol(triangle)  lcolor("`color6'") mcolor("`color6'") ///
			by(timescale getallpossible, col(2) legend(off) title("Change in `title'", size(medium)) subtitle("from nconf>sconf>bconf") note("") imargin(0 0 0 0 )) /// Panel by `getallpossible` and `ratio`
		xlabel(1 "1:1" 3 "3:1" 5 "5:1" 7 "10:1" 9 " ", labsize(medium) angle(45)) /// Confounding labels
		xtitle("") ///
		yline(0) ///
		ylabel(, format(%9.2g) labsize(medium)) /// Y-axis customization
		ytitle("") ///
		legend(order(1 "Cox-UA rep"  2 "Cox-UA norep" ///
		3 "Pois-UA rep"  4 "Pois-UA norep" ///
		5 "Cox-A rep"  6 "Cox-A norep" ///
		7 "Pois-A  rep"  8 "Pois-A  norep" /// 
		9 "Cox Strat rep"  10 "Cox Strat norep" /// 
		11 "Cond Pois  rep"  12 "Cond Pois norep") cols(2) ///
			   region(lcolor(white)) position(11)) name(slope_`t'_`c', replace) 

		graph export "slope_`t'_`c'.emf", replace
			  drop newconf
			   label drop conf_lbl
			   drop conf_offset
   restore				   
	}

		}
		
