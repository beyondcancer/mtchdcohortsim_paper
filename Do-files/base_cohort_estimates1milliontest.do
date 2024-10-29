
//I also think we can empirically estimate the matched estimands by generating eg a million exposed people, duplicating them and jittering the age of one (who becomes unexposed), and then generating outcomes assuming one is exposed and the other not.  We delete pairs where the exposed is simulated an outcome time prior to exposure and regenerate outcomes for unexposed that are before. I think that works. I'll try and see. It will be reassuring if these estimates then correspond to what we've been using

set seed 88

clear
set obs 100000


//generate a patid
 gen patid= _n

//generate age (based on UK population)
 gen age_start = rnormal(40, 12)
*drop if age_start <18 // KA note: tmporarily removed to check if removing less people works better


//generate gender. Recoded 1 & 2 to fit getmatched cohort program,
 gen gender = runiform()<=0.5
 replace gender = 2 if gender==0
label define sexlab 1 "male" 2 "female"
label values gender sexlab
recast byte gender 

//generate gp variable & group from 0 - 13 practices (equivalent to CPRD if size were 100000)
gen gengp = runiform(0, 13)
gen pracid = round(gengp)
replace pracid= 13 if pracid==0
drop gengp
tab pracid 

//gen min date and max-date of entering and leaving database
gen earliest_date = date("01apr1997", "DMY")
gen latest_date = date("30nov2018", "DMY")
format earliest_date %td
format latest_date %td

//generate random startdates within this period
//gen startdate= round(earliest_date + (latest_date - earliest_date) * runiform()) - ROUNDING REMOVED
gen startdate= earliest_date + (latest_date - earliest_date) * runiform()
format startdate %td

*gen yob
cap drop dob
gen dob = startdate - (age_start * 365.25)
format dob %td
gen yob= year(dob)


//Calculate total time between startdate and latest_date (in years)
gen total_time = (latest_date-startdate)/365.25
gen total_time_days = latest_date-startdate


//generate exposure variable 100000
gen exposed = 1

// Step 2: Duplicate each individual to create a matched unexposed pair
expand 2
gen setid = "10" + string(patid)
bysort setid: replace exposed = 0 if _n== _N

// Step 3: Jitter age for unexposed individuals (by 3 yrs)
replace age_start = age_start + rnormal(0, 3) if exposed == 0

// Generate a random integer between 0 and total_time
gen itime = round(runiform() * total_time_days)

gen indexdate = startdate + itime
format indexdate %td

replace indexdate=. if exposed==0
//assert startdate <= indexdate & indexdate <= latest_date if indexdate!=. 10 contradictions

* Tidy dataset
order startdate, before(latest_date)

//generate outcome data. maximum time is now stipulated by t_years (total time per row)
survsim stime outcome, dist(weib) lambda(1.19) gamma(0.012) cov(exposed 0 age_start 0 gender 0 pracid 0) maxtime(total_time_days) //ltruncated(t0_years) maxtime(t_years)


//generate the date the outcome occurs
gen outcomedate = (startdate + stime*365.25)
format outcomedate %td 
replace outcomedate=. if outcome==0


//step 3: delete pairs where outcome happens before exposure

bysort setid: drop if exposed ==1 & outcomedate < indexdate //drop if exposure is before outcome
bysort setid: drop if _N < 2 //drop matched pair


//step 4: regenerate outcomes for unexposed that are before
	// there is no indexdate for unexposed? so no need to regenerate?


// step 5: generate enddate
egen enddate = rowmin(latest_date outcomedate)

order patid dob yob age_start gender pracid ///
	startdate indexdate	///
	outcomedate  latest_date ///
	exposed outcome enddate 

	
	qui stset enddate, origin(startdate) enter(indexdate) fail(outcome) /// 
			id(patid) scale(365.25)
			
			stcox i.exposed, nohr
			stcox i.exposed age_start gender pracid, nohr
		
			
						capture { 
			qui stcox i.exposed, nohr
			return scalar base_uadj_est= r(table)[1,2] 
			return scalar base_uadj_se= r(table)[2,2]
			return scalar base_uadj_pval= r(table)[4,2]
			return scalar base_uadj_lci= r(table)[5,2]
			return scalar base_uadj_uci= r(table)[6,2]
			}
			
			if _rc!=0 {
				noi disp "No convergence achieved for base Cox unadjusted model in rep `i'"
		    return scalar base_uadj_est= .
			return scalar base_uadj_se= .
			return scalar base_uadj_pval= .
			return scalar base_uadj_lci= .
			return scalar base_uadj_uci= .
			}
			
			capture{
			qui stcox i.exposed age_start gender pracid, nohr
		
			return scalar base_adj_est= r(table)[1,2] 
			return scalar base_adj_se= r(table)[2,2] 
			return scalar base_adj_pval= r(table)[4,2]
			return scalar base_adj_lci= r(table)[5,2]
			return scalar base_adj_uci= r(table)[6,2]
			}
			
			if _rc!=0 {
				noi disp "No convergence achieved for base Cox adjusted model in rep `i'"
		    return scalar base_adj_est= .
			return scalar base_adj_se= .
			return scalar base_adj_pval= .
			return scalar base_adj_lci= .
			return scalar base_adj_uci= .
			}
			
			drop _*
		
	    bysort patid: gen tag = _n == 1
		count if tag 
		return scalar total_n = r(N)
		
		save "bc_`i'_`conf_type'", replace
		

		

				
///////

program define matched_sim, rclass
version 17

syntax [, getallpossible(string) ratio(int 5) i(int 1) conf_type(string)]
             
			qui use "bc_`i'_`conf_type'", clear

			if "`getallpossible'" == "noreplacement"{	
				
			qui return local state_`getallpossible'_`ratio' = c(rngstate)

			
			 getmatchedcohort, cprddb(gold) practice gender yob yobwindow(1) ctrlsperexp(`ratio') savedir($workdir) filesuffix(`i'`conf_type'`getallpossible') followup dontcheck
			
			qui use "getmatchedcohort`i'`conf_type'`getallpossible'"
			qui return scalar nomatch= r(nomtch)	
			qui merge 1:1 patid exposed using "bc_`i'_`conf_type'" 
			
			qui keep if _merge==3 
			
			qui gen age_index= (indexdate-dob)/365.25
			}
				
			if "`getallpossible'" == "replacement"{
				
			qui return local state_`getallpossible'_`ratio' = c(rngstate) 
			
			 getmatchedcohort, cprddb(gold) practice gender yob yobwindow(1) savedir($workdir) getallpossible filesuffix(`i'`conf_type'`getallpossible') followup dontcheck
				
			qui use "getmatchedcohort`i'`conf_type'`getallpossible'"
			qui return scalar nomatch= r(nomtch)
			qui merge m:1 patid exposed using "bc_`i'_`conf_type'" 
			
			qui keep if _merge==3 
							
				gen count=`ratio'+1
				gen rand = runiform()
				gsort setid -exposed  rand
				by setid, sort: keep if exposed==1| _n<=count
			    gen og_patid= patid
				replace patid = _n
			
			 gen age_index= (indexdate-dob)/365.25
			}
				
			*Step 3: Calculate the hr for each type of model// adjustment


			*Cox regressions stset data*
			qui stset enddate, origin(indexdate) enter(indexdate) fail(outcome) id(patid) scale(365.25)

			*unadjusted
			capture { 
			qui stcox i.exposed, nohr
			return scalar cox_uadj_est= r(table)[1,2] 
			return scalar cox_uadj_se= r(table)[2,2]
			return scalar cox_uadj_pval= r(table)[4,2]
			return scalar cox_uadj_lci= r(table)[5,2]
			return scalar cox_uadj_uci= r(table)[6,2]
			}
			
			if _rc!=0 {
				noi disp "No convergence achieved for Cox unadjusted model in rep `i'"
		    return scalar cox_uadj_est= .
			return scalar cox_uadj_se= .
			return scalar cox_uadj_pval= .
			return scalar cox_uadj_lci= .
			return scalar cox_uadj_uci= .
			}
			
			*adjusted for matching variables
			capture { 
			qui stcox i.exposed age_index gender pracid, nohr 
			return scalar cox_adj_est= r(table)[1,2] 
			return scalar cox_adj_se= r(table)[2,2]
			return scalar cox_adj_pval= r(table)[4,2]
			return scalar cox_adj_lci= r(table)[5,2]
			return scalar cox_adj_uci= r(table)[6,2]
			
			}
						if _rc!=0 {
				noi disp "No convergence achieved for Cox adjusted model in rep `i'"
			return scalar cox_adj_est= .
			return scalar cox_adj_se= .
			return scalar cox_adj_pval= .
			return scalar cox_adj_lci= .
			return scalar cox_adj_uci= .
			}
			
			*adjusted by matched set
			capture { 
			qui stcox i.exposed, strata(setid) nohr
			return scalar cox_mtch_est= r(table)[1,2]
			return scalar cox_mtch_se= r(table)[2,2]
			return scalar cox_mtch_pval= r(table)[4,2]
			return scalar cox_mtch_lci= r(table)[5,2]
			return scalar cox_mtch_uci= r(table)[6,2]
			}
			
						if _rc!=0 {
				noi disp "No convergence achieved for Cox adjusted model in rep `i'"
				return scalar cox_mtch_est= .
				return scalar cox_mtch_se= .
				return scalar cox_mtch_pval= .
				return scalar cox_mtch_lci= .
				return scalar cox_mtch_uci= .
				
			}
				
//split data for Lexis expansion
			//gen follow up time variable
	       stsplit survival, at(1 2 3 4 5)
	       gen y= _t - _t0
		
			*Poisson regressions
				*unadjusted
			capture {
			qui streg  i.exposed i.survival, dist(exp) nohr
			return scalar p_uadj_est= r(table)[1,2]
			return scalar p_uadj_se= r(table)[2,2]
			return scalar p_uadj_pval= r(table)[4,2]
			return scalar p_uadj_lci= r(table)[5,2]
			return scalar p_uadj_uci= r(table)[6,2]
			}
			
						if _rc!=0 {
				noi disp "No convergence achieved for Poisson unadjusted model in rep `i'"
			return scalar p_uadj_est= .
			return scalar p_uadj_se= .
			return scalar p_uadj_pval= .
			return scalar p_uadj_lci= .
			return scalar p_uadj_uci= .
			}
			
			*adjusted for matching variables 
			capture {
			 qui streg  i.exposed age_index gender pracid i.survival, dist(exp) nohr
			return scalar p_adj_est= r(table)[1,2]
			return scalar p_adj_se= r(table)[2,2]
			return scalar p_adj_pval= r(table)[4,2]
			return scalar p_adj_lci= r(table)[5,2]
			return scalar p_adj_uci= r(table)[6,2]
			}
									if _rc!=0 {
				noi disp "No convergence achieved for Poisson adjusted model in rep `i'"
				return scalar p_adj_est= .
				return scalar p_adj_se= .
				return scalar p_adj_pval= .
				return scalar p_adj_lci= .
				return scalar p_adj_uci= .
			}
			
			*conditional poisson*
			
			capture { 
	
			qui xtpois _d i.exposed i.survival, fe i(setid) exp(y)
			return scalar p_mtch_est= r(table)[1,2]
			return scalar p_mtch_se= r(table)[2,2]
			return scalar p_mtch_pval= r(table)[4,2]
			return scalar p_mtch_lci= r(table)[5,2]
			return scalar p_mtch_uci= r(table)[6,2]
			}
			
						if _rc!=0 {
				noi disp "No convergence achieved for Poisson conditional model in rep `i'"
				return scalar p_mtch_est= .
				return scalar p_mtch_se= .
				return scalar p_mtch_pval= .
				return scalar p_mtch_lci= .
				return scalar p_mtch_uci= .
			}	
			return list
			

end
