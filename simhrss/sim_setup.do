
* Define programs to run 1 repetitions
cap program drop matched_sim
cap program drop create_bc
program define create_bc, rclass
version 17
syntax [, nobs(int 100000) lambda_exp(real .00002909) gamma_out(real 1.19) lambda_out(real 0.012) loghr(real .0) gender_exp(real .0) age_exp(real .0) pracid_exp(real .0) gender_out(real .0) age_out(real .0) pracid_out(real .0) i(int 1) conf_type(string)] 

clear

return local bc_state = c(rngstate)

	
set obs `nobs'


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


//drop those with less than 1 yr of data in the dataset
//drop if total_time < 1

//generate exposure variable (as an exponential time to event as crc rates havent really changed since 1990s) and keep exposure date only for those with a cancer dx. 
* lamba = .00002909 gives ~20% exposed through whole study period, around 11% given staggered start times

survsim itime ever_exposed, dist(exp) lambda(`lambda_exp') cov(gender `gender_exp' age_start `age_exp' pracid `pracid_exp') maxtime(total_time_days)


* For descriptive purposes
gen itime_year = itime/365.25

gen indexdate = startdate + itime
format indexdate %td

replace indexdate=. if ever_exposed==0
assert startdate<= indexdate & indexdate <= latest_date if indexdate!=.

* Tidy dataset
drop total_time total_time_days itime itime_year
order startdate, before(latest_date)


****************
*  Split data  *
****************

* Aim: we want people who have a cancer diagnosis to have different hazards 
*      before and after cancer diagnosis


//generate dummy failure variable to split exposure data
gen fail=0

//Stset data to time to event data. 
stset latest_date, enter(startdate) origin(startdate) failure(fail) id(patid) 


* Keeping it in days: have removed option scale(365.25)
* Also, have changed name of x to pre_post_exposure to keep names explanatory
stsplit exposed, at(0) after(time=indexdate)
recode exposed 0=1 -1=0
label define pre_post 0 "Pre-exposure" 1 "Post-exposure"
label values exposed pre_post
replace exposed = 0 if ever_exposed==0


 * Remove vars just created to trick stata into doing what we want, to avoid confusion
drop fail
 
//keep _t0 and t from the stset as stored variables. Generate a variable which calculates the total_time within each row. 
gen t0=_t0
gen t=_t
gen ttime= t-t0

gen t0_years = t0/365.25
gen t_years = t/365.25
gen ttime_years= ttime/365.25


//Recategorise exposure as 0 in time period of startdate to indexdate & remove exposure date from the 1st row
//assert _t0==0 if exposed==0 // check pre-exposure time-period starts at t=0 


* Delete survival stuff to avoid confusion below
drop _*


//generate outcome data. maximum time is now stipulated by t_years (total time per row)
survsim stime outcome, dist(weib) lambda(`lambda_out') gamma(`gamma_out') cov(exposed `loghr' age_start `age_out' gender `gender_out' pracid `pracid_out') ltruncated(t0_years) maxtime(t_years)
* NB: when you add exposure effect, do it on time-varying exposure variable: exposed, not variable ever_exposed (which isn't time-varying) in this dataset

//generate the date the outcome occurs
gen outcomedate = (startdate + stime*365.25)
format outcomedate %td 
replace outcomedate=. if outcome==0


*** Retain one row per exposed person
bysort patid (t0): gen row = _n 
bysort patid: egen nout = sum(outcome)
bysort patid (t0): gen out_b4_dx = outcome[1]==1 

gen exclude = 0

* Necessary changes for the exposed pool:
replace indexdate=. if exposed==0
replace latest_date = latest_date - 1 if row==1 & ever_exposed==1 & exposed==0
drop if outcomedate>latest_date & outcome==1 //rare case where outcome is on indexdate
egen enddate = rowmin(latest_date outcomedate)
format enddate %td


// drop if outcome prior to indexdate
* Case (2)
drop if row==2 & ever_exposed==1 & out_b4_dx==1

// gen timein
gen timein = startdate 
* Case 1:
replace timein = startdate+t0 if row==2 & ever_exposed==1 & nout==0 
* Case 3:
replace timein = startdate+t0 if row==2 & ever_exposed==1 & out_b4_dx!=1 & nout==1

format timein %td

* Generate survival time (time to outcome or end study)
gen timeout = startdate + stime*365.25
format timeout %td


* Check dates are in correct order
assert startdate>=earliest_date
assert indexdate>=startdate
assert outcomedate>=startdate if exposed==0
assert outcomedate<=latest_date if outcome==1

*tidy dataset
drop t0 t ttime  t0_years t_years ttime_years stime row nout out_b4_dx earliest_date exclude timeout

order patid dob yob age_start gender pracid ///
	startdate indexdate timein 	///
	outcomedate  latest_date ///
	exposed outcome enddate 

	
	qui stset enddate, origin(startdate) enter(timein) fail(outcome) /// 
			id(patid) scale(365.25)
			
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
		

		
end 

				
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
	       qui stsplit survival, every(1)
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
