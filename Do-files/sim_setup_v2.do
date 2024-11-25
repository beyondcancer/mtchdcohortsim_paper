
////////////////////////////////////////////////////////
/////////////// 19NOV 2024/////////////////////////////
////////////////TWO TIMESCALES /////////////////////////
///////////////////////////////////////////////////////
quietly{

*Define programs to run 1 repetition
cap program drop matched_sim
cap program drop create_bc

*create_bc program
program define create_bc, rclass

//version of stata
version 17

//Set default values for simulation parameters

syntax [, nobs(int 100000) lambda_exp(real 0.0011) gamma_exp(real 1.5) gamma_out(real 1.9) lambda_out(real 0.0003) loghr(real .0) gender_exp(real .0) pracid_exp(real .0) gender_out(real .0) pracid_out(real .0) i(int 1) conf_type(string)] 

clear

//return random nuber state
return local bc_state = c(rngstate)

//set number of observations as per input to the program
set obs `nobs'

//generate unique patient id
gen patid= _n

//generate age based on a truncated normal distribution. Truncated at 18 but with a mean og 40.95 and a standard deviation 12.
gentrun double age_start, left(-1.8333333)
replace age_start= 12*age_start + 40

//Generate gender. Coded 1 & 2 to fit getmatched cohort program,
 gen gender = runiform()<=0.5
 replace gender = 2 if gender==0
label define sexlab 1 "male" 2 "female"
label values gender sexlab
recast byte gender 

//centre gender variable
gen centre_gender = (gender == 2) - 0.5

//generate gp variable & group from 0 - 13 practices (equivalent to CPRD if size were 100000)
gen gengp = runiform(0, 13)
gen pracid = round(gengp)
replace pracid=13 if pracid==0
drop gengp

// center the gp variable.
gen centre_pracid = (pracid - 7)

//gen min date and max-date of entering and leaving database so that the max time is 22 years
gen earliest_date = date("01jan2000", "DMY")
gen latest_date = date("31dec2021", "DMY")
format earliest_date %td
format latest_date %td

//generate random startdates within this period
gen startdate= earliest_date + (latest_date - earliest_date) * runiform()
format startdate %td

//generate year of birth - for use in the getmatchedcohort program
gen dob = startdate - (age_start * 365.25)
format dob %td
gen yob= year(dob)

//calculations for age as timescale
gen age_exit = (latest_date - dob)/365.25

//generate exposure variable changing exposure hazard for attained age by using age as timescale 
noi di("Exposure generation: lambda is: `lambda_exp', gamma is `gamma_exp'")

return local pre_state = c(rngstate)

//Need to use weibull to account for increasing hazards with age. age_itime (is age in years at exposure)
survsim age_itime ever_exposed, dist(weibull) 					///
	lambda(`lambda_exp') gamma(`gamma_exp') 					///
	cov(centre_gender `gender_exp' centre_pracid `pracid_exp') ltruncated(age_start) maxtime(age_exit)

// generate indexdate, and remove for those whose exposure dates are not within the study period
gen indexdate = dob + (age_itime * 365.25)
format indexdate %td "DDMMMYYYY"
replace indexdate=. if ever_exposed==0

* Tidy dataset
drop  age_itime
order startdate, before(latest_date)

//avoid issues with those that have a startdate on the same day as their enddate
//replace startdate = startdate - 1 if startdate==latest_date


****************
*  Split data  *
****************

// AIM: people who have a cancer diagnosis to have different hazards before and after cancer diagnosis

//generate dummy failure variable to split exposure data
gen fail=0

//Stset data with age as underlying timescale
gen late_date = latest_date
qui stset late_date, enter(startdate) origin(startdate) failure(fail) id(patid) 

// Split data at date of exposure (for those who become exposed)
stsplit exposed, at(0) after(time=indexdate)

//Recategorise exposure as 0 in time period of startdate to indexdate & mark as unexposed in the 1st row
recode exposed 0=1 -1=0
replace exposed = 0 if ever_exposed==0
label define pre_post 0 "Pre-exposure" 1 "Post-exposure"
label values exposed pre_post

*Generate future time since study start time scale

//recalculate age at entry and exit for each row
gen date_entry = (_origin + _t0)
gen date_exit  = (_origin + _t)
format date_entry date_exit %td
gen age_entry     = (date_entry - dob)/365.25
replace age_exit  = (date_exit - dob)/365.25

 * tidy unawanted variables
drop fail _* late_date
order age_start, before(age_entry)
order age_exit, after(age_entry) 

//generate outcome variable with age as underlying timescale

noi di("Outcome generation lambda is: `lambda_out', gamma is `gamma_out'")

survsim age_stime outcome, dist(weib) 				///
	lambda(`lambda_out') gamma(`gamma_out') 		///
	cov(exposed `loghr' centre_gender `gender_out'  ///
	centre_pracid `pracid_out') 					///
	ltruncated(age_entry) maxtime(age_exit)
	 

 // No longer need centred covariates
drop centre*

//generate the date the outcome occurs
gen outcomedate = dob + (age_stime * 365.25)
format outcomedate %td 
replace outcomedate=. if outcome==0

//Remove indexdate for the rows of unexposed time
replace indexdate=. if exposed==0

//Count the number of records and outcomes for each patient
bysort patid (age_entry): gen row = _n 
bysort patid: egen nout = sum(outcome)

//Identify patients with an outcome before indexdate
bysort patid (age_entry): gen out_b4_dx = outcome[1]==1 

//Necessary changes for the exposed pool:

//Make sure that end date is 1 day before indexdate for row with unexposed time-period
replace date_exit = date_exit - 1 if row==1 & ever_exposed==1 & exposed==0

// to avoid errors with in/out dates being equal for survival
replace date_exit = date_exit + 0.5 if row==2 & ever_exposed==1 & exposed==1

// drop time after exposure if outcome is prior to exposure (data is censored at first outcome so this row is meaningless)
noi disp("Dropping time after outcome")
drop if row==2 & ever_exposed==1 & out_b4_dx==1

//Generate censoring date: first of end of data or outcome date
egen enddate = rowmin(date_exit outcomedate)
format enddate %td

// Tidy data
drop row nout out_b4_dx
drop age_entry age_exit age_stime

order patid dob yob age_start gender pracid ///
	startdate indexdate ///
	outcomedate latest_date ///
	exposed outcome enddate 

	// ANALYTICAL MODELS FOR BASE COHORT
	
	// *AGE AS TIMESCALE	
	qui stset enddate, origin(dob) enter(date_entry) fail(outcome) /// 
			id(patid) scale(365.25)
				
						capture { 
			qui  stcox i.exposed, nohr
			return scalar abase_uadj_est= r(table)[1,2] 
			return scalar abase_uadj_se= r(table)[2,2]
			return scalar abase_uadj_pval= r(table)[4,2]
			return scalar abase_uadj_lci= r(table)[5,2]
			return scalar abase_uadj_uci= r(table)[6,2]
			}
			
			if _rc!=0 {
				noi disp "No convergence achieved for abase Cox unadjusted model in rep `i'"
		    return scalar abase_uadj_est= .
			return scalar abase_uadj_se= .
			return scalar abase_uadj_pval= .
			return scalar abase_uadj_lci= .
			return scalar abase_uadj_uci= .
			}
			
			return list 
			capture{
			qui stcox i.exposed i.gender i.pracid, nohr
		
			return scalar abase_adj_est= r(table)[1,2] 
			return scalar abase_adj_se= r(table)[2,2] 
			return scalar abase_adj_pval= r(table)[4,2]
			return scalar abase_adj_lci= r(table)[5,2]
			return scalar abase_adj_uci= r(table)[6,2]
			}
				
			if _rc!=0 {
				noi disp "No convergence achieved for abase Cox adjusted model in rep `i'"
		    return scalar abase_adj_est= .
			return scalar abase_adj_se= .
			return scalar abase_adj_pval= .
			return scalar abase_adj_lci= .
			return scalar abase_adj_uci= .
			}
			
			return list 
			
		// Create stsplit by age
		stsplit age_group_base, every(5)
			
		drop _*	
	
	
	// *TIME IN STUDY AS TIMESCALE

	qui stset enddate, origin(startdate) enter(date_entry) fail(outcome) /// 
			id(patid) scale(365.25)
			
				
						capture { 
			qui  stcox i.exposed, nohr
			return scalar tbase_uadj_est= r(table)[1,2] 
			return scalar tbase_uadj_se= r(table)[2,2]
			return scalar tbase_uadj_pval= r(table)[4,2]
			return scalar tbase_uadj_lci= r(table)[5,2]
			return scalar tbase_uadj_uci= r(table)[6,2]
			}
			
			if _rc!=0 {
				noi disp "No convergence achieved for tbase Cox unadjusted model in rep `i'"
		    return scalar tbase_uadj_est= .
			return scalar tbase_uadj_se= .
			return scalar tbase_uadj_pval= .
			return scalar tbase_uadj_lci= .
			return scalar tbase_uadj_uci= .
			}
			
			capture{
			qui stcox i.exposed i.age_group_base i.gender i.pracid, nohr
		
			return scalar tbase_adj_est= r(table)[1,2] 
			return scalar tbase_adj_se= r(table)[2,2] 
			return scalar tbase_adj_pval= r(table)[4,2]
			return scalar tbase_adj_lci= r(table)[5,2]
			return scalar tbase_adj_uci= r(table)[6,2]
			}
				
			if _rc!=0 {
				noi disp "No convergence achieved for tbase Cox adjusted model in rep `i'"
		    return scalar tbase_adj_est= .
			return scalar tbase_adj_se= .
			return scalar tbase_adj_pval= .
			return scalar tbase_adj_lci= .
			return scalar tbase_adj_uci= .
			}
			
		// Get rid of base cohort stsplit and stset variabes
		drop age_group_base
		stjoin
		drop _*
		
		//Count total number of patients in data
	    bysort patid: gen tag = _n == 1
		count if tag 
		return scalar total_n = r(N)
		
		save "bc_`i'_`conf_type'", replace
		
end 
				
///////

// Matching program
program define matched_sim, rclass
version 17

syntax [, getallpossible(string) ratio(int 5) i(int 1) conf_type(string)]
             
			qui use "bc_`i'_`conf_type'", clear

	//No replacement
	if "`getallpossible'" == "noreplacement"{	
				
			qui return local state_`getallpossible'_`ratio' = c(rngstate)

	// Count time it takes to run getmatchedcohort program
timer on 1
		
			 getmatchedcohort, cprddb(gold) practice gender yob yobwindow(1) ctrlsperexp(`ratio') savedir($workdir) filesuffix(`i'`conf_type'`getallpossible') followup dontcheck
			
			qui use "getmatchedcohort`i'`conf_type'`getallpossible'"
			qui return scalar nomatch= r(nomtch)	
			qui merge 1:1 patid exposed using "bc_`i'_`conf_type'" 
			
			qui keep if _merge==3 
			
			qui gen age_index= (indexdate-dob)/365.25
			
			save "temp_`i'`conf_type'`getallpossible'.dta", replace
			
timer off 1
noi disp "Elapsed time for getmatchedcohort in rep `i': "
timer list  
return scalar mtch_time = r(t1)
timer clear 1

			}
				
	//Replacement

				if "`getallpossible'" == "replacement"{
				
			qui return local state_`getallpossible'_`ratio' = c(rngstate) 
			
timer on 1			
			 getmatchedcohort, cprddb(gold) practice gender yob yobwindow(1) savedir($workdir) getallpossible filesuffix(`i'`conf_type'`getallpossible') followup dontcheck
				
			qui use "getmatchedcohort`i'`conf_type'`getallpossible'"
			qui return scalar nomatch= r(nomtch)
			qui merge m:1 patid exposed using "bc_`i'_`conf_type'" 
			
			qui keep if _merge==3 
			
			// only keep the amount of matches we need
				gen count=`ratio'+1
				gen rand = runiform()
				gsort setid -exposed  rand
				by setid, sort: keep if exposed==1| _n<=count
			    gen og_patid= patid
				replace patid = _n
				
			save "temp_`i'`conf_type'`getallpossible'.dta", replace
			
timer off 1
noi disp "Elapsed time for getmatchedcohort in rep `i': "
timer list  
return scalar mtch_time = r(t1)

timer clear 1

			gen age_index= (indexdate-dob)/365.25
			 

			}
				
			*ANALYTICAL MODELS FOR MATCHED COHORT ANALYS@IS

			*AGE AS TIMESCALE*
			qui stset enddate, origin(dob) enter(indexdate) fail(outcome) id(patid) scale(365.25)
			
			//split data by age 
			* Use stsplit to split the data (every 5 years of age)
			stsplit age_group, every(5)
			gen y_age = _t - _t0
			
			*unadjusted
			capture { 
			qui stcox i.exposed, nohr
			return scalar cox_uadj_a_est= r(table)[1,2] 
			return scalar cox_uadj_a_se= r(table)[2,2]
			return scalar cox_uadj_a_pval= r(table)[4,2]
			return scalar cox_uadj_a_lci= r(table)[5,2]
			return scalar cox_uadj_a_uci= r(table)[6,2]
			}
			
			if _rc!=0 {
				noi disp "No convergence achieved for Cox unadjusted model in rep `i' (age as timescale)"
		    return scalar cox_uadj_a_est= .
			return scalar cox_uadj_a_se= .
			return scalar cox_uadj_a_pval= .
			return scalar cox_uadj_a_lci= .
			return scalar cox_uadj_a_uci= .
			}
			
			*adjusted for matching variables
			capture { 
			qui stcox i.exposed i.gender i.pracid, nohr 
			return scalar cox_adj_a_est= r(table)[1,2] 
			return scalar cox_adj_a_se= r(table)[2,2]
			return scalar cox_adj_a_pval= r(table)[4,2]
			return scalar cox_adj_a_lci= r(table)[5,2]
			return scalar cox_adj_a_uci= r(table)[6,2]
			
			}
						if _rc!=0 {
				noi disp "No convergence achieved for Cox adjusted model in rep `i' (age as timescale)"
			return scalar cox_adj_a_est= .
			return scalar cox_adj_a_se= .
			return scalar cox_adj_a_pval= .
			return scalar cox_adj_a_lci= .
			return scalar cox_adj_a_uci= .
			}
			
			*adjusted by matched set
			capture { 
			qui stcox i.exposed, strata(setid) nohr
			return scalar cox_mtch_a_est= r(table)[1,2]
			return scalar cox_mtch_a_se= r(table)[2,2]
			return scalar cox_mtch_a_pval= r(table)[4,2]
			return scalar cox_mtch_a_lci= r(table)[5,2]
			return scalar cox_mtch_a_uci= r(table)[6,2]
			}
			
						if _rc!=0 {
				noi disp "No convergence achieved for Cox adjusted model in rep `i' (age as timescale)"
				return scalar cox_mtch_a_est= .
				return scalar cox_mtch_a_se= .
				return scalar cox_mtch_a_pval= .
				return scalar cox_mtch_a_lci= .
				return scalar cox_mtch_a_uci= .
				
			}
				

			*Poisson regressions
				*unadjusted
			capture {
			qui streg  i.exposed i.age_group, dist(exp) nohr
			return scalar p_uadj_a_est= r(table)[1,2]
			return scalar p_uadj_a_se= r(table)[2,2]
			return scalar p_uadj_a_pval= r(table)[4,2]
			return scalar p_uadj_a_lci= r(table)[5,2]
			return scalar p_uadj_a_uci= r(table)[6,2]
			}
			
						if _rc!=0 {
				noi disp "No convergence achieved for Poisson unadjusted model in rep `i' (age as timescale)"
			return scalar p_uadj_a_est= .
			return scalar p_uadj_a_se= .
			return scalar p_uadj_a_pval= .
			return scalar p_uadj_a_lci= .
			return scalar p_uadj_a_uci= .
			}
			
			*adjusted for matching variables 
			capture {
			 qui streg  i.exposed i.gender i.pracid i.age_group, dist(exp) nohr
			return scalar p_adj_a_est= r(table)[1,2]
			return scalar p_adj_a_se= r(table)[2,2]
			return scalar p_adj_a_pval= r(table)[4,2]
			return scalar p_adj_a_lci= r(table)[5,2]
			return scalar p_adj_a_uci= r(table)[6,2]
			}
									if _rc!=0 {
				noi disp "No convergence achieved for Poisson adjusted model in rep `i' (age as timescale)"
				return scalar p_adj_a_est= .
				return scalar p_adj_a_se= .
				return scalar p_adj_a_pval= .
				return scalar p_adj_a_lci= .
				return scalar p_adj_a_uci= .
			}
			
			*conditional poisson*
			
			capture { 
	
			qui xtpois _d i.exposed i.age_group, fe i(setid) exp(y_age)
			return scalar p_mtch_a_est= r(table)[1,2]
			return scalar p_mtch_a_se= r(table)[2,2]
			return scalar p_mtch_a_pval= r(table)[4,2]
			return scalar p_mtch_a_lci= r(table)[5,2]
			return scalar p_mtch_a_uci= r(table)[6,2]
			}
			
						if _rc!=0 {
				noi disp "No convergence achieved for Poisson conditional model in rep `i' (age as timescale)"
				return scalar p_mtch_a_est= .
				return scalar p_mtch_a_se= .
				return scalar p_mtch_a_pval= .
				return scalar p_mtch_a_lci= .
				return scalar p_mtch_a_uci= .
			}	
			
			//drop age as timescale stset vars
			drop _*
			
			*TIME IN STUDY AS TIMESCALE*
			qui stset enddate, origin(indexdate) enter(indexdate) fail(outcome) id(patid) scale(365.25)
			
			// create stsplit for time in study (every year of follow up)
			stsplit survival, every(1)
			gen y_surv = _t - _t0
			

			*unadjusted
			capture { 
			qui stcox i.exposed, nohr
			return scalar cox_uadj_t_est= r(table)[1,2] 
			return scalar cox_uadj_t_se= r(table)[2,2]
			return scalar cox_uadj_t_pval= r(table)[4,2]
			return scalar cox_uadj_t_lci= r(table)[5,2]
			return scalar cox_uadj_t_uci= r(table)[6,2]
			}
			
			if _rc!=0 {
				noi disp "No convergence achieved for Cox unadjusted model in rep `i' (survival as timescale)"
		    return scalar cox_uadj_t_est= .
			return scalar cox_uadj_t_se= .
			return scalar cox_uadj_t_pval= .
			return scalar cox_uadj_t_lci= .
			return scalar cox_uadj_t_uci= .
			}
			
			*adjusted for matching variables
			capture { 
			qui stcox i.exposed i.gender i.pracid i.age_group, nohr 
			return scalar cox_adj_t_est= r(table)[1,2] 
			return scalar cox_adj_t_se= r(table)[2,2]
			return scalar cox_adj_t_pval= r(table)[4,2]
			return scalar cox_adj_t_lci= r(table)[5,2]
			return scalar cox_adj_t_uci= r(table)[6,2]
			
			}
						if _rc!=0 {
				noi disp "No convergence achieved for Cox adjusted model in rep `i' (survival as timescale)"
			return scalar cox_adj_t_est= .
			return scalar cox_adj_t_se= .
			return scalar cox_adj_t_pval= .
			return scalar cox_adj_t_lci= .
			return scalar cox_adj_t_uci= .
			}
			
			*adjusted by matched set
			capture { 
			qui stcox i.exposed i.age_group, strata(setid) nohr
			return scalar cox_mtch_t_est= r(table)[1,2]
			return scalar cox_mtch_t_se= r(table)[2,2]
			return scalar cox_mtch_t_pval= r(table)[4,2]
			return scalar cox_mtch_t_lci= r(table)[5,2]
			return scalar cox_mtch_t_uci= r(table)[6,2]
			}
			
						if _rc!=0 {
				noi disp "No convergence achieved for Cox adjusted model in rep `i' (survival as timescale)"
				return scalar cox_mtch_t_est= .
				return scalar cox_mtch_t_se= .
				return scalar cox_mtch_t_pval= .
				return scalar cox_mtch_t_lci= .
				return scalar cox_mtch_t_uci= .
				
			}
				

			*Poisson regressions
				*unadjusted
			capture {
			qui streg  i.exposed i.survival i.age_group, dist(exp) nohr
			return scalar p_uadj_t_est= r(table)[1,2]
			return scalar p_uadj_t_se= r(table)[2,2]
			return scalar p_uadj_t_pval= r(table)[4,2]
			return scalar p_uadj_t_lci= r(table)[5,2]
			return scalar p_uadj_t_uci= r(table)[6,2]
			}
			
						if _rc!=0 {
				noi disp "No convergence achieved for Poisson unadjusted model in rep `i' (survival as timescale)"
			return scalar p_uadj_t_est= .
			return scalar p_uadj_t_se= .
			return scalar p_uadj_t_pval= .
			return scalar p_uadj_t_lci= .
			return scalar p_uadj_t_uci= .
			}
			
			*adjusted for matching variables 
			capture {
			 qui streg  i.exposed gender pracid i.survival i.age_group, dist(exp) nohr
			return scalar p_adj_t_est= r(table)[1,2]
			return scalar p_adj_t_se= r(table)[2,2]
			return scalar p_adj_t_pval= r(table)[4,2]
			return scalar p_adj_t_lci= r(table)[5,2]
			return scalar p_adj_t_uci= r(table)[6,2]
			}
									if _rc!=0 {
				noi disp "No convergence achieved for Poisson adjusted model in rep `i' (survival as timescale)"
				return scalar p_adj_t_est= .
				return scalar p_adj_t_se= .
				return scalar p_adj_t_pval= .
				return scalar p_adj_t_lci= .
				return scalar p_adj_t_uci= .
			}
			
			*conditional poisson*
			
			capture { 
	
			qui xtpois _d i.exposed i.age_group i.survival, fe i(setid) exp(y_surv)
			return scalar p_mtch_t_est= r(table)[1,2]
			return scalar p_mtch_t_se= r(table)[2,2]
			return scalar p_mtch_t_pval= r(table)[4,2]
			return scalar p_mtch_t_lci= r(table)[5,2]
			return scalar p_mtch_t_uci= r(table)[6,2]
			}
			
						if _rc!=0 {
				noi disp "No convergence achieved for Poisson conditional model in rep `i' (survival as timescale)"
				return scalar p_mtch_t_est= .
				return scalar p_mtch_t_se= .
				return scalar p_mtch_t_pval= .
				return scalar p_mtch_t_lci= .
				return scalar p_mtch_t_uci= .
			}	
			//drop stset vars 
			drop _*
			

end

}
