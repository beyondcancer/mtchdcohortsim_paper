
////////////////////////////////////////////////////////
/////////////// 22 OCT 2024/////////////////////////////
///////CHANGES AS PER MEETING WITH FIZZ/////////////////
///////////////////////////////////////////////////////

*Define programs to run 1 repetition

cap program drop matched_sim
cap program drop create_bc

program define create_bc, rclass

//version of stata
version 17

//Set default values for simulation parameters //LOOK AT STATS WITH JIM

syntax [, nobs(int 100000) lambda_exp(real 1.0) gamma_exp(real 1.0) gamma_out(real 1.19) lambda_out(real 0.012) loghr(real .0) gender_exp(real .0) age_exp(real .0) pracid_exp(real .0) gender_out(real .0) age_out(real .0) pracid_out(real .0) i(int 1) conf_type(string)] 

//lambda_exp(real 100) gamma_exp(real 0.5)

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
sum age_start

//center age variable
di 40 + normalden((18-40)/12)*12/(1-normal((18-40)/12)) // calculates the mean of the truncated variable. RESULT: 40.922525
gen age_start_centered = age_start - 40.922525

//Generate gender. Coded 1 & 2 to fit getmatched cohort program,
 gen gender = runiform()<=0.5
 replace gender = 2 if gender==0
label define sexlab 1 "male" 2 "female"
label values gender sexlab
recast byte gender 

//centre gender variable
gen gender_centered = (gender == 2) - 0.5

//generate gp variable & group from 0 - 13 practices (equivalent to CPRD if size were 100000)
gen gengp = runiform(0, 13)
gen pracid = round(gengp)
replace pracid=13 if pracid==0
drop gengp
tab pracid 

// center the gp variable. Center is 8, and the maximum amounf of variation is half the level of confounding in the binary variable
gen pracid_centered = (pracid - 8) * ((0.39/2)/6)

//gen min date and max-date of entering and leaving database so that the max time is 22 years
gen earliest_date = date("01jan2000", "DMY")
gen latest_date = date("31dec2021", "DMY")
format earliest_date %td
format latest_date %td

//generate random startdates within this period
gen startdate= earliest_date + (latest_date - earliest_date) * runiform()
format startdate %td

//generate year of birth - for use in the getmatchedcohort program
cap drop dob
gen dob = startdate - (age_start * 365.25)
format dob %td
gen yob= year(dob)

//Calculate total time between startdate and latest_date (in years)
gen total_time = (latest_date-startdate)/365.25
gen total_time_days = latest_date-startdate

//calculations for age as timescale
gen age_exit = (latest_date - dob)/365.25

//generate exposure variable canging exposure hazard for attained age by using age as timescale // QUESTION: WOULD AGE BE REMOVED HERE TOO IF AGE IS TIMESCALE
// Using age at timescale allows for increasing hazards with attainded age

noi di("lambda is: `lambda_exp', gamma is `gamma_exp'")

//Need to use weibull to account for increasing hazards with age. age_itime (is age in years at exposure)
survsim age_itime ever_exposed, dist(weibull) lambda(`lambda_exp') gamma(`gamma_exp') cov(gender_centered `gender_exp' age_start_centered `age_exp' pracid_centered `pracid_exp') ltruncated(age_start) maxtime(age_exit)


// generate indexdate, and remove for those whose exposure dates are not within the study period
gen indexdate = dob + (age_itime * 365.25)
format indexdate %td

replace indexdate=. if ever_exposed==0

// descriptive statistics for exposed
tab ever_exposed
sum age_itime

local mean = r(mean)

hist age_itime, name("exposure", replace) 
hist age_itime, by(ever_exposed) name("exposure_exp", replace)

br
hist total_time, name("total_time", replace)
assert startdate<= indexdate & indexdate <= latest_date if indexdate!=.

* Tidy dataset
drop total_time total_time_days age_itime
order startdate, before(latest_date)


****************
*  Split data  *
****************

// AIM: people who have a cancer diagnosis to have different hazards before and after cancer diagnosis


//generate dummy failure variable to split exposure data
gen fail=0

//Stset data with age as underlying timescale
stset latest_date, enter(startdate) origin(startdate) failure(fail) id(patid) 


stsplit exposed, at(0) after(time=indexdate)

//Recategorise exposure as 0 in time period of startdate to indexdate &  mark as unexposed in the 1st row
recode exposed 0=1 -1=0
//Also, have changed name of x to pre_post_exposure to keep names explanatory
label define pre_post 0 "Pre-exposure" 1 "Post-exposure"
label values exposed pre_post
replace exposed = 0 if ever_exposed==0
assert _t0==0 if exposed==0 //check pre-exposure time-period starts at t=0 

*Generate future time since study start time scale

//keep _t0 and t from the stset as stored variables. Generate a variable which calculates the total_time within each row. 

gen t0=_t0
gen t=_t
gen ttime= t-t0

gen t0_years = t0/365.25 
gen t_years = t/365.25
gen ttime_years= ttime/365.25

 * tidy unawanted stset generated variables
drop fail
drop _*

//recalculate age_exit for the each row
replace age_exit = (latest_date - dob)/365.25

noi di("lambda is: `lambda_out', gamma is `gamma_out'")

//generate outcome variable with age as underlying timescale
 survsim age_stime outcome, dist(weib) lambda(`lambda_out') gamma(`gamma_out') cov(exposed `loghr' gender_centered `gender_out' pracid_centered `pracid_out') ltruncated(age_start) maxtime(age_exit)

 
 hist age_stime, name("outcome", replace)
 hist age_stime, by(ever_exposed) name("outcome_evexp", replace)
 hist age_stime, by(exposed) name("outcome_exp", replace)

//generate the date the outcome occurs
gen outcomedate = dob + (age_stime * 365.25)
format outcomedate %td 
replace outcomedate=. if outcome==0

noi disp("% of sample with and without outcome (before and after exposure)")
tab outcome 


//Order and track the number of records for each patient
bysort patid (t0): gen row = _n 

//Generate the number of outcomes per patient
bysort patid: egen nout = sum(outcome)

//Identify patients with an outcome before indexdate
bysort patid (t0): gen out_b4_dx = outcome[1]==1 

noi disp("% of sample with outcome and exposure ")
tab out_b4_dx if ever_exposed ==1


//Necessary changes for the exposed pool:

//Remove indexdate for the rows of unexposed time
replace indexdate=. if exposed==0

//Make sure that end date is 1 day before indexdate for row with unexposed time-period
replace latest_date = latest_date - 1 if row==1 & ever_exposed==1 & exposed==0

//rare case where outcome is on indexdate:
noi di ("These people have outcome on index date: ")
list if outcomedate>=latest_date & outcome==1 
noi di ("Dropping them from analysis")
drop if outcomedate>=latest_date & outcome==1 

//Generate censoring date: first of end of data or outcome date
egen enddate = rowmin(latest_date outcomedate)
format enddate %td

//replcae age_exit using enddate
replace age_exit = (enddate - dob)/365.25

noi disp ("Tabulation of number of exposed with outcome")
tab exposed outcome

//tab how many have an outcome = 1 in row 1 and row 2 by exposed
noi disp ("Data after outcome has happened - 535 people have an outcome before and after exposure ")
 tab exposed outcome if row==2 & ever_exposed==1 & out_b4_dx==1 

// drop time after exposure if outcome is prior to exposure (data is censored at first outcome so this row is meaningless)
noi disp("Dropping time after outcome")
drop if row==2 & ever_exposed==1 & out_b4_dx==1

// Ascribe study startdate: 
//Same as startdate for those that are unexposed for the whole timeperiod.
gen timein = startdate 

//Enters study at date of exposure if:
//individual is exposed but no outcome
replace timein = startdate+t0 if row==2 & ever_exposed==1 & nout==0 
//or
//individual is exposed and outcome happens after exposure
replace timein = startdate+t0 if row==2 & ever_exposed==1 & out_b4_dx!=1 & nout==1

format timein %td

*tidy dataset
drop t0 t ttime  t0_years t_years ttime_years age_stime row nout out_b4_dx 

* Check dates are in correct order
assert startdate>=earliest_date
assert indexdate>=startdate
assert outcomedate>=startdate if exposed==0
assert outcomedate<=latest_date if outcome==1


order patid dob yob age_start gender pracid ///
	startdate indexdate ///
	outcomedate latest_date ///
	exposed outcome enddate 

	
	

	 stset enddate, origin(startdate) enter(timein) fail(outcome) /// 
			id(patid) scale(365.25)
			
						capture { 
			noi  stcox i.exposed 
			noi stcox i.ever_exposed
			 //, nohr
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
			noi stcox i.exposed age_start gender pracid

			 //, nohr
		
			return scalar base_adj_est= r(table)[1,2] 
			return scalar base_adj_se= r(table)[2,2] 
			return scalar base_adj_pval= r(table)[4,2]
			return scalar base_adj_lci= r(table)[5,2]
			return scalar base_adj_uci= r(table)[6,2]
			}
			
						capture{
			 noi stcox i.ever_exposed age_start gender pracid

			 //, nohr
		
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
