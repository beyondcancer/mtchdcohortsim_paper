//version of stata
version 17

* Parameters
local nobs 			= 100000
local lambda_exp 	= 0.001 
local gamma_exp 	= 1.5
local gamma_out 	= 1.9 
local lambda_out 	= 0.0003 
local loghr 		= 0.7 
local gender_exp 	= 2*0.39 
local pracid_exp 	= 2*0.01625
local gender_out 	= 2*0.28 
local pracid_out 	= 2*0.012 


clear

//set number of observations as per input to the program
set obs `nobs'

//generate unique patient id
gen patid = _n

//generate age based on a truncated normal distribution. 
// Truncated at 18
// Underlying normal mean of 40 and a standard deviation 12.
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
//Need to use weibull to account for increasing hazards with age. age_itime (is age in years at exposure)
survsim age_itime ever_exposed, dist(weibull) 					///
	lambda(`lambda_exp') gamma(`gamma_exp') 					///
	cov(centre_gender `gender_exp' centre_pracid `pracid_exp') 	///
	ltruncated(age_start) maxtime(age_exit)

// generate indexdate, and remove for those whose exposure dates are not within the study period
gen indexdate = dob + (age_itime * 365.25)
format indexdate %td
replace indexdate=. if ever_exposed==0

* Check order of individual study entry, exposure, and study end date
assert startdate<= indexdate & indexdate <= latest_date if indexdate!=.

* Tidy dataset
drop age_itime
order startdate, before(latest_date)




****************
*  Split data  *
****************

// AIM: people who have a cancer diagnosis to have different hazards before and after cancer diagnosis

//generate dummy failure variable to split exposure data
gen fail=0

//Stset data with days since start date as underlying timescale
gen late_date = latest_date
stset late_date, enter(startdate) origin(startdate) failure(fail) id(patid) 
// Split data at date of exposure (for those who become exposed)
stsplit exposed, at(0) after(time=indexdate)

//Recategorise exposure as 0 in time period of startdate to indexdate &  mark as unexposed in the 1st row
recode exposed 0=1 -1=0
replace exposed = 0 if ever_exposed==0
label define pre_post 0 "Pre-exposure" 1 "Post-exposure"
label values exposed pre_post
assert _t0==0 if exposed==0 //check pre-exposure time-period starts at t=0 

//recalculate age at entry and exit for each row
gen date_entry = (_origin + _t0)
gen date_exit  = (_origin + _t)
format date_entry date_exit %td
gen age_entry     = (date_entry - dob)/365.25
replace age_exit  = (date_exit - dob)/365.25

 * Tidy unawanted stset generated variables
assert date_exit==late_date
drop fail _* late_date
order age_start, before(age_entry)
order age_exit, after(age_entry) 

//generate outcome variable with age as underlying timescale
survsim age_stime outcome, dist(weib) 				///
	lambda(`lambda_out') gamma(`gamma_out') 		///
	cov(exposed `loghr' centre_gender `gender_out'  ///
	centre_pracid `pracid_out') 					///
	ltruncated(age_entry) maxtime(age_exit)

// No longer need centred covariates
drop centre_*

//generate the date the outcome occurs
gen outcomedate = dob + (age_stime * 365.25)
format outcomedate %td 
replace outcomedate=. if outcome==0

//Remove indexdate for the rows of unexposed time
replace indexdate=. if exposed==0


//Count the number of records and outcomes for each patient
bysort patid (age_entry): gen row = _n 
bysort patid: egen nout = sum(outcome)

//Identify patients with an outcome before indexdate (exposure)
bysort patid (age_entry): gen out_b4_dx = outcome[1]==1 

//Reset the end date to 1/2 day before/after indexdate for the two rows per person 
// to avoid errors with in/out dates being equal for survival   
replace date_exit = date_exit - 0.5 if row==1 & ever_exposed==1 & exposed==0
replace date_exit = date_exit + 0.5 if row==2 & ever_exposed==1 & exposed==1

// drop time after exposure if outcome is prior to exposure (data is censored at first outcome so this row is meaningless)
drop if row==2 & ever_exposed==1 & out_b4_dx==1
// rare case where outcome is on indexdate:
drop if outcomedate>=latest_date & outcome==1 



//Generate censoring date: first of end of data or outcome date
egen enddate = rowmin(date_exit outcomedate)
format enddate %td

// Tidy data
drop row nout out_b4_dx
drop age_entry age_exit age_stime


* Check dates are in correct order
assert startdate>=earliest_date
assert indexdate>=startdate
assert outcomedate>=startdate if exposed==0
assert outcomedate<=latest_date if outcome==1


order patid dob yob age_start gender pracid ///
	startdate indexdate ///
	outcomedate latest_date ///
	exposed outcome enddate 
	

stset enddate, origin(startdate) enter(date_entry) fail(outcome) /// 
			id(patid) scale(365.25)
			
stcox i.exposed, nohr 
stcox i.exposed age_start i.gender i.pracid, nohr

