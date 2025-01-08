/*
* Define programs to run 1 repetitions
cap program drop matched_sim
cap program drop create_bc
program define create_bc, rclass
version 17
syntax [, nobs(int 100000) lambda_exp(real .00002909) gamma_out(real 1.19) lambda_out(real 0.012) loghr(real .0) gender_exp(real .0) age_exp(real .0) pracid_exp(real .0) gender_out(real .0) age_out(real .0) pracid_out(real .0) i(int 1) conf_type(string)] 

clear

return local bc_state = c(rngstate)
*/
cd "C:\Users\emsuewil\Documents\Work\PhD\Kirsty\mtchdcohortsim_paper\Do-files"
adopath + "C:\Users\emsuewil\Documents\Work\PhD\Kirsty\mtchdcohortsim_paper\Do-files"

set seed 2190890

clear
local nobs 1000000
local lambda_exp = 0.00003
local gender_exp = 0
local age_exp    = 0
local pracid_exp = 0
local lambda_out = 0.012
local gamma_out = 1.19
local loghr  = 1.098
local age_out = 0
local gender_out  = 0
local pracid_out = 0



********************************************
*  Obtain sample of nobs of exposed only   *
********************************************


forvalues i = 1 (1) 10 {
	clear
	set obs `nobs'

	//generate a patid
	 gen patid= (`i'-1)*`nobs' + _n

	//generate age (based on UK population)
	 gen age_start = rnormal(40, 12)
	 drop if age_start <18 

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

	* Keep exposed only
	drop if ever_exposed == 0
	
	gen indexdate = startdate + itime
	format indexdate %td

	replace indexdate=. if ever_exposed==0
	assert startdate<= indexdate & indexdate <= latest_date if indexdate!=.

	* Tidy dataset
	drop total_time total_time_days itime 
	order startdate, before(latest_date)

	* Save
	save exp_`i', replace
}

use exp_1, clear
forvalues i = 2 (1) 10 {
	append using exp_`i'
}
keep if _n<=`nobs'
assert _N>=`nobs'
forvalues i = 1 (1) 10 {
	erase  exp_`i'.dta
}


********************
*  Duplicate data  *
********************

drop age_start yob
gen age_index = (indexdate-dob)/365.25

rename ever_exposed exposed
rename patid setid
isid setid

gen f=2
expand f
drop f

bysort setid: replace exposed = 0 if _n==2
bysort setid (exposed): replace age_index = age_index+(uniform()-0.5) if _n==2
bysort setid (exposed): replace dob = indexdate - age_index*365.25  if _n==2
bysort setid (exposed): replace startdate = indexdate - uniform()*(indexdate-earliest_date) if _n==2

sort setid exposed
gen patid = _n
isid patid
isid setid exposed



**************************************
*  Generate outcomes from indexdate  *
**************************************

gen age_db_end = age_index + (latest_date - indexdate)/365.25

//generate outcome data. maximum time is now stipulated by t_years (total time per row)
survsim age_stime outcome, dist(weib) lambda(`lambda_out') gamma(`gamma_out') cov(exposed `loghr' age_index `age_out' gender `gender_out' pracid `pracid_out') ltruncated(age_index) maxtime(age_db_end)


//generate the date the outcome occurs
gen outcomedate = (dob + age_stime*365.25)
format outcomedate %td 
replace outcomedate=. if outcome==0

assert outcomedate<=latest_date | outcome==0


egen enddate = rowmin(latest_date outcomedate)
format enddate %td


******************
*  Analyse data  *
******************

	
stset enddate, origin(indexdate) enter(indexdate) fail(outcome) /// 
			id(patid) scale(365.25)

stcox i.exposed, nohr
stcox i.exposed age_index i.gender i.pracid, nohr 
stcox i.exposed, strata(setid) nohr 

stsplit survival, every(1)
gen y= _t - _t0
streg  i.exposed i.survival, dist(exp) nohr
streg  i.exposed age_index i.gender i.pracid i.survival, dist(exp) nohr  
xtpois _d i.exposed i.survival, fe i(setid) exp(y)


