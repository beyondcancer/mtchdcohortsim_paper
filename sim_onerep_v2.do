///SETUP///

global resultsDir "H:\mtchdcohortsim"

cap log using "$resultsDir\simsum_onerep.log", replace

clear

cd $resultsDir

 //SET-UP the sim states	
set rng mt64s
set rngstream 41
set seed 88

//number of simulations
global nsim_start = 1 // where the job will start from
global nsim = 1 //number of simulations in the loop for the number of base cohorts to be created

//SET-UP the DGM parameters
global n = 1000000  //number of observations in each cohort
global loghr = 0  // true estimate

//SIM loop setup
global ratio 1 3 5 10  //matching ratios 
global implementation replacement //noreplacement
 
// SET-UP exposure parameters 
global lambda_exp .0011 .0011 .0011 .0001 .0001 .0001
global gamma_exp 1.5 1.5 1.5 1.5 1.5 1.5

// SET-UP outcome parameters 
global lambda_out 0.0003 0.0003 0.0003 0.0003 0.0003 0.0003
global gamma_out 1.9 1.9 1.9 1.9 1.9 1.9

// SET-UP covariate  parameters 
global conf_test 1 2 3 4 5 6

global conf_expbi 0 0.93 1.609 0 0.93 1.609 
global conf_expmult 0 0.039 0.223 0 0.039 0.223 
global conf_outbi 0 0.672 1.386 0 0.672 1.386 
global conf_outmult	0 0.288 0.405 0 0.288 0.405 
global conf_type 10_nconf 10_sconf 1_bconf 1_nconf 1_sconf 1_bconf	

//SET-UP the simulation programs
do sim_setup_v2.do

///RUN SIM///

foreach l in $conf_test {
		
	forval i= $nsim_start/$nsim {	


		local exp_ithbi = word("$conf_expbi", `l')
		
		local out_ithbi = word("$conf_outbi", `l')
		local type_ith = word("$conf_type", `l')
		local lambda_ith = word("$lambda_exp", `l')
		local gamma_ith = word("$gamma_exp", `l')
		local lambda_o_ith = word("$lambda_out", `l')
		local gamma_o_ith = word("$gamma_out", `l')
		local exp_ithm = word("$conf_expmult", `l')
		local out_ithm = word("$conf_outmult", `l')

	create_bc, nobs($n) loghr($loghr) i(`i') lambda_exp(`lambda_ith') gamma_exp(`gamma_ith') lambda_out(`lambda_o_ith') gamma_out(`gamma_o_ith')gender_exp(`exp_ithbi') gender_out(`out_ithbi') pracid_exp(`exp_ithm') pracid_out(`out_ithm') conf_type("`type_ith'")

	//do sim_estimates_v2.do

		qui stset enddate, origin(dob) enter(date_entry) fail(outcome) /// 
				id(patid) scale(365.25)

	//KMPLOT

	sts graph, survival by(exposed) name(km_graph_`type_ith'_age, replace) title("origin(dob) enter(date_entry)")
	
	drop _*
	
			qui stset enddate, origin(startdate) enter(date_entry) fail(outcome) /// 
				id(patid) scale(365.25)
	
	sts graph, survival by(exposed) name(km_graph_`type_ith'_time, replace) title("origin(startdate) enter(date_entry)")
	
	graph combine km_graph_`type_ith'_age km_graph_`type_ith'_time, name(km_base_`type_ith', replace )
	
	
	}
}

log close