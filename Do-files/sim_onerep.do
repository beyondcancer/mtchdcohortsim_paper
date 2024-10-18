///SETUP///

global resultsDir "H:\mtchdcohortsim"

cap log using "$resultsDir\simsum_onerep", replace

clear

cd $resultsDir

 //SET-UP the sim states	
set rng mt64s
set rngstream 1
set seed 88

 //SET-UP the simulation programs
qui do sim_setup.do


//SIM loop setup
global conf_test 1 2 3 4 5 6 //this is to set up the loop for different confounders
global ratio 1 3 5 10 //matching ratios 
global implementation //noreplacement // replacement //type of matching 


//number of simulations
global nsim_start = 1 // where the job will start from
global nsim = 1 //number of simulations in the loop for the number of base cohorts to be created



//SET-UP the DGM parameters
global n = 100000  //number of observations in each cohort
global loghr = 0 // true estimate
  


global lambda_exp .00003 .000015 .000021 .0000145 .000007 .0000097

global conf_expbi 0 0 0 0.39 0.39 0.39
global conf_outbi 0 0 0 0.1 0.1 0.1 
global conf_expcont 0 0 0 .0048 .0048 .0048
global conf_outcont 0 0 0 .00029 .00029 .00029

///Confounding
global conf_type 10_noconf 5_noconf 7_noconf 10_smallconf 5_smallconf 7_smallconf

/*
//global lambda_exp  .000029 .000014
//global lambda_exp .000015 .000008 .000004 .0000013
//global lambda_exp .000015 .000013 .000007 .0000025
global conf_expbi 0 0.4 0.7 1.60
global conf_outbi 0 0.28 0.4 1.09
global conf_expcont 0 0.0048 0.0095 0.018
global conf_outcont 0 0.0003 0.006 0.014
*/


///RUN SIM///
do sim_estimates

log close