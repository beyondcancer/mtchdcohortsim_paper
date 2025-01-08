
clear


cd $workdir

//SET-UP the sim states
set rng mt64s
set rngstream 401
set seed 88

//number of simulations
global nsim_start = 1 // where the job will start from
global nsim = 2000 //number of simulations in the loop for the number of base cohorts to be created

//SET-UP the DGM parameters
global n = 100000  //number of observations in each cohort
global loghr = 0  // true estimate

//SIM loop setup
global ratio 1 3 5 10  //matching ratios
global implementation noreplacement replacement
// SET-UP exposure parameters
global lambda_exp .0001 
global gamma_exp 1.5

// SET-UP outcome parameters
global lambda_out 0.0003 
global gamma_out 1.9 

// SET-UP covariate parameters
global conf_test 1
global conf_expbi 0.93 
global conf_expmult 0.039  
global conf_outbi 0.672
global conf_outmult	0.288 
global conf_type 1_sconf

//SET-UP the simulation programs
do sim_setup_v2.do

///RUN SIM///
do sim_estimates_v2.do

