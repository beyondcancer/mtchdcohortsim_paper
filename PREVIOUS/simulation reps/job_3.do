
clear





cd $workdir

//SET-UP the sim states
set rng mt64s
set rngstream 201
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
global lambda_exp .0011 
global gamma_exp 1.5

// SET-UP outcome parameters
global lambda_out 0.0003 
global gamma_out 1.9 

// SET-UP covariate parameters
global conf_test 1

global conf_expbi 1.609 
global conf_expmult 0.223 
global conf_outbi 1.386 
global conf_outmult	0.405 
global conf_type 10_bconf

//SET-UP the simulation programs
do sim_setup_v2.do

///RUN SIM///
do sim_estimates_v2.do

