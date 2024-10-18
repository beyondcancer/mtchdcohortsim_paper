///SETUP///

clear

global workdir "/home/lsh1804992/simhrss"

cd $workdir

 //SET-UP the sim states	
set rng mt64s
set rngstream 47
set seed 88

//number of simulations
global nsim_start = 301 // where the job will start from
global nsim = 350 //number of simulations in the loop for the number of base cohorts to be created

 //SET-UP the simulation programs
qui do sim_setup.do

//SIM loop setup
global conf_test 1 2 3 4 5 6 //this is to set up the loop for different confounders
global ratio 1 3 5 10 //matching ratios 
global implementation noreplacement replacement //type of matching 

//SET-UP the DGM parameters
global n = 100000  //number of observations in each cohort
global loghr = 0 // true estimate
  
global lambda_exp .00003 .000015 .000021 .0000145 .000007 .0000097

global conf_expbi 0 0 0 0.39 0.39 0.39
global conf_outbi 0 0 0 0.10 0.10 0.10
global conf_expcont 0 0 0 0.0048 0.0048 0.0048
global conf_outcont 0 0 0 0.00029 0.00029 0.00029

///Confounding
global conf_type 10_noconf 5_noconf 7_noconf 10_smallconf 5_smallconf 7_smallconf

///RUN SIM///
do sim_estimates
