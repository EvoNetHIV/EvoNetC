setwd("~/EvonetC")

options(warn = 2)
source("ProcessDefaultParams.R") # Reads in default parameters from a text file and puts them into R format 

# User-specific parameters.  These will over-write any default parameters with the same name.
user_base_params <- c(
  tfinal = 365+100, ##The duration of the run
  N0 = 100,
  Num_Loci = 5,
  replicates = 1,
  Infected0 = 50,
  fast_metab_prop = 0.60,
  inter_metab_prop = 0.10,
  slow_metab_prop = 0.30,
  fast_decay3_change=1.0, # value of decay change for drug3
  inter_decay3_change=1.0, # value of decay change for drug3
  slow_decay3_change=1.0, # value of decay change for drug3
  DrugDose1 = 50.0,
  DrugDose2 = 50.0,
  DrugDose3 = 50.0,
  BaseIC50Drug1 = 100,  BaseIC50Drug2 = 100,  BaseIC50Drug3 = 100, BaseIC50Drug4 = 100,
  FC_D1_Mut1 = 200.0, FC_D1_Mut2 =   1.0, FC_D1_Mut3 =   1.0, FC_D1_Mut4 =  10.0, FC_D1_Mut5 =   1.0, 
  FC_D2_Mut1 =   1.0, FC_D2_Mut2 = 300.0, FC_D2_Mut3 =   1.0, FC_D2_Mut4 =   2.0, FC_D2_Mut5 =  20.0, 
  FC_D3_Mut1 =   1.0, FC_D3_Mut2 =   1.0, FC_D3_Mut3 = 400.0, FC_D3_Mut4 =   1.0, FC_D3_Mut5 =   1.0, 
  FC_D4_Mut1 =   1.0, FC_D4_Mut2 =   1.0, FC_D4_Mut3 =   1.0, FC_D4_Mut4 =   1.0, FC_D4_Mut5 =   1.0, 
  drug_decay1 = 0.1,  drug_decay2 = 0.1,  drug_decay3 = 0.1,
  Adherence1 = 1.0,  Adherence2 = 1.0,  Adherence3 = 1.0,
  Start_Spontaneous_Treatment = 100,
  Start_TasP_Campaign =100,
#  stockout = 11250,
#  restart = 11300,
  AverageLogSP0 = 4.0,
  VarianceLogSP0 = 1.0,
  plt_AD = 0,
  plt_NS = 1,
  plt_AM = 0,
  mu = 1e-4,
  print_frequency = 10
)

# Run program with different drug decay rates
par(mfrow=c(2,2))

V0_vec<-c(0.001)

for (ii in 1: length(V0_vec)) {
  user_run_params <- list(V0 = V0_vec[ii]) # Run-specific parameters
 
  #Add loops for task 1 - Aim 2
  
  source("UpdateParams.R") # Create file with user_base_params and user_run_params
  # Compile program: CSDE server C:/Rtools ... / mingw64...gcc,  MacOS: system("gcc -O3 EvoNetC.c")
  if (system("C:/Rtools40/mingw64/bin/gcc EvoNetC_Aim3.c") != 0) {stop("Halting program due to complier error\n!")}
  #if (system("gcc -O3 EvoNetC_Aim3.c -o a.exe") != 0) {stop("Halting program due to complier error\n!")}
  system("./a.exe")
  if (plt_GS == 1) source("Plot_Time_Series.R")
  if (plt_NS == 1) source("Plot_Network_Stats.R")
  if (plt_AM == 1) source("Plot_Age_Matches.R")
  if (plt_AD == 1) source("Plot_Age_Distribution.R")
  source("Plot_Percent_Resistant.R")
  
}

##export AgentHistory.txt file
agenthistory <- read.delim("AgentHistory.txt")
write.csv(agenthistory,"agenthistory_8192021_1129_fastdecayhigh_slowdecaylow__adh100_2reps_n500_inf175_1year.csv", row.names=TRUE)

##export CheckSSCMax.txt file
checksscmax <- read.delim("CheckSSCMax.txt")
write.csv(checksscmax,"CheckSSCMax_8192021_1129_fastdecayhigh_slowdecaylow__adh100_2reps_n500_inf175_1year.csv", row.names=TRUE)
source("Plotting_SSCMax_VerifyOutput.R")
