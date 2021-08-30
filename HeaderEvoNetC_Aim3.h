/*
 *  NetworkHeader.h
 *
 *  Created by John Mittler on 12/6/12.
 *  Copyright 2012 University of Washington. All rights reserved.
 *  Treatment-as-prevention parameters added Winter of 2018
 *  Aim 3 (drug resistance) parameters added Winter of 2021
 */


// Files 
FILE *ParFile; // Contains parameters
FILE *ParOut;  // Output file that contains parameters only (used to verify that parameters were read in correctly) 
FILE *OutputFile; // Records average and standard deviation of viral loads (and other parameters) at each time step
FILE *NumPartnersFile; // Histogram num partners.  1st col = 0 partners, 2nd = 1 part, 21st col = 20+ partners. rows: total, women, men, concur-women, concur-men 
FILE *HeritOutputFile; // Records relationship between donors and recipients at a single timepoint (tprintHerit)
FILE *PatientRegistryOutputFile; // Records characteristics of each patient who died and state of living patients at tfinal
FILE *VLOutputFile; // Gives a histogram of log viral loads at different times
FILE *spVLOutputFile; // Gives a histogram of (log-transformed) set point viral loads at different times
FILE *NetworkStatsOutputFile; // Gives summaries of mean degree, etc. 
FILE *AgeDistFile; // Gives age distributions at various times during the simulation
FILE *RelLengthFile; // Gives average relationship lengths (uncensored) for people in different age groups
FILE *AgeMatchFile; // Gives list ages of all couples to test age-related homophily terms
FILE *RelationshipRegistryFile; // Gives list ages of all couples to test age-related homophily terms
FILE *GrandSummaryFile; // Gives final statistics from each run
FILE *AgentHistoryFile; // Detailed time course for selected agents
FILE *WarningsFile; // List of warning messages from each run
FILE *CheckSSCMaxFile; // Re-calculating the SSCMax for a couple of agents(i.e., patients) to verify drug concentration is correct

// Functions
void GetParameters(void), PrintStats(void);
void UserParams(void);
void InitializationRoutines(void);
void InitializeVariables(void);
void StartOfDayPrinting(void);
void EndOfDayPrinting(void);
void EndOfSimulationPrinting(void);
void OpenFiles(void);
void GetParameters(void);
void InitializePopulation(void);
void InitialAgeDistribution(void);
void CalculateLifeExpectancies(void);
void Births(void);
void Aging(void);
void UpdateDALYs(void);
void AddInfecteds(void);
void DefineSexualContactNetwork(void);
void PrintEdgeList(void);
void PrintHeaders(void);
void PrintStats(void);
void PrintHeritabilityStats(void);
void PrintPartnershipDistributions(void);
void InitiateTreatment(void);
void InitiateTreatmentGradual(void);
void TreatmentDropout(void);
void AdherenceModule(void);
void CollectStatistics(void);
void InitiateGradualTreatment(void);
void GetTreatmentProbs(double baseline_treat_prob);
void GetTreatmentProbsMultiplicative(double baseline_treat_prob);
void TreatAgents(int targeted_treatments, long num_newly_treated);
void GetRandomTreatmentProbs(double baseline_treat_prob);
void UpdateCD4Counts(void);
void UpdateCD4CountsStochastic(void);
void UpdateViralLoads(void);
void UpdateViralLoadsAim3(void);
void SimulateTransmission(void);
void BirthOfNewSusceptibles(void);
void NaturalDeath(void);
void DeathOfAIDSPatients(void);
void RemoveLinksToDeadPeople(void);
void SimulatePartnerShipDissolution(void);
void SimulateNewPartnershipFormation(void);
void CheckEdgeListForImpossibleValues(void);  // Does some elementary error checking
void RecordStatusOfLivingHIVPatientsForRegistryFile(void); 
void PrintFirstAndLastContacts(void);
void PrintAgeDistributions(void);
void PrintTrackedAgents(void);
void ConstantInputNewSusceptibles(void);
void VaccinatePopulation(void);
void SimulateVaccineDecay(void);
void UsualInitialValues(long i);
void PrintVLdistribution(void);
void AddNewLinks(long LinksToAdd);
void AddNewLinks2(long LinksToAdd);
void AddLinkForSinglePerson(long person1);
void PrintAgeMatches(void);
void PrintRelationshipRegistry(void);
void RandomizeMetabRate(void);

// Target treatment parameters
double PercentNonTargetedYear1;
double ActualNonTargetedYear1;
double baseline_treat_prob, total_prob;
double treatment_goal_start; 
double Start_Spontaneous_Treatment;
long repl;
long RandomIndex[5001];
long old_random_number_seed;
int first_run = 1; // used to ensure that big outputs are made only on first of many runs
int gradual_tx_incr = 0;  // If 1 the number of newly treated agents increases at a constant rate per year.  If 0, tx increases suddenly
                          // though the model allows for a yearly incr afterwards and a pre_TasP ramp-up.

long replicates = 1;
long N0, Infected0;

// Run options
int randomized_parameters = 0; // if 1 multiple parameters get randomized before the simulation

// Testing parameters
void HIVTesting(void);
int Diagnosed[5001];
long TimeDiag[5001];
double daily_prob_diagnosis = 1.0/365.0;
long delay_diag_tx = 14.0;
long delay_test_antibodies = 35.0;

// Treatment strategies
long under25 = 0, under30 = 0, under35 = 0, cd4_low = 0, cd4_super_low = 0, vl_high = 0, vl_very_high = 0, vl_super_high = 0, men = 0, women = 0, men27women23 = 0;
long men23women27 = 0, men33women27 = 0, men27women33 = 0;
double lowest_tx = 0.0, highest_tx = 1.0, tx_incr = 0.1;
double fold_incr_prob_with_targeting = 10.0;
double yearly_incr_tx = 0.01;
int TasP_Strategy = 1;

// Tracking variables for statistics
long susc_days = 0, susc_days_u50 = 0, susc_days_u25 = 0, susc_days_Mid = 0, susc_days_o50 = 0;
long new_infections = 0, new_infections_u50 = 0, new_infections_u25 = 0, new_infections_Mid = 0, new_infections_o50 = 0;;
long Time_All_Treated = 100000000, Last_Time_Not_All_Treated = 0;
long susc_days_last_5_years = 0, susc_days_last_5_years_u50 = 0, new_infections_last_5_years = 0, new_infections_last_5_years_u50 = 0;
long InfectedTreated = 0;
long pill_count, TotalTreated, NonTargetedTx, NonTargetedTxYear1, TotalTreatedYear1, TotalDiedAIDSAfterTasP=0;
long NumPartners[5001], NumPartnersLastTime[5001], Time_AIDS_Death[5001];
double ActTxStartTasP = 0.0; // Percent actually treated at the start of the TasP campaign (mainly useful for gradual treatment increase)
long TrackedAgentsGroup1start = 10000000, TrackedAgentsGroup1end = 10000000, TrackedAgentsGroup2start = 10000000, TrackedAgentsGroup2end = 10000000; 
long tx_person_days = 0, cd4_cat1_person_days = 0, cd4_cat2_person_days = 0, cd4_cat3_person_days = 0, cd4_cat4_person_days = 0, cd4_cat5_person_days = 0; // For DALY calcs

// DALY estimates
double DALY00, DALY01, DALY02, DALY03, DALY04, DALY05, DALY06, DALY07, DALY08, DALY10, DALY15, DALY20; // DALY03 means DALYs with 3% discount rate
double cost_hiv_treated = 0.003;
double cost_cd4_gt_350 = 0.02;
double cost_cd4_200_350 = 0.05;
double cost_cd4_lt_200 = 0.5;
double cost_died_AIDS = 1.0;

// Plotting options
int plt_GS=1;
int plt_NS= 1; 
int plt_out=0;
int plt_H=0;
int plt_R=0;
int plt_AH=0;
int plt_AM=1;
int plt_AD=1;
int plt_VL=0;
int plt_spVL=0;
int plt_warn=1;
long print_frequency = 1;

// Infection risk factors
double MaxAgeSex = 75.0;
double circum_prob = 0.4;
double RR_female_recipient = 1.5;
double RR_circum = 0.5;
double RR_youth = 2.0;
double RR_condom = 0.2;
double condom_prob16 = 0.69;
double age_condom_use_halves = 34.0; // Approx fit 2012 SA Nat Prev, Incid & Behav Surveyjj
double number_tx_spontaneous = 0.0;
double prob_sex_in_AIDS = 0.1;
double male_concurrency = 0.40;
double female_concurrency = 0.04;
double prob_conc_more_partners = 0.0;
double LongevityFudgeFactor = 1.0;
double increased_duration_per_year = 60.0; // In days.  If 50, duration will increase 500 days in 10 years, 2500 days (=6 years) in 50 years
double rate_stops_being_concurrent = 0.04 * 1.0/365; // 2.5% chance per year ==> ~50% drop over 25 years (as in Maugh-Brown et al. (2015 or 2015)  
double perc_tx_start = 0.80; 
double reduced_linkage_low_risk = 1.0;

// Aim3 variables

  long  Immune_Resp;
  long Chronic_Phase;
 
  double StochasticCut = 1e-6;
  double cut = 1e-7;
 
  double closeness_to_cutoff;
  double deltaDrug[5], deltaGIcnc[5];
  double ran_val;
  double base_h;
  double Time;
  double new_val, deltaX;
  double float_part;
  double plusterm, minusterm, probability;
  long count, problems;
  double r_base_orig; // r_base set by parent program.  Allows for diff r prior to peak.
  double  K_adj; // Adjusted carrying capacity that accounts for costs of drug resistance
  double  K_term; // Multiplication term that forces V to rapidly converge on K_adj
  double  deltaI[2][2][2][2][2];
  double  deltaM[2][2][2][2][2];
  double  deltaL[2][2][2][2][2];
  double  deltaV[2][2][2][2][2];
  double  r[2][2][2][2][2];
  double V_total;
  double I_total;
  double  deltaCD4;
  double fM = 0.019999, fL = 1.0e-6;
  double d = 0.6, dM = 0.04, dL = 0.001;
  
  long TherapyStarted, SecondLineTherapyStarted, Dosing_Interval;
  long Therapy_Type = 1;
  double start_drug = 100.0;

 
  double intrinsic_fitness[2][2][2][2][2];
  double Drug_effect;
  double GIcnc[5];
  double Drug[5001][5];
  double DrugDose[5];
  double BaseIC50[5];
  double  IC50[5][2][2][2][2][2];
  int i1,i2,i3,i4, i5, drug;
  double mu = 3e-5;
  double I_one_step = 0.0; // Sum of types that are one mutational step away
  double M_one_step = 0.0; // Sum of types that are one mutational step away
  double L_one_step = 0.0; // Sum of types that are one mutational step away
  double jrand;
  long Num_Loci = 3;
  int alleles[17];
  double Effect[6][5]; // Effect of mutation i on drug j (note array size one more than number to avoid subscripts of 0)
    
  double BaselineAdherence1, BaselineAdherence2, BaselineAdherence3, BaselineAdherence4, BaselineAdherence5;
  double stockout = 10000000.0;
  double restart = 10000000.0;
  double LifeExpectancy[101];

  double time_0;

  double DrugDose1 = 100.0, DrugDose2 = 100.0, DrugDose3 = 100.0, DrugDose4 = 100.0;

  double BaseIC50Drug1 = 100;
  double BaseIC50Drug2 = 100;
  double BaseIC50Drug3 = 100;
  double BaseIC50Drug4 = 100;

  double Adherence1 = 1.0;
  double Adherence2 = 1.0;
  double Adherence3 = 1.0;
  double Adherence4 = 1.0;
  double Adherence5 = 1.0;


  // Pharmacogenomic Parameters
  
  long TransDR[5001];// transmitted drug resistance status
  long AcqDR[5001];// acquired drug resistance status
  
  long metab_type[5001]; //agent attribute, metabolizer type 1=intermediate, 2=slow, 3=fast
  
  long mut_locus1[5001]; // Mutation at locus position 1, GenericTDF
  long mut_locus2[5001]; // Mutation at locus position 2, GenericEFV
  long mut_locus3[5001]; // Mutation at locus position 3, K103N
  long mut_locus4[5001]; // Mutation at locus position 4, M184V
  long mut_locus5[5001]; // Mutation at locus position 5, K65R
  
  double fast_metab_prop=0.0; //Proportion of fast metabolizers in the population
  double inter_metab_prop=1.0; //Proportion of intermediate metabolizers in the population
  double slow_metab_prop=0.0; //Proportion of slow metabolizers in the population
  
  long prop_mut_locus1; //Proportion of people with a mutant at locus position 1
  long prop_mut_locus2; //Proportion of people with a mutant at locus position 2
  long prop_mut_locus3; //Proportion of people with a mutant at locus position 3
  long prop_mut_locus4; //Proportion of people with a mutant at locus position 4
  long prop_mut_locus5; //Proportion of people with a mutant at locus position 5
  
  double fast_decay1_change=1.0; // value of decay change w/ fast metabolizers for drug1 
  double inter_decay1_change=1.0; // value of decay change w/ inter metabolizers for drug1
  double slow_decay1_change=1.0; // value of decay change w/ slow metabolizersfor drug1
  
  double fast_decay2_change=1.0; // value of decay change w/ fast metabolizers for drug2
  double inter_decay2_change=1.0; // value of decay change w/ inter metabolizers for drug2
  double slow_decay2_change=1.0; // value of decay change w/ slow metabolizers for drug2
  
  double fast_decay3_change=1.0; // value of decay change w/ fast metabolizers for drug3
  double inter_decay3_change=1.0; // value of decay change w/ inter metabolizers  for drug3
  double slow_decay3_change=1.0; // value of decay change w/ slow metabolizers for drug3
  
  double fast_decay4_change=1.0; // value of decay change w/ fast metabolizers for drug4
  double inter_decay4_change=1.0; // value of decay change w/ inter metabolizers  for drug4
  double slow_decay4_change=1.0; // value of decay change w/ slow metabolizers for drug4
  
  double cost1 = 0.01, cost2 = 0.01, cost3 = 0.02, cost4 = 0.01, cost5 = 0.0;

  double FC_D1_Mut1 = 10.0;
  double FC_D1_Mut2 = 1.0;
  double FC_D1_Mut3 = 1.0;
  double FC_D1_Mut4 = 1.0;
  double FC_D1_Mut5 = 1.0;
  
  double FC_D2_Mut1 = 1.0;
  double FC_D2_Mut2 = 10.0;
  double FC_D2_Mut3 = 1.0;
  double FC_D2_Mut4 = 1.0;
  double FC_D2_Mut5 = 1.0;
  
  double FC_D3_Mut1 = 1.0;
  double FC_D3_Mut2 = 1.0;
  double FC_D3_Mut3 = 10.0;
  double FC_D3_Mut4 = 1.0;
  double FC_D3_Mut5 = 1.0;
  
  double FC_D4_Mut1 = 1.0;
  double FC_D4_Mut2 = 1.0;
  double FC_D4_Mut3 = 1.0;
  double FC_D4_Mut4 = 1.0;
  double FC_D4_Mut5 = 1.0;


  long Interaction_Model_Drugs12;
  double cost_reduct5on1, cost_reduct4on2;
  double StopDrug1, RestartDrug1;
  double StopDrug2, RestartDrug2;
  double StopDrug3, RestartDrug3;
  double StopDrug4, RestartDrug4;
  double M_act;
  double L_act; 
  
  double drug_decay1;
  double drug_decay2;
  double drug_decay3;
  double drug_decay4;
  double drug_2nd_decay1;
  double drug_2nd_decay2;
  double drug_2nd_decay3;
  double drug_2nd_decay4;
  double conc_2nd_phase1;
  double conc_2nd_phase2;
  double conc_2nd_phase3;
  double conc_2nd_phase4;
  long additive_fitness;
 
  double growth_and_mut_terms_I; 
  int do_calcs; // Only do calculations when there is a possibility for growth or death


// Aim 3 functions
double define_growth_rate(long agent, long i1, long i2, long i3, long i4, long i5, int Interaction_Model_Drugs12);
void Initialize_values(void);
void init_alleles(long Num_Loci);
void drug_effects(void);
double GetAverageFitness(long additive_fitness, long i);
double UpdateCarryingCapacity(double K, double V_total, double V_peak, double CD4, double V_AIDS, 
                           double vl_increase_AIDS, double h, double prog_rate,
                           double SPVL, double Time_Inf, double t_peak, double t_acute);
double Get_I_one_step(long i, int i1, int i2, int i3, int i4, int i5);
double Get_L_one_step(long i, int i1, int i2, int i3, int i4, int i5);
double Get_M_one_step(long i, int i1, int i2, int i3, int i4, int i5);
void init_alleles(long Num_Loci);
void init_effect_array(void);
void initialize_values(void);
//void initialize_values(double h, double Drug1, double Drug2, double Drug3, double Drug4,
//                       double BaseIC50Drug1,double BaseIC50Drug2, double BaseIC50Drug3, double BaseIC50Drug4,
//                       double DrugDose1, double DrugDose2, double DrugDose3, double DrugDose4,
//                       double time_0, long Aim3RoundingErrors, double r_base, int additive_fitness);
void integerize_values(double cut, double StochasticCut, int i1, int i2, int i3, int i4);
double round_algo(double X,double dX, double cut, double StochasticCut);
void SetVtotal_ZeroOutDrugsAndDelta(void);
double StochasticRoutine(double X, double plusterm, double minusterm, double StochasticCut, double cut, double time_0, double h);
double adjust_carrying_capacity(long i, double K_in, double h);
void update_drug_concs_pill_combos(long agent, double h, double Time,double time_0,
                                   double drug_decay1, double drug_decay2, double drug_decay3, double drug_decay4,
                                   long TherapyStarted, long SecondLineTherapyStarted,
                                   double Adherence1,double Adherence2, double Adherence3, double Adherence4,
                                   long Therapy_Type, long Dosing_Interval);
void update_drug_concs_pill_combos_stockout(long agent, double h, double Time,double time_0,
                                   double drug_decay1, double drug_decay2, double drug_decay3, double drug_decay4,
                                   long TherapyStarted, long SecondLineTherapyStarted,
                                   double Adherence1,double Adherence2, double Adherence3, double Adherence4,
                                   long Therapy_Type, long Dosing_Interval);
void update_drug_concs_pill_combos_two_phase(long agent, double h, double Time,double time_0,
                                   double drug_decay1, double drug_decay2, double drug_decay3, double drug_decay4,
                                   double fast_decay1_change, double inter_decay1_change, 
                                   double slow_decay1_change, double fast_decay2_change, 
                                   double inter_decay2_change, double slow_decay2_change,
                                   double fast_decay3_change, double inter_decay3_change,
                                   double slow_decay3_change, double fast_decay4_change, 
                                   double inter_decay4_change, double slow_decay4_change,  
                                   long TherapyStarted, long SecondLineTherapyStarted,
                                   double Adherence1,double Adherence2, double Adherence3, double Adherence4,
                                   long Therapy_Type, long Dosing_Interval);

// General use functions
double norm_rand(double rn_mean, double rn_variance);
double jmpow(double a, double b);
double GetGammaDelay(double k, double theta);
long GetPoisson(long num, double pois_prob);
long max(long a, long b);
long rpois(double lamdba);

// Special case parameters
long Flat_Viral_Load = 0; // Setting this to 1 forces VL to be the same value during the entire primary infection period
double PrimaryInfectionLevel = 1.0e7; // Average viral load during primary infection
double expected_time_treated;
long Progression_Model;

// Printing parameters
long VL_print_time; // Instructs program to print out VL distributions every VL_print_time days.
double VL_print_lower, VL_print_upper, VL_print_interval;  // Lower and upper end of the distribution together with bin size
long VL_print_count = 0;
long printcount;
int printOutput;

// Network Parameters
double OrigAverageDegree;
long BreakUpList1[15001], BreakUpList2[15001], BreakUpList3[15001], BreakUpList4[15001]; // These arrays track indices of partnerships to be broken up
long RevisedNumLinks[5001], TimeLastBreakUp[5001];
long num_cases_no_link_found = 0;
int time_varying_mean_degree = 1;
double d_rate = 0.00025;
long NetStatsCount;

// New lines for vaccination
int Vaccinated[5001], ResistantToVaccine[5001];
long LengthVaccinated[5001];
double PercentVaccinated, VaccineDuration, PercentResistantToVaccine, RR_Vaccinated;

// Condom-campaign stuff
double MutationVariance, H, h;
long OrigN;
long Start_Condom_Campaign;
double Percent_using_condums_after_condom_campaign;
double orig_prob_sex;


/*Epidemiological parameters */
double Inf_Prob; // Probability of infection (calculated from quantities below)
double MaxInfRate, VHalfMaxInfRate, HillCoeffInfRate, beta; // Maximum probability of infection,  viral load at which infection prob is 0.5 of the max
double BaselineDeathRate, MaxAIDSDeathRate, VHalfMaxAIDSDeathRate, HillCoeffAIDSDeathRate, alpha; // Rate of inputs and deaths of people
double PopGrowth; // Rate that population increases each year
double yearly_rate, daily_rate; // rates of progression to AIDS using John's new formula
double prog_rate; // prog_rate is rate of increase of viral load per year
long Inf_Rate_Function; // 1 = power function from Lingappa et al, 0 = hill function  
double InfRateBaseline, InfRateExponent; // Prob(Infmission) = InfRateBaseline * [log(VL)] ^ InfRateExponent
double BirthRate; // Rate of entry of new susceptibles into the population.
double NaturalDeathRate; // Rate at which people normally die (of causes other than AIDS)
long MaxN = 10000; // NOTE: This is the extent of the remaining arrays minus 1
long MaxNReached = 0; // Flag to indicate if population was about to exceed the maximum array size
double tMaxNReached = 0.0; // Time that births or reincarnations were halted
long NewAIDSDeaths = 0; // Number of deaths in each time period (only used when reincarnation is allowed)
long Widowed_partner;
long newbirths; // Internal counter to count number of births each day
double Active; // Internal -- number of sexually active people
double asmr[101]; // Age-specific mortality rate.
 
/* Virological parameters */
double V[5001], I[5001], M[5001], L[5001];
double K[5001];
double V_vec[5001][2][2][2][2][2];
double I_vec[5001][2][2][2][2][2];
double M_vec[5001][2][2][2][2][2];
double L_vec[5001][2][2][2][2][2];
double Donors_V[5001], s[5001];
double V0, r0, V_peak, Max_VL_AIDS, VL_Incr_AIDS, t_peak, t_acute, Vss;
double d_acute[5001], Donors_d_acute[5001];
double death_rate_constant, death_rate_exponent; // Death rate parameters for daily probability death function (used only when Gamma_Death !=1)
double VL_Start_AIDS[5001]; 
long   Time_Start_AIDS[5001];
long   TimeInAIDS[5001];

double ExpectedDelayTime; // Assume gamma distributed delays
double Dmax, Dk, D50; // Parameters governing Fraser's gamma distributed times to AIDS
double shape_parameter,theta; // These parameters control the mean and variance of the gamma delay

/* Infection status parameters */
long N, Infected, Susceptible; // Initial number of people, Initial number infected, Number who could be infected
long currN;  // Value of N used to calculate the number of births that will occur
long Dead, DiedAIDS, DiedNat; // Cumulative Number Died of AIDS or natural causes
double prob_sex; // Probability that partners will have sex each day
int age_dependent_sex = 1; // Probablity of sex drops linearly until age 75 if set to 1
int age_dependent_mixing = 1; // People more likely to partner with someone their own age if set to 1
double Time_Inf[5001], Donors_Total_Time_Inf_At_Trans[5001]; // Time infected
int Status[5001]; // Status of each person (Infected=1, Suscp = 0, Died = -1, DiedAIDS = -2)
long Generation[5001], Donors_Generation[5001], Donors_Index[5001], NumRecipients[5001];
long HoneymoonDays;
double AgeAIDS[5001];
double IncreasedProbSexHoneymoonDays;

/* Risk group status */
int HadSexAlready[5001], HighRisk[5001];
double percent_high_risk;
long Sex[5001], Concurrent[5001], FSW[5001];
double Age[5001];
int Circumcised[5001];
 
/* Parameters governing heritability of VL setpoint */
double AverageLogSP0, VarianceLogSP0, Heritability; // Primary inputs
double SetPoint[5001], Donors_SetPoint[5001]; // Patient specific values (drawn from random number generator)
double LogSetPoint[5001], Donors_LogSetPoint[5001]; // Patient specific values (drawn from random number generator)
double ViralContribToLogSP0[5001],EnvirContribToLogSP0[5001]; // Viral and environmental deviations that contribute to final setpoint
double Donors_ViralContribToLogSP0[5001], Donors_EnvirContribToLogSP0[5001]; // Viral and environmental deviations that contribute to final setpoint

/* CD4 Dynamics */
int CD4[5001];
int CD4_nadir[5001];
int CD4_initial_value[5001];
double CD4time[5001];
int SpvlCat[5001];
int CD4_Exp_Flag = 1;
double CD4_TimeToAIDS[5001];
double CD4_TimeToAIDS_exp[5001][5];
int enhanced_progression_age = 0;
long CD4_treatment_delay_index[5001];
double prob_cd4_rebound, prob_cd4_rebound_further;

/* Look up table for CD4 counts */
double CD4_lookup[10][5];//cd4 as function of SPVL and time since infection
void CD4_Category();
int SPVL_Category(double SPVL_param);
double cd4_initProbs[10][4];
void cd4InitProbs();
void CreateCD4Table();
int initialCD4(int spvl_level);
int CD4_After_Treatment; // Viral load after treatment starts
int CD4_Determines_Treatment; //if 0, vl determines treament, if 1, then cd4

/* Statistical results (and internal params needed to get those statistics) */
double Vsum, Vave, Vstd, Vstdsum; // Statistics of viral load
double set_ave, set_sum, set_std, set_stdsum; // Statistics on logSetPoint
double d_ave, d_sum, d_std, d_stdsum; // Statistics on decay rate of infected cells after peak (our way of setting setpoints)
double G_ave, G_sum, G_std, G_stdsum; // Statistics on viral generation
double vcount; // Needed for averaging above quantities
double AveLinks,SumLinks;  // Average number of links per person

/* Treatment Parameters */
int first_tx_strategy = 1, last_tx_strategy = 1;  // Strategy 1 = random, strategy 2 = under age 25, ... stategy 9 = under age 25 and CD4 > 500
long UnderCare[5001]; // Pre-defined list of patients who would get treated after infection (e.g., they live near a clinic)
double ProbDropout[5001]; // Pre-defined list of patients who would get treated after infection (e.g., they live near a clinic)
long Treated[5001]; // Patients currently being treated
long TakingDrug[5001]; // Treated patients who are actually taking their drugs (i.e., they are adhering to their prescribed regimens)
long Time_Treated[5001]; // Day that patient started therapy
long TargetedTreated[5001]; // Patients currently being treated specifically b/c of a targeted TasP campaign
double treatment_prob[5001], rel_treatment_prob[5001], cum_treatment_prob[5001]; // Used in Initiate Treatment routine
long tx_days[5001]; // Number of days that each agent has been treated. (Useful for adjusting VL after a treatment interruption)
long tx_days_aids[5001]; // Number of days agent treated after onset of AIDS  (Useful for adjusting VL after a treatment interruption)
int ever_aids[5001]; // Agent ever progressed to AIDS (even with therapy)
double Start_TasP_Campaign; // Time (day of epidemic) when public health authorities start new treatment campaign
double prob_dropout; //  Per year probability of a treated person discontinuing therapy
double percent_under_care; // Percent of people who get treated after being infected
double Increment_Prob_Under_Care_High_CD4; // Increase in probability of going under care for people with low CD4 counts (= high CD4 categories)
double VL_After_Treatment; // Viral load after treatment starts
double Start_SexReduction_Campaign; // # Day of epidemic when public health authorities start treatment campaign
double Reduction_Mean_Degree; //  Reduction in mean degree resulting from this campaign
double Start_Faithfulness_Campaign; // Day of epidemic when public health authorities start campaign to reduce partnership turnover
double Increased_Duration; // Increase in partnership durations resulting from this campaign
double dropout_prob; // Per per probability of discontinuing therapy
/* Time-keeping variables */
long time, tfinal, tprintHerit;

/* Stopping paramters */
int StopEarly = 0;   // This variable is set to 1 when the program is unable to find a new link

/* Internal variables */
long i, j, k, daycount, deadcount;
double epsilon = 1e-10;
char junk[12],descript[30];

/* Variables required for random number generation*/
double nrand(void); // Normal random numbers with mean 0 and stdev 1
double randnum, normalrand; // internal variables
long random_number_seed; // This is the randum number seed
void init_genrand(unsigned long SSS); // All of the rest are existing MT functions
long genrand_int31(void);
unsigned long genrand_int32(void);
double genrand_real1(void);
double genrand_real2(void);
double genrand_real3(void);
double genrand_res53(void);

/* Sexual network parameters */
long NumLinks[5001], Links[5001][10]; // Number of sexual links for each person, Identities of those links (up to 5)
long TimeLinked[5001][10];
long TotalEdges; // Total number of sexual links in the entire population (note that one link links two individuals
double AverageDegree; // Average number of sexual partners person
long MaxLinks; // Maximum number of links per person (MaxLinks = 1 implies monogamy)
long NewLinks; // Temporary variable giving the number of links for each new person entering the population
long Link_back; // Internal variable to find links back from widow's of the dead person back to the dead person
double dp_TotalEdges; // double precision version of TotalEdges (for reading in only, immediately converted to long)
double Duration[5001]; // Length that person i tends to stay in a relationship
double MinDuration, MaxDuration; // Mins and Maxs used to draw a uniform random number for each person
double PartnershipDuration;  // Duration of each partnership (mean of the underlygin partnership duration of each partner)
double ProbPartnersBreakUp;  // Prob of partnerships breaking up 
long PartnersOfPerson1; // Internal parameter: index of each person's partners 
long LinksToRemove; // Temporary variable indicating number of links per dead person that should be purged
long AlreadyLinked; // Bookkeeping variable to prevent same persons from linking up more than once (but still counting this as two links)
long failed_to_find_a_new_link = 0; // Set to one if algorithm is unable to find a new pair (indicates some kind of programming problem)
double tr; // Turnover rate = Links_Broken_Per_Day/TotalEdges;
long escape_counter; // Used to prevent endless loops (in case of programming error)
long random_partner; // Which of person i's partners (if more than one) does person i have sex with today
long index_to_person2, found_link, index_to_person1; // Used internally to add and break links
long link_break; // Used internally to count up links that need to broken
long breakups; // Internal parameter: how many people broke up each day
double NLinks; // Number of people linked after deaths, births, and breakups
double dpLinksToAdd; // Number of links to add after accounting for deaths, births, and breakups
long LinksToAdd; // Number of links to add after accounting for deaths, births, and breakups
double ExpectedLinks; // Expected number of links that you would have had had their been now births, deaths, or breakups 
long HadDifficutlyFindingLinks = 0; // Flag set when the program had difficulty finding new links

