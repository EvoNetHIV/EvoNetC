options(error=NULL)

parameters <- read.table("DefaultParameters.txt",header <- FALSE, row.names <- NULL)

ns_with_burn_in <- read.table("NetworkStats.txt",header=TRUE)
#ns <- ns_with_burn_in[11:nrow(ns_with_burn_in),]
ns <- ns_with_burn_in[11:nrow(ns_with_burn_in),]
Prev50 <- ns[,"Prev50"]
Time <- ns[,"Time"]/365
Time <- Time - Time[1] 
final_time <- length(Time)
Sus <- ns[,"Sus"] ; HIV_pos <- ns[,"Pos"] ;Alive <- ns[,"Alv"]
Age <- ns[,"Age"] ; AgeSusc <- ns[,"AgeS"]; AgeInf <- ns[,"AgeI"]
Females <- ns[,"Fem"]; f_u25 <- ns[,"f25"]; f_Mid <- ns[,"fMid"]; f_o50 <- ns[,"f50"] 
Males <- ns[,"Male"];  m_u25 <- ns[,"m25"]; m_Mid <- ns[,"mMid"]; m_o50 <- ns[,"m50."] 
Inf_f <- ns[,"if."]; Inf_f_u25 <- ns[,"if25"]; Inf_f_Mid <- ns[,"ifMid"]; Inf_f_o50 <- ns[,"if50"]
Inf_m <- ns[,"im"];  Inf_m_u25 <- ns[,"im25"]; Inf_m_Mid <- ns[,"imMid"]; Inf_m_o50 <- ns[,"im50"]

Inf_u25 <- ns[,"i25"]; Inf_Mid <- ns[,"iMid"]; Inf_o50 <- ns[,"i50"] 
Inf_Tx_u25 <- ns[,"Inf_Tx_u25"]; Inf_Tx_Mid <- ns[,"Inf_Tx_Mid"]; Inf_Tx50 <- ns[,"Inf_Tx50"]
Inf_Tx_f <- ns[,"Inf_Tx_f"]; Inf_Tx_m <- ns[,"Inf_Tx_m"]

md <- ns[,"md"] ;   md_m <- ns[,"md_m"] ; md_f <- ns[,"md_f"]
conc <- ns[,"conc"]; c_m <- ns[,"c_m"]  ; c_f <- ns[,"c_f"]
md_u25 <- ns[,"md_u25"] ; md_Mid <- ns[,"mdMid"] ;md50 <- ns[,"md50"]
c_u25 <- ns[,"c_u25"] ; c_Mid <- ns[,"c2550"] ;c50 <- ns[,"c50"]
Tx <- ns[,"Tx"] ; TarTx <-  ns[,"TarTx"]; dAIDS <- ns[,"dAIDS"]; dNat <- ns[,"dNat"]

TotTx <- ns[,"TotTx"] ; NonTarg <- ns[,"NoTarg"]

MeanDur <- ns[,"Dur"]/365; MeanDur25 <- ns[,"Dur25"]/365; MeanDurMid <- ns[,"DurMid"]/365
MeanDur50 <- ns[,"Dur50"]/365; MeanDur18 <- ns[,"Dur18"]/365

susc_d <- ns[,"susc_d"]; new_inf <- ns[,"new_inf"]
new_infs_last_year <- new_inf[2:length(new_inf)] - new_inf[1:(length(new_inf)-1)]
incid <- 100*365.0*new_infs_last_year/(susc_d[2:length(susc_d)] - susc_d[1:(length(susc_d)-1)])

susc_d_u50 <- ns[,"susc_d_u50"]; new_inf50 <- ns[,"new_inf50"]

new_infs_last_year50 <- new_inf50[2:length(new_inf50)] - new_inf50[1:(length(new_inf50)-1)]
incid50 <- 100*365.0*new_infs_last_year50/(susc_d_u50[2:length(susc_d_u50)] - susc_d_u50[1:(length(susc_d_u50)-1)])

susc_d25 <- ns[,"susc_d25"]; new_inf25 <- ns[,"newinfs25"]
new_infs_last_year25 <- new_inf25[2:length(new_inf25)] - new_inf25[1:(length(new_inf25)-1)]
incid25 <- 100*365.0*new_infs_last_year25/(susc_d25[2:length(susc_d25)] - susc_d25[1:(length(susc_d25)-1)])

susc_dMid <- ns[,"susc_dMid"]; new_infMid <- ns[,"newinfsMid"]
new_infs_last_yearMid <- new_infMid[2:length(new_infMid)] - new_infMid[1:(length(new_infMid)-1)]
incidMid <- 100*365.0*new_infs_last_yearMid/(susc_dMid[2:length(susc_dMid)] - susc_dMid[1:(length(susc_dMid)-1)])

susc_d_o50 <- ns[,"susc_d_o50"]; new_inf_o50 <- ns[,"newinfs_o50"]
new_infs_last_year_o50 <- new_inf_o50[2:length(new_inf_o50)] - new_inf_o50          [1:(length(new_inf_o50)-1)]
incid_o50 <- 100*365.0*new_infs_last_year_o50/(susc_d_o50[2:length(susc_d_o50)] - susc_d_o50[1:(length(susc_d_o50)-1)])

par(mfrow=c(3,2))
plot(Time,HIV_pos/Alive,col="black",type="l", ylab = "Prevalence",
     ylim=c(0,max(Inf_f/Females,Inf_m/Males,HIV_pos/Alive,Prev50)),
     main = "Prevalence [Total: black, Males: blue, Females: red, Under 50: grey]",lwd=4)
lines(Time,Inf_f/Females,col="red",type="l",lwd=2)
lines(Time,Inf_m/Males,col="blue",type="l",lwd=2)
lines(Time,Prev50,col="gray",type="l",lwd=3)

plot(Time[-1],incid,ylim=c(0,max(1.1*incid50,1.1*incid25)),lwd=2,xlab="Time (years)",ylab= "Incidence rate / 100 person-years", 
     main="Incidence [Total: black, Under 50: blue, Under 25: red, 25-49: orange, Over 50: yellow]",type="l")
lines(Time[-1],incid50,col="blue",lwd=2)
lines(Time[-1],incid25,col="red",lwd=2)
lines(Time[-1],incidMid,col="orange",lwd=2)
lines(Time[-1],incid_o50,col="yellow",lwd=2)

plot(Time,HIV_pos,ylim=c(0,max(HIV_pos)),col="red",type="l", ylab = "Number", main = "Number of people [Infected: red, Treated: blue, Targeted Tx: purple]",lwd=2)
lines(Time,Tx,col="blue",lwd=2)
lines(Time,TarTx,col="purple",lwd=2)


non_zero_incid <- which(incid50 > 0)
prev_incid_ratio <- Prev50[non_zero_incid]/incid50[non_zero_incid]
treat_ratio <- Tx/HIV_pos
plot(Time,50*treat_ratio,type="l",lwd=0.5,ylab="Prevalence / Incidence (gry = %tx)",
     main="Prevalence / Incidence")
lines(Time[1:28],100*Prev50[1:28]/incid50[1:28],type="l",lwd=2,col="red")

plot(Time,dAIDS,ylim=c(0,max(dAIDS,dNat)),col="red",type="l", ylab = "Number", main = "Number of deaths [AIDS: red, Natural: blue]",lwd=2)
lines(Time,dNat,col="blue",lwd=2)

plot(Time,Alive,ylim=c(0,max(Alive)),type="l", ylab = "Number", lwd=4,
     main = "Number of people  [Total: black, Males: blue, Females:red]",col="black")
lines(Time,Females,col="red",lwd=2)
lines(Time,Males,col="blue",lwd=2)

plot(Time,Inf_f_u25/f_u25,col="red",type="l", ylab = "Prevalence", ylim = c(0,max(Inf_f_u25/f_u25)),
     main = "Prevalence in people under 25 [Females : red, males: blue]",lwd=2)
lines(Time,Inf_m_u25/m_u25,col="blue",lwd=2)

plot(Time,Inf_f_Mid/f_Mid,col="red",type="l", ylab = "Prevalence", ylim = c(0,max(Inf_f_Mid/f_Mid)),
     main = "Prevalence in people between 25 and 50 [Females : red, males: blue]",lwd=2)
lines(Time,Inf_m_Mid/m_Mid,col="blue",lwd=2)

plot(Time,Inf_f_o50/f_o50,col="red",type="l", ylab = "Prevalence", ylim = c(0,max(Inf_f_o50/f_o50,Inf_m_o50/m_o50)),
     main = "Prevalence in people over 50 [Females : red, males: blue]",lwd=2)
lines(Time,Inf_m_o50/m_o50,col="blue",lwd=2)

plot(Time,md,ylim=c(0.5*min(md,md_m,md_f),1.5*max(md,md_m,md_f)),type="l", lwd=4, ylab = "Mean Degree",
     main="Mean Degree [Total: black, Males: blue, Females: red]", col="black")
lines(Time,md_m,col="blue",lwd=2)
lines(Time,md_f,col="red",lwd=2)

plot(Time,md_u25,ylim=c(0,max(md_u25,md_Mid,md50)),type="l",lwd=2, ylab = "Mean Degree", main="Mean Degree[<25: black, 25-50: blue, >50: green]",col="black")
lines(Time,md_Mid,col="blue",lwd=2)
lines(Time,md50,col="green",lwd=2)

partner_dist <- read.table("NumPartners.txt",header=FALSE)
num_partners <- c(0,1,2,3,4,5,6,7,8,9)
plot(num_partners,100*partner_dist[2,1:10]/sum(partner_dist[2,1:10]),type="b",lwd=2,col="red", ylim=c(0,100),
     xlab="Number of partners",ylab="Percent",
     main="Partner Distribution [Females: red, Males: blue]")
lines(num_partners,100*partner_dist[3,1:10]/sum(partner_dist[3,1:10]),type="b",lwd=2,col="blue")

plot(Time,MeanDur,col="black",type="l", ylab = "Relationship Duration (years)", ylim = c(0, max(MeanDur50)),
     main = "Duration tendancies [Total: black, <18: purple, <25: red, 25-50: orange, >50: yellow]",lwd=4)
lines(Time,MeanDur18,col="purple",lwd=2)
lines(Time,MeanDur25,col="red",lwd=2)
lines(Time,MeanDurMid,col="orange",lwd=2)
lines(Time,MeanDur50,col="yellow",lwd=2)

plot(Time,Inf_u25,ylim=c(0,max(Inf_u25)),col="red",type="l", ylab = "Number", main = "Treatment: under age 25 [Infected: red, Treated: blue]",lwd=2)
lines(Time,Inf_Tx_u25,col="blue",lwd=2)

plot(Time,Inf_Mid,ylim=c(0,max(Inf_Mid)),col="red",type="l", ylab = "Number", main = "Treatment: Ages 25-50 [Infected: red, Treated: blue]",lwd=2)
lines(Time,Inf_Tx_Mid,col="blue",lwd=2)

plot(Time,Inf_o50,ylim=c(0,max(Inf_o50)),col="red",type="l", ylab = "Number", main = "Treatment: Over age 50 [Infected: red, Treated: blue]",lwd=2)
lines(Time,Inf_Tx50,col="blue",lwd=2)

plot(Time,Inf_f,ylim=c(0,max(Inf_f)),col="red",type="l", ylab = "Number", main = "Treatment: Females [Infected: red, Treated: blue]",lwd=2)
lines(Time,Inf_Tx_f,col="blue",lwd=2)

plot(Time,Inf_m,ylim=c(0,max(Inf_m)),col="red",type="l", ylab = "Number", main = "Treatment: Males [Infected: red, Treated: blue]",lwd=2)
lines(Time,Inf_Tx_m,col="blue",lwd=2)

plot(Time,c_u25,ylim=c(0,max(c_u25,c_Mid,c50)),type="l",lwd=2, ylab = "Concurrency", main="Concurrency [<25: black, 25-50: blue, >50: green]",col="black")
lines(Time,c_Mid,col="blue",lwd=2)
lines(Time,c50,col="green",lwd=2)

plot(Time,conc,ylim=c(0,max(conc,c_m,c_f)),type="l",ylab="Concurrency", 
     main="Concurrency [Total: black, Males: blue, Females: red]",col="black",lwd=4)
lines(Time,c_m,col="blue",lwd=2)
lines(Time,c_f,col="red",lwd=2)


AgeRange <- c(20,37.5,55)
MP_m25 <- ns[,"MP_m25"]; MP_mMid <- ns[,"MP_mMid"]; MP_m50 <- ns[,"MP_m50"];
MP_f25 <- ns[,"MP_f25"]; MP_fMid <- ns[,"MP_fMid"]; MP_f50 <- ns[,"MP_f50"];
sTime <- 20
Male_MP <- 100*c(MP_m25[sTime], MP_mMid[sTime], MP_m50[sTime])
Female_MP <- 100*c(MP_f25[sTime], MP_fMid[sTime], MP_f50[sTime])
plot(AgeRange,Male_MP,col="blue",lwd=2,type="b",xlab = "Age", ylab = "Percent Multiple Partners",ylim=c(0,max(40,Female_MP,Male_MP)),
     main = paste(">1 Partner previous year  (year ",sTime,")  [Men: blue, Females: red]\n(Thin lines SA HIV survey 2012, Fig 3.8)",sep=""))
lines(AgeRange,Female_MP,col="red",lwd=2,type="b")
lines(AgeRange,c(37.5,18.3,6.5),col="blue",type="l",lwd=0.5)
#lines(AgeRange,c(30.8,14.8,3.7),col="blue",type="l",lwd=0.5)
lines(AgeRange,c(8.2,4.0,0.8),col="red",type="l",lwd=0.5)
#lines(AgeRange,c(6,3.0,0.8),col="red",type="l",lwd=0.5) 

AvePartMen20 <- ns[,"AvePartMen20"]; AvePartMen25 <- ns[,"AvePartMen25"]; AvePartMen30 <- ns[,"AvePartMen30"]; 
AvePartMen40 <- ns[,"AvePartMen40"]; AvePartMen50 <- ns[,"AvePartMen50"]; 
AvePartWomen20 <- ns[,"AvePartWomen20"]; AvePartWomen25 <- ns[,"AvePartWomen25"]; AvePartWomen30 <- ns[,"AvePartWomen30"]; 
AvePartWomen40 <- ns[,"AvePartWomen40"]; AvePartWomen50 <- ns[,"AvePartWomen50"]; 
AgeRange <- c(15, 17.5,22.5,27.5,45,55)
sTime <- 40
if (nrow(ns) < sTime) sTime <- nrow(ns)
MaleLifetime <-  c(0, AvePartMen20[sTime],AvePartMen25[sTime],AvePartMen30[sTime],AvePartMen40[sTime],AvePartMen50[sTime])
FemaleLifetime <-  c(0, AvePartWomen20[sTime],AvePartWomen25[sTime],AvePartWomen30[sTime],AvePartWomen40[sTime],AvePartWomen50[sTime])
plot(AgeRange,MaleLifetime,col="blue",lwd=2,type="b",xlab = "Age", ylab = "Cumulative Number of Partners",
     ylim=c(0,max(MaleLifetime,FemaleLifetime,20)),
     main = paste("Number of partners end of simulation (year ",sTime,") [Men: blue, Females: red]\n(Thin lines SA Demo survey 2016)",sep=""))
lines(AgeRange,FemaleLifetime,col="red",lwd=2,type="b")
lines(AgeRange,c(0,5.7,11,15,19.3,17.1),col="blue",lwd=0.5,type="l")
lines(AgeRange,c(0,1.9,3.4,4.0,4.8,3.7),col="red",lwd=0.5,type="l")

plot(Time,Age,col="black",lwd=2,type="b",xlab = "Time (years)", ylab = "Age",
     main = "Average Age [Total: black, Susc: blue, Infected: red")
plot(Time,Age,col="black",lwd=2,type="b",xlab = "Time (years)", ylab = "Age",
     ylim=c(0,max(Age,AgeInf,AgeSusc)),
     main = "Average Age [Total: black, Susc: blue, Infected: red")
lines(Time,AgeInf,col="red",lwd=2,type="l")
lines(Time,AgeSusc,col="blue",lwd=2,type="l")



if (parameters[which(parameters[,2]=="plt_AD"),1] == 1) {
  par(mfrow=c(3,2))
  age_dist <- read.table("AgeDistribution.txt",header=TRUE,check.names=FALSE)
  total <- 0
  for (ii in 1:nrow(age_dist)) {
    total <- sum(age_dist[ii,])
    total_over65 <- sum(age_dist[ii,(65-16):(100-16)])
    plot(as.numeric(colnames(age_dist)),age_dist[ii,],xlab="Age",ylab="Number",ylim=c(0,max(age_dist)),
      main=paste("Age distribution at year ",rownames(age_dist)[ii]," (",round(100*total_over65/total,1),"% over 65)",sep=""))
  }
}

if (parameters[which(parameters[,2]=="plt_AM"),1] == 1) {
  par(mfrow=c(3,2))
  age_matches <- read.table("AgeMatches.txt")
  diff_years <- unique(age_matches[,1])
  for (ii in 1:length(diff_years)) {
    year_ii_rows <- which(age_matches[,1]==diff_years[ii])
    plot(age_matches[year_ii_rows,2],age_matches[year_ii_rows,3],xlim=c(15,80),ylim=c(15,80),xlab="Female Age",ylab="Male Age",
       main=paste("Age-based homophily at year",diff_years[ii]))
    abline(a=0,b=1)
  }
}

if (1==2) {
 if (parameters[which(parameters[,2]=="plt_R"),1] == 1 & length(TimeToAIDS) >=1) {
  registry <- read.table("PatientRegistryOutput.txt",row.names=NULL,header=TRUE)
  par(mfrow=c(2,2))
  Time <- registry[,"Time"]/365 
  SPVL <- registry[,"LogSP0"]
  Age <- registry[,"Age"] 
  TimeInf <- registry[,"TimeInf"]/365 
  Status <- registry[,"Status"]
  DiedAIDS <- which(Status == "DiedAIDS")
  AgeInfected <- Age - (Time - TimeInf)
  DiedAIDS_inf_5_20 <- which(Status == "DiedAIDS" & TimeInf > 2 & TimeInf < 20)
  TimeToAIDS = Time[DiedAIDS_inf_5_20]- TimeInf[DiedAIDS_inf_5_20]
  plot(TimeInf[DiedAIDS_inf_5_20],TimeToAIDS,xlab="Time Infected",ylab="Time To AIDS")
  plot(SPVL[DiedAIDS_inf_5_20],TimeToAIDS,xlab="SPVL",ylab="Time to AIDS",main=paste("Time To AIDS (median:",summary(TimeToAIDS)[3],", mean: ",summary(TimeToAIDS)[4],")",sep=""))
  plot(AgeInfected[DiedAIDS_inf_5_20],TimeToAIDS,xlab="Age Infected",ylab="Time To AIDS")
  lines(c(16,31),c(19,4))
  hist(TimeToAIDS,xlab="Time To AIDS")
  print(summary(TimeToAIDS))
 }
}
# Graph individual agent trajectories
if (1 == 2 & parameters[which(parameters[,2]=="plt_AH"),1] == 1) {
 Group1start = parameters[which(parameters[,2]=="TrackedAgentsGroup1start"),1]
 Group1end = parameters[which(parameters[,2]=="TrackedAgentsGroup1end"),1]
 Group2start = parameters[which(parameters[,2]=="TrackedAgentsGroup2start"),1]
 Group2end = parameters[which(parameters[,2]=="TrackedAgentsGroup2end"),1]

 if ( (Group1end - Group1start >0 ) || (Group2end - Group2start > 0)) {
  ah <- read.table("AgentHistory.txt",row.names=NULL,header=TRUE)
  ah_time <- ah[,"time"]/365
  ah_agent <- ah[,"Agent"]
  ah_spvl <- ah[,"SPVL"]
  ah_drop <- ah[,"Pdrop"]
  ah_age <- ah[,"Age"]
  ah_tInf <- ah[,"tInf"]/365
  ah_tx <- ah[,"tx"]
  ah_V <- ah[,"V"]
  ah_cd4 <- ah[,"CD4"]
  unique_agents <- unique(ah_agent)
  
  par(mfrow=c(3,2))
  for (ii in 1:length(unique_agents)) {
     ii_rows <- which(ah[,"Agent"] == unique_agents[ii]) # ii_rows is a list of all rows within AgentHistory for agent ii
     plot(ah_time[ii_rows],log10(ah_V[ii_rows]),col="red",type="l",lwd=2, xlab = "Time (Years)", ylab="VL (red), CD4 (blue)",
          main = paste("Treatment history for agent ",unique_agents[ii],
                        " [ Infected at year ",round(ah_tInf[ii_rows][1],1),
                        " when ",round(ah_age[ii_rows][1]),
                       " yrs old, P_drop = ",
                       round(100*365*ah_drop[ii_rows][1],1),"% ]",
                       sep=""),
          ylim=c(0,7))
     lines(ah_time[ii_rows],ah_cd4[ii_rows],type="l",lwd=2,col="blue")
     lines(ah_time[ii_rows],log10(ah_spvl[ii_rows]),type="l",lwd= 0.4,col="grey")
     lines(ah_time[ii_rows],4*ah_tx[ii_rows]-1,type="p",lwd="3",pch=15,col="purple")
   }
 }
}




