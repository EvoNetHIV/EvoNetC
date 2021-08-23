# Graph viral load, drug concentrations, mutants and CD4 T-cell categories #######
agent_data <- read.csv("AgentHistory.txt",sep="\t",header = TRUE)
inf_agents <- unique(agent_data$Agent)
repls <- unique(agent_data$repl)
for (jj in 1: length(repls)) {
  for (kk in 1:length(inf_agents)) {
    # Organize data from agent history file into a form that can be graphed
    agent_lines <- which(agent_data$Agent == inf_agents[kk] & agent_data$repl == repls[jj])
    if (length(agent_lines) >= 1) {
      time <- agent_data$time[agent_lines]
      vl <- agent_data$V[agent_lines]
      cd4 <- agent_data$CD4[agent_lines]
      D1 <- agent_data$D1[agent_lines]
      D2 <- agent_data$D2[agent_lines]
      D3 <- agent_data$D3[agent_lines]
      cd4 <- agent_data$CD4[agent_lines]
      WT <- agent_data$V00000[agent_lines]
      Single_muts <- agent_data$V10000[agent_lines] + agent_data$V01000[agent_lines] + agent_data$V00100[agent_lines] + agent_data$V00010[agent_lines] + agent_data$V00001[agent_lines] 
      Double_muts <-               agent_data$V11000[agent_lines] + agent_data$V10100[agent_lines] + agent_data$V10010[agent_lines] + agent_data$V10001[agent_lines]
      Double_muts <- Double_muts + agent_data$V01100[agent_lines] + agent_data$V01010[agent_lines] + agent_data$V01001[agent_lines] 
      Double_muts <- Double_muts + agent_data$V00110[agent_lines] + agent_data$V00101[agent_lines] + agent_data$V00011[agent_lines]
      Triple_muts <-               agent_data$V11100[agent_lines] + agent_data$V11010[agent_lines] + agent_data$V11001[agent_lines] 
      Triple_muts <- Triple_muts + agent_data$V01110[agent_lines] + agent_data$V01101[agent_lines] + agent_data$V00111[agent_lines]
      Quad_muts  <- agent_data$V11110[agent_lines] + agent_data$V11101[agent_lines] + agent_data$V11011[agent_lines] + agent_data$V10111[agent_lines] + agent_data$V01111[agent_lines] 
      Quint_muts <- agent_data$V11111[agent_lines]
      # Create the plots
      plot(time,log10(WT),ylim=c(-8,7),xlab="Time",ylab="Log Concentration",col="grey",type="l",
         main = paste("Agent ",inf_agents[kk]," (repl ",repls[jj],")",sep=""))
      lines(time,log10(vl),lwd=3,col="black")
      lines(time,log10(Single_muts),col="gold")
      lines(time,log10(Double_muts),col="orange")
      lines(time,log10(Triple_muts),lwd=2,col="red")
      lines(time,log10(Quad_muts),lwd=3,col="purple")
      lines(time,log10(Quint_muts),lwd=4,col="brown")
      lines(time,cd4,col="green")
      #lines(time,cd4,lwd=3,col="blue")
      lines(time,log10(D1),col="violet")
      lines(time,log10(D2),col="purple")
      lines(time,log10(D3),col="pink")
    }
  }
}
