# Graph viral load, drug concentrations, mutants and CD4 T-cell categories #######
par(mfrow=c(1,1))
sample_times <- unique(agent_data$time)
graph_threshold_highly_resistant <- 0.1
graph_threshold_partially_resistant <- 0.1
for (jj in 1: length(repls)) {
  par(mfrow=c(1,1))
  for (kk in 1:length(sample_times)) {
    # Organize data from agent history file into a form that can be graphed
    sample_lines <- which(agent_data$time == sample_times[kk] & agent_data$repl == jj)
    Single_muts <-               agent_data$V10000[sample_lines] + agent_data$V01000[sample_lines] + agent_data$V00100[sample_lines] + agent_data$V00010[sample_lines] + agent_data$V00001[sample_lines] 
    Double_muts <-               agent_data$V11000[sample_lines] + agent_data$V10100[sample_lines] + agent_data$V10010[sample_lines] + agent_data$V10001[sample_lines]
    Double_muts <- Double_muts + agent_data$V01100[sample_lines] + agent_data$V01010[sample_lines] + agent_data$V01001[sample_lines] 
    Double_muts <- Double_muts + agent_data$V00110[sample_lines] + agent_data$V00101[sample_lines] + agent_data$V00011[sample_lines]
    Triple_muts <-               agent_data$V11100[sample_lines] + agent_data$V11010[sample_lines] + agent_data$V11001[sample_lines] 
    Triple_muts <- Triple_muts + agent_data$V01110[sample_lines] + agent_data$V01101[sample_lines] + agent_data$V00111[sample_lines]
    Quad_muts  <- agent_data$V11110[sample_lines] + agent_data$V11101[sample_lines] + agent_data$V11011[sample_lines] + agent_data$V10111[sample_lines] + agent_data$V01111[sample_lines] 
    Quint_muts <- agent_data$V11111[sample_lines]
    resistant <- length(which(Triple_muts + Quad_muts + Quint_muts >= graph_threshold_highly_resistant))
    partially_resistant <- length(which(Double_muts + Triple_muts + Quad_muts + Quint_muts >= graph_threshold_partially_resistant))
    total <- length(sample_lines)
    if (kk == 1) {
      plot(sample_times[kk],resistant/total,xlim=c(0,max(sample_times)),
           main = paste("Average Resistance in the Population (repl ",jj,")"), 
           ylim = c(0,1),xlab="Time",
           ylab=paste("Percent agents with at least ",100*graph_threshold_highly_resistant,"% resistant virus",sep=""),col="red",cex=1.5,type="p")
    } else {
      lines(sample_times[kk],resistant/total,xlab="Time",ylab="Percent Resistant",col="red",type="p",cex=1.5)
    }
    lines(sample_times[kk],partially_resistant/total,xlab="Time",ylab="Percent Resistant",col="orange",cex=0.7,pch=8,type="p")
  }
}