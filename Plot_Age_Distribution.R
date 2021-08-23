par(mfrow=c(3,2))
age_dist <- read.table("AgeDistribution.txt",header=TRUE,check.names=FALSE)
total <- 0
for (ii in 1:nrow(age_dist)) {
  total <- sum(age_dist[ii,])
  total_over65 <- sum(age_dist[ii,(65-16):(100-16)])
  plot(as.numeric(colnames(age_dist)),age_dist[ii,],xlab="Age",ylab="Number",ylim=c(0,max(age_dist)),
       main=paste("Age distribution at year ",rownames(age_dist)[ii]," (",round(100*total_over65/total,1),"% over 65)",sep=""))
}