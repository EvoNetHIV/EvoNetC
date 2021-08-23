par(mfrow=c(3,2))
age_matches <- read.table("AgeMatches.txt")
diff_years <- unique(age_matches[,1])
for (ii in 1:length(diff_years)) {
  year_ii_rows <- which(age_matches[,1]==diff_years[ii])
  plot(age_matches[year_ii_rows,2],age_matches[year_ii_rows,3],xlim=c(15,80),ylim=c(15,80),xlab="Female Age",ylab="Male Age",
       main=paste("Age-based homophily at year",diff_years[ii]))
  abline(a=0,b=1)
}