# Convert parameter data into R format for subsequent manipulation in R
params <- read.table("DefaultParameters.txt",header=FALSE,colClasses = c("numeric","character"))
param_vals <- as.numeric(as.character(params[,1]))
param_names <- params[,2]

# Create a dataframe that includes all of the parameters for anyone who wants to see the complete list in one place
param.df <- data.frame(param_names, param_vals)

# Create R variables that have the same name and value as the 
for (ii in 1:nrow(params)) {
  assign(param_names[ii],param_vals[ii]) 
}



