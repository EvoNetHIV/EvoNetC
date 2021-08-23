# Create C header file called "UserParams.h" that is called by the EvoNetC.c
#   that updates user-specific parameters
# Also creates non-list versions of parameters in user_base_params and user_run_params
# Note: While takes an unconventional approach (R script that creates a C file), it solves a bunch of  
#       problems that cropped up in our older, more conventional approach.

# Create first line of the C header
cat("void UserParams(void) {",file="UserParams.h",sep="\n",append=FALSE)

# First add the user's base parameters
if (length(user_base_params) >= 1) {
  u_pars <- user_base_params
  for (kk in 1:length(user_base_params)) {
    cat(paste("  ",names(user_base_params[kk])," = ",u_pars[[kk]],";",sep=""),file="UserParams.h",sep="\n",append=TRUE)
    assign(names(user_base_params[kk]),u_pars[[kk]])
  }
}

# Second add any run time parameters (parameter that might be systematically varied inside a run loop)
if (length(user_run_params) >= 1) {
  u_pars2 <- user_run_params #unlist(user_run_params)
  for (kk in 1:length(user_run_params)) {
    cat(paste("  ",names(user_run_params[kk])," = ",u_pars2[[kk]],";",sep=""),file="UserParams.h",sep="\n",append=TRUE)
    assign(names(user_run_params[kk]),u_pars2[[kk]])
  }
}

#cat("printf(\"end of user params\");",file="UserParams.h",sep="\n",append=TRUE)
cat("}",file="UserParams.h",sep="\n",append=TRUE)

#
Sys.sleep(1)
