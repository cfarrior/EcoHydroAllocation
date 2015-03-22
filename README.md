# EcoHydroAllocation

#to use follow these steps: 
#(1) Create two folders: "~/Documents/Farrior_etal_EcoHydroAllocation/" and "~/Documents/Farrior_etal_EcoHydroAllocation/Figures"
#(2) Save this text as a text document, named "~/Documents/Farrior_etal_EcoHydroAllocation/EcoHydroAllocation.r"
#(3) Copy and paste the following into R (without the first hashmarks)
#   setwd("~/Documents/Farrior_etal_EcoHydroAllocation/")
#	  source("EcoHydroAllocation.r")
#(4) Change any parameter values you want to, by writing their desired values into R (i.e. "t = 180" changes the growing season length to 180 days). 
#   See below and the paper for defaults.
#(5) Pick a name for the set of runs; here "runName". Running the following functions will find the ESSs and then produce the plots from the paper. 
#   FindESSs("runName")  #this will take a couple of hours; use a name to designate your run
# 	### for some reason, sometimes you have to exit the RGui and open again before the next step
#	  PaperPlots("runName") #this will produce the figures of the paper. 

