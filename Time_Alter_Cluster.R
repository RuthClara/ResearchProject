### time alter for cluster ###

### Altering time period probability of extinction is measured over ###
rm(list = ls())
graphics.off()

iter <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))

#source("/rds/general/user/rcb22/home/Research_Proj/ExtTimeHorizon/multtrees.R")
#source("/rds/general/user/rcb22/home/Research_Proj/input_files")
source("/rds/general/user/rcb22/home/Research_Proj/ExtTimeHorizon/Time_Alter.R")
load("/rds/general/user/rcb22/home/Research_Proj/input_files/Mammal_trees.rda")

set.seed(iter %% 100)

horizons <- c(1, 10, 50, 100, 500, 1000, 1000000)
s <- ceiling(iter/100)

Time_results(iter, tree_list, IUCNvals, horizons[s], 5000, 500, 5010, 10)
