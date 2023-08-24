####### Time Alter Condensed ######

require(phylobase)
require(data.table)
require(ape)

source("/rds/general/user/rcb22/home/Research_Proj/ExtTimeHorizon/Time_Alter_Functions.R")
source("/rds/general/user/rcb22/home/Research_Proj/ExtTimeHorizon/multtrees.R")

tree_list <- mam[[1]]
IUCNvals <- mam[[2]]

######### For Single tree #########
Time_Alter <- function(tree, P_50, pext_time, time_tot, number_protected, decision_interval, markov_pext, measure_interval = 1){
  
  names(P_50) <- c("Species", "Pext")
  tree_original <- tree
  time_points <- ceiling(time_tot/measure_interval)
  
  Protected <- matrix(nrow = time_points, ncol = number_protected)
  PDs <- matrix(nrow = 1, ncol = time_points + 1)
  colnames(PDs) <- as.character(measure_interval*(0:time_points))  
  
  # alter probabilities for EDGE scores
  Pext_EDGE <- change_pext_time(Pext, 50, pext_time)
  Pext_EDGE <- data.frame(P_50$Species, Pext_EDGE)
  
  # probabilities for Ext Sim
  Pext_Sim <- change_pext_time(markov_Pext$Pext, 50, measure_interval)
  Pext_Sim <- data.frame(markov_Pext$Species, Pext_Sim)
  
  extinctions <- list()
  
  PDs[1] <- PD(tree)

  for (r in 1:time_points){
    if ((r-1) %% decision_interval == 0){
    prot <- protect_decision(number_protected, Pext_EDGE, tree)
    }
    ext <- Ext_sim(tree_original, Pext_Sim)
    extinctions[[r]] <- ext
    
    tree <- drop_branches(tree, ext, prot)
    PDs[r+1] <- PD(tree)
    if (r %% 100 == 0){
    save(PDs, file = paste("/rds/general/user/rcb22/home/Research_Proj/results/single_trees/saving_", number_protected, "/decision_", decision_interval, "/time_alter_", iter, ".rda", sep = ""))     
    }
  }
  return(PDs)
}

Time_Alter_Mult <- function(tree_list, P_50, pext_time, time_tot, number_protected, decision_interval, markov_pext, measure_interval = 1){
  
  names(P_50) <- c("Species", "Pext")
  tree_original <- tree_list[[1]]
  time_points <- ceiling(time_tot/measure_interval)
  
  Protected <- matrix(nrow = time_points, ncol = number_protected)
  PDs <- matrix(nrow = 1, ncol = time_points + 1)
  colnames(PDs) <- as.character(measure_interval*(0:time_points))  
  
  # alter probabilities for EDGE scores
  Pext_EDGE <- change_pext_time(P_50$Pext, 50, pext_time)
  Pext_EDGE <- data.frame(P_50$Species, Pext_EDGE)
  
  # probabilities for Ext Sim

  Pext_Sim <- change_pext_time(markov_pext$Pext, 50, measure_interval)
  Pext_Sim <- data.frame(markov_pext$Species, Pext_Sim)
  
  extinctions <- list()
  
  PDs[1] <-mean_PD(tree_list)
  
  for (r in 1:time_points){
    if ((r-1) %% decision_interval == 0){
      prot <- mean_protect_decision(number_protected, Pext_EDGE, tree_list)
    }
    ext <- Ext_sim(tree_original, Pext_Sim)
    extinctions[[r]] <- ext
    
    tree_list <- drop_branches_all(tree_list, ext, prot)
    PDs[r+1] <- mean_PD(tree_list)
    if (r %% 100 == 0){
    save(PDs, file = paste("/rds/general/user/rcb22/home/Research_Proj/results/mark_prob/saving_", number_protected, "/decision_", decision_interval, "/time_alter_", iter, ".rda", sep = ""))
    }
  }
  return(PDs)
}
  
  
Time_results <- function(iter, tree_list, P_50, pext_time, time_tot, number_protected, decision_interval, markov_pext, measure_interval = 1){
  PhyDivs <- Time_Alter_Mult(tree_list, P_50, pext_time, time_tot, number_protected, decision_interval, markov_pext, measure_interval)
  save(PhyDivs, file = paste("/rds/general/user/rcb22/home/Research_Proj/results/long0/time_alter_", iter, ".rda", sep = ""))
  return ("Success")
}


Time_results_s <- function(iter, tree, P_50, pext_time, time_tot, number_protected, decision_interval, measure_interval = 1){
  PhyDivs <- Time_Alter(tree, P_50, pext_time, time_tot, number_protected, decision_interval, measure_interval)
  save(PhyDivs, file = paste("/rds/general/user/rcb22/home/Research_Proj/results/single_trees/saving_", numer_protected, "/decision_", decision_interval,"/time_alter_s", iter, ".rda", sep = ""))
  return ("Success")
}
