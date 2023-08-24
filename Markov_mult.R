

########### function to update IUCN statuses using stochastic poisson dist ############

require(phylobase)

#source("../Markov_Model/Markov.R")
#source("../Markov_Model/Update_Freq.R")
#source("../Extinction_time_horizon/Time_Alter_Functions.R")
#source("../EDGE/Median_GE.R")
#source("../EDGE/EDGE2_Opt.R")

# RL_Status: list of Species and accompanying IUCN statuses 

#load("../Data/lambda_dict.rda")   # dictionary of poisson rates
#load("../Data/msm_model.rda")     
#load("../Data/transition_mat.rda")# M is the transition matrix describing the probabilities of the IUCN
 
#load("../Data/Mammal_trees.rda") # data.frame
 


source("/rds/general/user/rcb22/home/Research_Proj/Markov/Markov.R")
source("/rds/general/user/rcb22/home/Research_Proj/ExtTimeHorizon/Time_Alter_Functions.R")
source("/rds/general/user/rcb22/home/Research_Proj/input_files/Median_GE.R")
source("/rds/general/user/rcb22/home/Research_Proj/input_files/EDGE2_Opt.R")

load("/rds/general/user/rcb22/home/Research_Proj/data/lambda_dict.rda")   # dictionary of poisson rates
load("/rds/general/user/rcb22/home/Research_Proj/data/transition_mat.rda")# M is the transition matrix describing the probabilities of the IUCN
load("/rds/general/user/rcb22/home/Research_Proj/data/Mammal_trees.rda")  # data.frame


tree_list <- mam[[1]]
IUCNvals <- mam[[2]]
IUCNstat <- mam[[3]]
                                # status of a species changing from one to another within a year

statuses <- c("LC", "NT", "VU", "EN", "CR", "EX")

theoretical_IUCN <- IUCNstat      # Reported IUCN
real_life_IUCN <- IUCNstat        # IUCN statuses if we had all relevant knowledge

#theoretical_IUCN[sample(1:82,5),2] <- sample(statuses[1:5], 5, replace = TRUE) # this is wrong


###### provides random update to theoretical_IUCN over a year ########
status_update <- function(IUCN_t, IUCN_real, lambdas){
  for (i in statuses[-6]){
    
    number_in_status <- table(IUCN_t$Status)[i]
    lambda <- lambdas[i]*number_in_status
    number_changed <- rpois(1, lambda)
    if(is.na(number_changed)){
    number_changed <- 0
    }
    state_sps <- IUCN_t[IUCN_t$Status == i,]
    if (length(rownames(state_sps)) < number_changed){number_changed <- length(rownames(state_sps))}
    species_changed <- as.numeric(sample(rownames(state_sps), number_changed))
    if (max(species_changed) > 6251) {browser()}
    IUCN_t[species_changed,'Status'] <- IUCN_real[species_changed, 'Status']
    
      }
  return (IUCN_t)
}

########## Simulation ##########

### In one year ###

######## Year Long Sim #########
Knowledge_delay_year <- function(year, tree_list_t, tree_list_r, IUCN_t, IUCN_real, M, poisson_rates, number_protected, pext_time, decision_time){
  #browser()
#
  # Theoretical Changes
  IUCN_t <- status_update(IUCN_t, IUCN_real, poisson_rates)
  
  for (tre in 1:10){
    tree_list_t[[tre]] <- keep.tip(tree_list_t[[tre]], IUCN_t$Species[IUCN_t$Status != "EX"])
  }
  
  # Conservation Action
  pext <- IUCN_t
  names(pext) <- c("Species", "RedListEval")
  pext <- MedGE(pext)[c("Species", "p50")]
  p <- change_pext_time(pext$p50, 50, pext_time)
  pext <- data.frame(pext$Species, p)
  
  if ((year-1) %% decision_time == 0){
  protected <<- mean_protect_decision(number_protected, pext, tree_list_t) # does this change indices of species?
  }
  
  # Real World Changes
  IUCN_dist <- state_distr(IUCN_real)

  state_movements <- stoch_state_change_step(M, IUCN_dist)[[2]]
  IUCN_real <- allocate_changes(state_movements, IUCN_real, protected)

  for (tre in 1:10){
  tree_list_r[[tre]] <- keep.tip(tree_list_r[[tre]], IUCN_real$Species[IUCN_real$Status != "EX"])
  }
  PhyDiv <- mean_PD(tree_list_r)
  l <- list(PhyDiv, tree_list_t, tree_list_r, IUCN_t, IUCN_real)
  return(l)
}

########### How to recalculate extinction probability?? ############



##################################

# Function to repeat Knowledge_delay_year over a specified number of years
Knowledge_delay_sim <- function(tree_list, IUCN_t, IUCN_real, M, poisson_rates, number_protected, time_tot, pext_time, decision_time){
  start <- Sys.time()
  tree_list_t <- tree_list
  tree_list_r <- tree_list
  PhyDivs <- matrix(nrow = 1, ncol = time_tot + 1)
  PhyDivs[1] <- mean_PD(tree_list_r)
  for (t in 1:time_tot){
    #browser()
    results <- Knowledge_delay_year(t, tree_list_t, tree_list_r, IUCN_t, IUCN_real, M, poisson_rates, number_protected, pext_time, decision_time)
    tree_list_t <- results[[2]]
    tree_list_r <- results[[3]]
    IUCN_t <- results[[4]]
    IUCN_real <- results[[5]]
    PhyDivs[t+1] <- mean_PD(tree_list_r)
    if (difftime(Sys.time(),start) > 49*60*60) {save(PhyDivs, file = paste("n", filename, sep =""))}
  }
  return(PhyDivs)
}


# Function saves results of knowledge_delay_sim
Knowledge_results <- function(trees, IUCN_t, IUCN_real, M, poisson_rates, number_protected, time_tot, pext_time, decision_time, filename){
    
  PhyDivs <- Knowledge_delay_sim(trees, IUCN_t, IUCN_real, M, poisson_rates, number_protected, time_tot, pext_time, decision_time)
  print(PhyDivs)    
  save(PhyDivs, file = filename)
  return("Sucess")
}


  
###############################################################################################
