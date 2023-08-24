

########### function to update IUCN statuses using stochastic poisson dist ############


source("/rds/general/user/rcb22/home/Research_Proj/Markov/Markov.R")
#source("/rds/general/user/rcb22/home/Research_Proj/Markov/Update_Freq.R")
source("/rds/general/user/rcb22/home/Research_Proj/ExtTimeHorizon/Time_Alter_Functions.R")
source("/rds/general/user/rcb22/home/Research_Proj/input_files/Median_GE.R")
source("/rds/general/user/rcb22/home/Research_Proj/input_files/EDGE2_Opt.R")

require(phylobase)
require(ape)
# RL_Status: list of Species and accompanying IUCN statuses 

load("/rds/general/user/rcb22/home/Research_Proj/data/Mammal_trees.rda") # data.frame
tree_list <- mam[[1]]
IUCNvals <- mam[[2]]
IUCNstat <- mam[[3]]
load("/rds/general/user/rcb22/home/Research_Proj/data/lambda_dict.rda")   # dictionary of poisson rates
#load("../data/msm_model.rda")     
load("/rds/general/user/rcb22/home/Research_Proj/data/transition_mat.rda")# M is the transition matrix describing the probabilities of the IUCN
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
Knowledge_delay_year <- function(year, tree, IUCN_t, IUCN_real, M, poisson_rates, number_protected, pext_time, decision_time){
  
  tree <- keep.tip(tree, IUCN_real$Species[IUCN_real$Status != "EX"])
 
  # Theoretical Changes
  IUCN_t <- status_update(IUCN_t, IUCN_real, poisson_rates)
  
  # Conservation Action
  pext <- IUCN_t
  names(pext) <- c("Species", "RedListEval")
  pext <- MedGE(pext)[c("Species", "p50")]
  p <- change_pext_time(pext$p50, 50, pext_time)
  pext <- data.frame(pext$Species, p)
  
  if ((year-1) %% decision_time == 0){
  protected <<- protect_decision(number_protected, pext, tree) # does this change indices of species?
  }
  
  # Real World Changes
  IUCN_dist <- state_distr(IUCN_real)

  state_movements <- stoch_state_change_step(M, IUCN_dist)[[2]]
  IUCN_real <- allocate_changes(state_movements, IUCN_real, protected)

  tree <- keep.tip(tree, IUCN_real$Species[IUCN_real$Status != "EX"])
  PhyDiv <- PD(tree)
  l <- list(PhyDiv, tree, IUCN_t, IUCN_real)
  return(l)
}

########### How to recalculate extinction probability?? ############



##################################

# Function to repeat Knowledge_delay_year over a specified number of years
Knowledge_delay_sim <- function(tree, IUCN_t, IUCN_real, M, poisson_rates, number_protected, time_tot, pext_time, decision_time){

  PhyDivs <- matrix(nrow = 1, ncol = time_tot + 1)
  PhyDivs[1] <- PD(tree)
  for (t in 1:time_tot){
    #browser()
    results <- Knowledge_delay_year(t, tree, IUCN_t, IUCN_real, M, poisson_rates, number_protected, pext_time, decision_time)
    tree <- results[[2]]
    IUCN_t <- results[[3]]
    IUCN_real <- results[[4]]
    PhyDivs[t+1] <- PD(tree)
  }
  return(PhyDivs)
}


# Function saves results of knowledge_delay_sim
Knowledge_results <- function(tree, IUCN_t, IUCN_real, M, poisson_rates, number_protected, time_tot, pext_time, decision_time, filename){
    
  PhyDivs <- Knowledge_delay_sim(tree, IUCN_t, IUCN_real, M, poisson_rates, number_protected, time_tot, pext_time, decision_time)
    save(PhyDivs, file = filename)
    return("Sucess")
}


  
###############################################################################################
  

