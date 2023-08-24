#### Knowledge - Action Delay ####
source("/rds/general/user/rcb22/home/Research_Proj/ExtTimeHorizon/Time_Alter_Functions.R")

# transition matrix

RedListDict <- c("Critically Endangered" = 0.97, "Endangered" = 0.97/2, "Vulnerable" = 0.97/4,
                 "Near Threatened" = 0.97/8, "Least Concern" = 0.97/16)

### Example tree and Red List Statuses ###
#tre <- read.tree("../Data/anole.tre")

# status order
statuses <- c("LC", "NT", "VU", "EN", "CR", "EX")
#RL_Status <- data.frame(tre$tip.label, sample(statuses,length(tre$tip.label), replace = TRUE))
#names(RL_Status) <- c("Species", "Status")

## get the 1 year extinction probabilities for each status ##

LC <- change_pext_time(0.97/16, 50, 1)
NT <- change_pext_time(0.97/8, 50, 1)
# vulnerable: >10% in 100 years
VN <- change_pext_time(0.1, 100, 1)
# endangered: > 20% in 20 years
EN <- change_pext_time(0.2, 20, 1)
# critically endangered: >50% in 10 years
CR <- change_pext_time(0.50, 10, 1)

## function to calculate distribution of red list statuses
# table is dataframe displaying Red List Status of Each Species

state_distr <- function(status_table){
  freqs <- table(status_table$Status)
  freq_vec <- rep(0,6)
  for (i in 1:6){
    freq_vec[i] <- freqs[statuses[i]]
  }
  freq_vec[is.na(freq_vec)] <- 0
  return(freq_vec)
}

# expected change of distributions between states
Exp_state_change <- function(t_matrix, vec, time = 1){
  for (i in 1:time){
   vec <- t_matrix%*%vec 
  }
  return(as.vector(vec))
}

# carries out stochastic state change over a year, returning matrix of movements from states
stoch_state_change_step <- function(t_matrix, vec){
  movement <- matrix(0, 6, 6)
  v <- rep(0, 6)
  for (j in 1:6){
    tot <- vec[j]
    for (i in 1:5){
      change <- rbinom(1, tot, M[i,j]/(1 - sum(M[-(i:6),j])))
      v[i] <- v[i] + change
      tot <- tot - change
      movement[i,j] <- change
    }
    v[6] <- v[6]+tot
    movement[6,j] <- tot
  }

  return(list(v,movement))
}

# state change for specified number of time steps (years)
# need to modify for movement matrix
stoch_state_change <- function(t_matrix, vec, time = 1){
  for (i in 1:time){
    vec <- stoch_state_change_step(t_matrix, vec)
  }
  return(vec)
}

# allocates changes in Markov chain to specific species.
allocate_changes <- function(change_M, real_IUCN, protected){
  for (i in 1:5){       # from state 
    state <- statuses[i]
    for (j in 1:6){     # to state
      if (i != j){
        state_IUCN <- real_IUCN[real_IUCN$Status == state,]
        no_state <- length(state_IUCN$Species)
        #if (is.na(state_IUCN)){
        #  no_state <- 0
        #}
        #browser()
        inds <- sample(length(state_IUCN$Species), change_M[j,i]) # to change
        sp <- real_IUCN[real_IUCN$Status == state,1][inds]
        if (any(sp %in% protected) && j>i){
          #browser()
          # remove index
          #browser()
          for_rm <- which(sp %in% protected)
          inds <- subset(inds, !(inds %in% inds[for_rm]))
        }
        #browser()
          real_IUCN[real_IUCN$Status == state,2][inds] <- statuses[j]
      #browser()
      }
    }
  }
  return(real_IUCN)
}

################################################################################
