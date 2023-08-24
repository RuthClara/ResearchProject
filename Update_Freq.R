######### Update Frequency #########

source("../Markov_Model/Markov.R")
updates <- read.csv("../Data/Corrected_SpeciesHistory_June222022.csv")

max(updates$year) - min(updates$year) # = 27

lambda_dict <- {}


for (i in statuses){
  l = list()
  sp_in_state <- which(updates$category == i)
  for (j in sp_in_state){
    if (j == 1){next}
    if(updates[j-1, 'scientific_name'] == updates[j, 'scientific_name']){
      l <- append(l, updates[j-1,'year'] -  updates[j,'year'])
    }
  }
 lambda_dict[i] <- mean(unlist(l))
 
}

lambda_dict <- 1/lambda_dict

#average_no_yrs <- mean(table(updates$scientific_name[updates$scientific_name %in% sp_in_state]))
#lambda_dict[i] <- average_no_yrs/27 # rate of updates











######################### ROUGH ###############################
######### or would best way be to take number of LC species/unique(LC)
for ( i in statuses){
  state_updates <- updates[updates$category == i,]
  dim(state_updates)[1]/length(unique(state_updates$scientific_name))
}

###############################################################