############# Assigning median GE and defining Tree ##############

# idea: investigate better complementary sets of taxa - i.e saving the species with the highest EDGE score may not always be optimal

require(ape)

# load data
#load("/rds/general/user/rcb22/home/Research_Proj/input_files/100_mammal_trees_with_2020_RL.RData")
abbrev_dict <- c("Critically Endangered" = "CR", "Endangered" = "EN", "Vulnerable" = "VU",
"Near Threatened" = "NT", "Least Concern" = "LC", "Data Deficient" = "DD")

MedGE <- function(IUCNdat){
  IUCNvals <- data.frame(IUCNdat$Species, IUCNdat$RedListEval)
  names(IUCNvals) <- c("Species","RedListEval")
  
  # remove excess columns
  IUCNvals$RedListEval[is.na(IUCNvals$RedListEval)] <- "Data Deficient"
  IUCNvals <- IUCNvals[!(IUCNvals$RedListEval %in% c("Extinct in the Wild","EX")),]
  
  ## add median pext column ## (Data Deficient goes to mean)
  RedListDict <- c("Critically Endangered" = 0.97, "Endangered" = 0.97/2, "Vulnerable" = 0.97/4,
                   "Near Threatened" = 0.97/8, "Least Concern" = 0.97/16, "Data Deficient" = 0.165846,
                   "CR" = 0.97, "EN" = 0.97/2, "VU" = 0.97/4,
                   "NT" = 0.97/8, "LC" = 0.97/16, "DD" = 0.165846)
  
  IUCNvals$p50 <- unname(RedListDict[IUCNvals$RedListEval])
  rownames(IUCNvals) <- seq(1, nrow(IUCNvals))  
  return(IUCNvals)
}

  ############# Define Tree ##############
  
  # Mammal tree
##defTree <- function(i){
#  Tree <- phy.block.1000[[1]]
  # remove NAs
#  Tree <- keep.tip(Tree, IUCNvals$Species)
#  return(Tree)
#}

#IUCN <- MedGE(Species)
#rownames(IUCN) <- seq(1, nrow(IUCN))

#IUCNstat <- IUCN[c("Species", "RedListEval")]
#IUCNstat$RedListEval <- unname(abbrev_dict[IUCNstat$RedListEval])
#names(IUCNstat) <- c("Species", "Status")

#IUCNvals <- IUCN[c("Species", "p50")]
#names(IUCNvals) <- c("Species", "P_ext")

#Tree <- defTree()

#to make tree list:
#tree_list <- phy.block.1000[sample(1:100, 10)]
#for (tr in 1:length(tree_list)){
#  tree_list[[tr]] <- keep.tip(tree_list[[tr]], IUCNvals$Species)
#}


