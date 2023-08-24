########## Time Alter Functions ##########

require(phylobase)
require(data.table)
require(ape)

#source("/rds/general/user/rcb22/home/Research_Proj/input_files/EDGE2_Opt.R")

### Change probability time period from x to y ###

change_pext_time <- function(x_pext, x, y, mat = F){
  y_pext <- 1 - (1 - x_pext)^(y/x)
  #print(typ)
  if (mat == T){   # this isn't correct #
    for (i in 1:ncol(y_pext)){
      y_pext[i,i] <- y_pext[i,i] + 1 - sum(y_pext[,i])
    } 
  }
  
  return(y_pext)
}

### Decide protection ###
# reevaluates after every species chosen to protect
protect_dependent <- function(number_protected, P_ext, tree){
  protected <- matrix(nrow = 1, ncol = number_protected)
  for (i in 1:number_protected){

    EDGE_dat <- EDGE2_mod(tree, P_ext)[1][[1]]
    protected[i] <- EDGE_dat[EDGE_dat$EDGE ==max(EDGE_dat$EDGE),'Species']
    P_ext[P_ext$Species == protected[[i]], 'P_ext'] <- 0
  }
  return(protected)
}

# simple top N species
protect_decision <- function(number_protected, P_ext, tree){
     
      EDGE_dat <- EDGE2_mod(tree, P_ext)[1][[1]]
      protected <- EDGE_dat[EDGE_dat$EDGE %in% tail(sort(EDGE_dat$EDGE), number_protected),]$Species[1:number_protected]

      return(protected)
}
###### Simulate Extinctions ######
# over given time period with appropriate P_ext given for time period.
Ext_sim <- function(tree, P_ext){   # P_ext 2 columns, Species and Pext
  names(P_ext) <- c("Species", "Pext")
  extinctions <- list()
    # simulate extinctions for all trials
  for (sp in tree$tip.label){
    if (runif(1) < P_ext[which(P_ext$Species == sp), 'Pext']){
        extinctions <- append(extinctions, sp)
    }
  }
  return (unlist(extinctions))
}    

########### Drop Branches if not protected ##############

drop_branches <- function(tree, extinctions, protected){
  for (sp in extinctions){
    if (sp %in% protected == F){
      tree <- drop.tip(tree, sp)
    }
  }
  return(tree)
}
      
# returns phylogenetic diversity (PD) of tree
PD <- function(tree){
  if (class(tree) == "phylo"){
  return(sum(tree$edge.length))}
  else {
    return(sum(tree@edge.length[!is.na(tree@edge.length)]))
  }
}

# returns mean PD of multiple trees
mean_PD <- function(tree_list){
  print(tree_list)
  tot <- 0
  for (tr in tree_list){
    tot <- PD(tr) + tot
  }
  return(tot/length(tree_list)) 
}

# drops branches of all trees in species not protected
drop_branches_all <- function(tree_list, extinctions, protected){
  for (sp in extinctions){
    if (sp %in% protected == F){
      for (tr in 1:length(tree_list)){
        tree_list[[tr]] <- drop.tip(tree_list[[tr]], sp)
      }
    }
  }
  return(tree_list)
}

## protection decision from mean EDGE
#mean_protect_decision <- function(number_protected, P_ext, tree_list){
#
#      EDGE_dat <- mean_EDGE(tree_list, P_ext)[1][[1]]
#      protected <- EDGE_dat[EDGE_dat$EDGE %in% tail(sort(EDGE_dat$EDGE), number_protected),]$Species[1:number_protected]
#
#      return(protected)
#}

## protection decision from mean EDGE, reevaluating after each species chosen to protect
mean_protect_decision <- function(number_protected, P_ext, tree_list){
  
  EDGE_dat <- mean_EDGE(tree_list, P_ext)
  names(EDGE_dat)[1] <- "Species"
  protected <- EDGE_dat[EDGE_dat$EDGE %in% tail(sort(EDGE_dat$EDGE), number_protected),]$Species[1:number_protected]
  
  return(protected)
}

mean_protect_dependent <- function(number_protected, P_ext, tree){
  protected <- matrix(nrow = 1, ncol = number_protected)
  for (i in 1:number_protected){
    
    EDGE_dat <- mean_EDGE(tree_list, P_ext)
    protected[i] <- EDGE_dat[EDGE_dat$EDGE ==max(EDGE_dat$EDGE),'Species']
    P_ext[P_ext$Species == protected[[i]], 'P_ext'] <- 0
  }
  return(protected)
}

## Returns mean EDGE score of multiple trees
mean_EDGE <- function(tree_list, pext){
  names(pext) <- c("Species", "Pext")
  #rownum <- length(keep.tip(tree_list[[1]], pext$Species)$tip.label)
  rownum <- length(tree_list[[1]]$tip.label)
  EDGE_score <- data.frame(matrix(rep(0, 4*rownum),ncol = 4, nrow=rownum))
  for (tr in tree_list){
    #Tree <- keep.tip(tr, IUCNvals$Species)
    EDGE_score <- EDGE2_mod(tr, pext)[[1]][,-1] + EDGE_score
  }
  EDGE_score <- EDGE_score/length(tree_list)
  return(cbind(EDGE2_mod(tr, pext)[[1]][,1], EDGE_score))
}

## Add column determining IUCN state ##


RedListDict <- c("Critically Endangered" = 0.97, "Endangered" = 0.97/2, "Vulnerable" = 0.97/4,
                 "Near Threatened" = 0.97/8, "Least Concern" = 0.97/16)


#IUCN_status <- function(dat){
#  pext <- dat$Pext
#  IUCN <- matrix(nrow = 1, ncol = length(pext))
  
#}


########## IUCN_t still crecords extinctions which are not on IUCN_real ##
############ ???????????~~ ~~~~~~~~~####################################
