#Plotting Results from Cluster
require(ggplot2)
horizons <- c(1, 10, 50, 100, 500, 1000, 1000000)
# boxplots

# function producing boxplot data for one tree.
# For average of all trees, change files read and set t = 1
# lambda: the amount the knowledge delay is scaled
# save: number of species protected
# decision: decision frequency - number of years before the next decision is made
# final: year at which PD is recorded

# file <- 

require(ggplot)
require(reshape2)
################## data #####################
box_data <- function(t, final, decision, save, lambda) {
  box_dat <- matrix(ncol = 2, nrow = 700)
  k <- 0
  for (j in horizons){
    li <- paste("l100_", j, sep="")
    assign(li, list())
    for (i in k +1:100){
    
      fil <- paste("file", "/lambda_", lambda, "/saving_", save, "/decision_", decision, "/markov_", i, ".rda", sep ="")
      load(fil)
      #browser()
      box_dat[i,1] <- j
      box_dat[i,2] <- PhyDivs[final + 1]
    }
    k <- k + 100
  }
  box_dat <- data.frame(box_dat)
  names(box_dat) <- c("Time_Period", "value_100")
  box_dat$Time_Period <- as.factor(box_dat$Time_Period)
  
return(box_dat)
}


################### graphing #####################
graph_dat <- function(decision, save, lambda){
  for (f in c(100, 500)){
    box_dat <- box_data(1,f, decision, save, lambda)
    graph <- ggplot(box_dat, aes(x=Time_Period, y=value_100, fill = Time_Period), show.legend = F) + geom_boxplot(show.legend = F)
    graph <- graph + xlab("Extinction Time Horizon") + ylab ("Phylogenetic Diversity") + theme_minimal()
    #graph <- #graph + scale_fill_manual(values = c("No Conservation" = "grey"))
    #graph <- graph + scale_x_discrete(labels = c("Period 1" = "New Label", "Period 2", "Period 3"))
    #graph <- graph + ggtitle(paste("Decision Frequency: ", decision, ". Conservation capacity: "), save, sep = "")+ theme(plot.title = element_text(hjust = 0.5))
    
    assign(paste("tree", f, sep = "_"), graph)
    filename <- paste ("file", 500, "/save", save, "_decision_", decision, ".png",sep ="")
    #ggsave(filename, plot = graph, height = 3.808948, width =  4.38)
  }
  return(graph)
}

######################### Time Data #############################

time_data <- function(t, decision, save, lambda) {
  time_dat <- matrix(rep(0, 5002*7), ncol = 5002, nrow = 7)
  k <- 0
  for (j in 1:length(horizons)){
    li <- paste("l100_", j, sep="")
    assign(li, list())
    for (i in k +1:100){
      fil <-  paste("~/Markov_long/1lam/markov_", i + (t-1)*700, ".rda", sep ="")
      load(fil)
      time_dat[j,1:5001] <- PhyDivs + time_dat[j,1:5001]
      time_dat[j, 5002] <- horizons[j]
    }
    k <- k + 100
  }
  time_dat <- data.frame(time_dat)
  colnames(time_dat) <- c(1:5001, "Time_Period")
  time_dat[1:7,1:5001] <- time_dat[1:7,1:5001]/100
  time_dat$Time_Period <- as.factor(time_dat$Time_Period)
  
  return(time_dat)
} 

require(ggplot2)
require(reshape2)


#########################################






