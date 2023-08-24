#Plotting Results from Cluster
require(ggplot2)
horizons <- c(1, 10, 50, 100, 500, 1000, 1000000)
# boxplots

# function producing boxplot data for one tree.
# For average of all trees, change files read and set t = 1

box_data <- function(t, final, decision, save, lambda) {
  box_dat <- matrix(ncol = 2, nrow = 700)
  k <- 0
  for (j in horizons){
    li <- paste("l100_", j, sep="")
    assign(li, list())
    for (i in k +1:100){
      
      #fil <-  paste("~/Documents//markov_", i + (t-1)*700, ".rda", sep ="")
      #fil <-  paste("~/decision_10/time_alter_", i + (t-1)*700, ".rda", sep ="")
      fil <- paste("~/Documents/Research/Results/Model2/Markov/lambda_1/saving_0/decision_600/markov_", i, ".rda", sep ="")
      load(fil)
      #browser()
      box_dat[i,1] <- j
      box_dat[i,2] <- PhyDivs[final + 1]
      print(dim(PhyDivs))
    }
    #browser() 
    k <- k + 100
  }
  box_dat <- data.frame(box_dat)
  names(box_dat) <- c("Time_Period", "value_100")
  box_dat$Time_Period <- as.factor(box_dat$Time_Period)
  
return(box_dat)
}

box_list <- list()
for (t in 1:10){
  box_dat <- box_data(t, 100)
  assign(paste("box_500_", t ,sep = ""), box_dat)
  box_list <- append(box_list, list(box_dat))
}


require(ggplot2)
for (i in 1:10){
  assign("box_dat", as.data.frame(box_list[[i]]))
  graph <- ggplot(box_dat, aes(x=Time_Period, y=value_100, fill = Time_Period), show.legend = F) + geom_boxplot(show.legend = F)
  graph <- graph + xlab("Extinction Time Horizon") + ylab ("Phylogenetic Diversity")
  #graph <- graph + ggtitle(paste("Decision Frequency: ", decision, ". Conservation capacity: "), saved, sep = "")+ theme(plot.title = element_text(hjust = 0.5))
  
  assign(paste("tree_500", i, sep = "_"), graph)
  ggsave(paste("box_100_", i, ".png"), plot = graph)
}

####### For Averaged Trees #######
filename <- "saving_100/decision_100"

graph_dat <- function(decision, save, lambda){
  for (f in c(100, 500)){
    box_dat <- box_data(1,f, decision, save, lambda)
    graph <- ggplot(box_dat, aes(x=Time_Period, y=value_100, fill = Time_Period), show.legend = F) + geom_boxplot(show.legend = F)
    graph <- graph + xlab("Extinction Time Horizon") + ylab ("Phylogenetic Diversity") + theme_minimal()
    #graph <- #graph + scale_fill_manual(values = c("No Conservation" = "grey"))
    #graph <- graph + scale_x_discrete(labels = c("Period 1" = "New Label", "Period 2", "Period 3"))
    #graph <- graph + ggtitle(paste("Decision Frequency: ", decision, ". Conservation capacity: "), save, sep = "")+ theme(plot.title = element_text(hjust = 0.5))
    
    assign(paste("tree", f, sep = "_"), graph)
    filename <- paste ("~/Documents/Research/Results/Model1/graphs/f", 500, "/save", save, "_decision_", decision, ".png",sep ="")
    #ggsave(filename, plot = graph, height = 3.808948, width =  4.38)
  }
  return(graph)
}

for (s in c(100, 500)){
  for (d in c(1, 10, 100, 500))
  {
    graph_dat(d, s)
  }
}
for (d in c(1, 100, 500)){
  
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

mean_PD <- time_data(1,1, 100, 1 )

rownames(mean_PD) <- horizons
mean_PD <- melt(mean_PD)
names(mean_PD) <- c("pext_time", "timepoints", "PD")
graph <- ggplot(mean_PD, aes(x = timepoints, y = PD, colour = as.factor(pext_time))) + geom_line(aes(group = pext_time))
graph <- graph + labs(color='Time Horizons (years)'); graph
graph <- graph + scale_x_discrete(breaks = c(4000, 4200, 4400, 4600, 4800,5000))
graph <- graph + scale_x_discrete(breaks = c(1000, 2000, 3000, 4000,5000))


mean_PD$timepoints <- as.numeric(mean_PD$timepoints)
mean_PD <- mean_PD[mean_PD$timepoints < 500,]

save(graph, file = "500_years.rda")
graph <- graph + labs(colour = "Time Horizon (years)") + xlab("Time (years)") + ylab("Phylogenetic Diversity")
graph <- graph + theme_minimal()


graph <- graph + ggtitle(graphname)
nam <- paste("../Data/", filename, ".png", sep ="")
print(nam)
graph
ggsave(nam, plot = graph, width = 8, height = 4)




#########################################








graph_dat()
######## this is wrong ########
ggplot(combined_data, aes(x = Tree, y = N)) +
  geom_boxplot() +
  theme_minimal() +
  facet_wrap(~ Tree, scales = "free") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
