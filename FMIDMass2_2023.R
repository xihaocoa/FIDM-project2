#install.packages("Matrix")
#install.packages("fields")
#install.packages("igraph")
#install.packages("matrixStats")

#This script contains all functions needed for FMIDM ASSESSMENT 2 (2023). 
#Once you have installed the relevant packages above, comment out lines 1-4. 
#DO NOT MODIFY THIS CODE, OR ADD AND NEW LINES! 
#Create a separate script to run codes and answer the questions. 
#Save it now, in the same place as the other practical materials. 
#Source this file by clicking "source" in the top-right corner of this window. 

library(Matrix)
library(fields)
library(igraph)
library(matrixStats)
library(ggplot2)


################ Network analysis ################
netMS1<- readRDS("netMS1.rds")
netMS2<- readRDS("netMS2.rds")
#Visualize an adjacency matrix:
visAdjMat <- function(A){
  A <- Matrix(A, sparse=FALSE)
  n <- nrow(A)
  cmap <- hcl.colors(11, "Greens3", rev = TRUE)
  cmap <- cmap[c(1,11)]
  image(A)
  #par(mar=c(0, 0, 0, 0))
  #image(A, useRaster=TRUE, axes=FALSE)
}

#Mean/variance of degree distribution:
calcDeg <- function(A){
  k <- rowSums(A)
  hist(k)
  k <- data.frame(mean(k), var(k))
  cnames <- c("mean", "var")
  colnames(k) <- cnames
  print(k)
}

#Mean degree of a neighbour:
calcDegNbr <- function(A){
  k <- rowSums(A)#All degrees of nodes
  k <- sum(k^2)/sum(k)
  print(k)
}

degree1<- calcDegNbr(netMS1)
degree2<- calcDegNbr(netMS2)



#Characteristic path length:
calcCPL <- function(A){
  #Note: this uses igraph - you won't need to look at the code
  cpl <- mean_distance(graph_from_adjacency_matrix(A, mode="undirected"))
  print(cpl)
}
calcCPL(netMS1)
calcCPL(netMS2)

#Local clustering coefficients:
calcCluster <- function(A){
  A[A>1] <- 1
  n <- nrow(A)
  X <- matrix(0, n)
  for (i in 1:n){
    vi <- A[i, ]
    vfind <- which(vi==1)
    lv <- length(vfind)
    Wi <- A[vfind, vfind]
    X[i, 1] <- sum(Wi)/lv/(lv-1)
  }
  #Summary stats:
  y <- data.frame(matrix(ncol = 9, nrow=1))#mean=numeric(), var=numeric(), p0=numeric(), p5=numeric(), p25=numeric(), p50=numeric(), p75=numeric(), p95=numeric(), p100=numeric())
  y[, 1] <- mean(X, na.rm=TRUE)
  y[, 2] <- var(X, na.rm=TRUE)
  y[, 3:9] <- quantile(t(X), probs=c(0, .05, .25, .5, .75, .95, 1), na.rm=TRUE)
  cnames <- c("mean", "var", "0%", "5%", "25%", "50%", "75%", "95%", "100%")
  colnames(y) <- cnames
  examClusterStats <- y
  print(y)
}

################ Epidemic simulations ################

#SIR model:
epiSIR <- function(A,R0,plot=FALSE){
  #R0 <- 2.2
  denom <- rowSums(A)
  denom <- sum(denom^2)/sum(denom)-1#Mean nbr degree -1
  gamma <- 1/5
  beta <- R0*gamma/denom
  n <- nrow(A)
  seed <- seq(1,2,1)#sample(n,seednum)
  deltat <- 1/12
  tstart <- 0
  tend <- 100
  S <- rep(1,n)
  I <- rep(0,n)
  I[seed] <- 1
  S <- S-I
  #FYI, the following is deliberately simple, therefore somewhat inefficient:
  tvec <- seq(tstart,tend,deltat)
  lt <- length(tvec)
  x <- rep(0,lt)
  x[1] <- length(seed)
  Sout <- matrix(0,n,lt)
  Sout[,1] <- S
  Iout <- matrix(0,n,lt)
  Iout[,1] <- I
  for (t in 2:lt){
    lambda <- beta*deltat*A%*%I
    probsS <- (1-exp(-lambda))*S
    Schoose <- S*runif(n)-probsS
    infect <- which(Schoose<0)
    S[infect] <- 0
    I[infect] <- 1
    probsI <- (1-exp(-gamma*deltat))*I
    Ichoose <- I*runif(n)-probsI
    recover <- which(Ichoose<0)
    I[recover] <- 0
    x[t] <- sum(I)/n
    Sout[,t] <- S
    Iout[,t] <- I
  }
  Xout <- 2*Sout+Iout
  
  if (plot==TRUE){
    plot(tvec,x)
    cmap <- hcl.colors(11, "YlOrRd", rev = TRUE)
    cmap <- cmap[c(1,6,11)]
    image(tvec,seq(1,n,1),t(Xout),
          xlab="Time", ylab="Node",
          col = cmap)
    legend(x="topright", legend=c("Sus","Inf","Rec"), fill=rev(cmap))
  }
  
  #print(paste("Final size: ",as.character(1-sum(S+I)/n)))
  return(Xout)
}
epiSIR(netMS1,2)
epiSIR(netMS2,2)
finalSizes <- function(x,n=200,h=5,nchild=2){
  #x is the output from the function "epiSIR"
  #nchild is the number of children per household
  #Output: proportion of children and proportion of adults infected
  xend<- x[,dim(x)[2]]
  totalchildren <- nchild*n/h
  totaladults <- n-totalchildren
  children_end <- xend[seq(1, totalchildren, 1)]
  adults_end <- xend[seq((totalchildren+1), n, 1)]
  out <- c(length(which(children_end==0))/totalchildren,length(which(adults_end==0))/totaladults)
  return(out)
}

propHHinfected <- function(x,n=200,h=5,nchild=2){
  #x is the output from the function "epiSIR"
  #nchild is the number of children per household
  #Output: number of households infected (irrespective of age)
  xend<- x[,dim(x)[2]]
  totalchildren <- nchild*n/h
  totaladults <- n-totalchildren
  children_end <- xend[seq(1, totalchildren, 1)]
  adults_end <- xend[seq((totalchildren+1), n, 1)]
  children_end <- matrix(children_end, nchild, n/h, byrow=F)
  adults_end <- matrix(adults_end, (h-nchild), n/h, byrow=F)
  hh_end <- rbind(children_end, adults_end)
  hh_end[which(hh_end>0)] <- 1
  hh_end <- 1-hh_end#swap 1's and 0's so >0 is infected
  hh_end <- colSums(hh_end)
  out <- length(which(hh_end>0))/n*h
  return(out)
}

################################################################

#Skeleton for loop for running multiple simulations
#Note: you will need to edit this if you want to use it for anything!
multiSim <- function(A,R0vec){
  numSims <- 10
  numParam <- length(R0vec)
  output <- matrix(0, numParam, 3)
  for (j in 1:numParam){
  R0j <- R0vec[j]
  hold <- matrix(0, numSims, 3)
  print(paste("Processing parameter value: ",as.character(R0j)))
    for (i in 1:numSims){
      xi <- epiSIR(A,R0j)#xi is the output of "epiSIR"
      hold[i,] <- c(propHHinfected(xi), finalSizes(xi))
      output
    }
  output[j, ] <- colMeans(hold)
  }
  return(output)
}

R0vec1<-seq(1,4,0.1)
result1<-multiSim(netMS1,R0vec1)
result2<-multiSim(netMS2,R0vec1)
plot(R0vec1, result1[, 1], type = "l", col = "black",ylim=c(0,1), xlab = "R0", 
     ylab = "Proportion of households infected", cex.axis = 1.5, cex.lab = 1.5,cex.main = 1.5,
     main = "Effect of R0 and different network topology on the proportion of households infected")
lines(R0vec1, result2[, 1], col = "red")
legend("topleft", legend = c("netMS1", "netMS2"), col = c("black", "red"),lty = 1,cex= 0.65)

plot(R0vec1, result1[, 2], type = "l", col = "black", ylim=c(0,1),xlab = "R0", 
     ylab = "Final size of epidemic", cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.5,
     main = "Effect of R0 and  different network topology on final size of epidemic among children and among adult")
lines(R0vec1, result1[, 3], col = "black", lty = 3)
lines(R0vec1, result2[, 2], col = "red")
lines(R0vec1, result2[, 3], col = "red", lty = 3)
legend("topleft", legend = c("Children (netMS1)", "Adults (netMS1)", "Children (netMS2)", "Adults (netMS2)"), 
       col = c("black", "black", "red", "red"), lty = c(1, 2, 1, 2), cex = 0.65)
a<-sd(result1[, 1])/mean(result1[, 1])
b<-sd(result2[, 1])/mean(result2[, 1])
################################################################

################ Make networks ################

#Ring lattice (used to make summer holiday networks):
netRing <- function(n,m){
  Aout <- matrix(0,n,n)
  Aout[upper.tri(Aout)] <- 1
  Am1 <- matrix(0,n-m,n-m)
  Am1[upper.tri(Am1)] <- -1
  Aout[seq(1,n-m,1),seq(m+1,n,1)] <- Aout[seq(1,n-m,1),seq(m+1,n,1)]+Am1
  Ap1 <- matrix(0,m+1,m+1)
  Ap1[upper.tri(Ap1)] <- 1
  Aout[seq(1,m+1,1),seq(n-m,n,1)] <- Aout[seq(1,m+1,1),seq(n-m,n,1)]+Ap1
  Aout <- Matrix(Aout,sparse=TRUE)
  Aout <- Aout+t(Aout)
}

#Summer holiday network:
netSummer <- function(case,m,n=200,h=5,nchild=2){
  #nchild is the number of children per household
  nadult <- h-nchild
  nextra <- 2*m#Number of extra households to add to children's contacts
  #
  #Re-define n if n/h not an integer:
  nhh <- floor(n/h)
  n <- h*nhh
  totalchildren=nhh*nchild
  totaladults=nhh*nadult
  A1 <- kronecker(diag(nhh), matrix(1,nchild,nchild))
  #
  if (case==1){
    diagpm <- netRing(nhh,m)
    A1b <- kronecker(diagpm, matrix(1,nchild,nchild))
    A1 <- A1+A1b
  } else if (case==2){
    newcontacts <- round(2*m*nchild*nhh)
    pool1 <- which(upper.tri(A1)==TRUE)#Indices of 0s in upper triangular part of A1
    pool2 <- which(A1==0)
    pool <- intersect(pool1,pool2)
    numposs <- length(pool)
    contacts <- sample(numposs, newcontacts ,replace = FALSE, prob = NULL)
    contacts <- pool[contacts]
    A1[contacts] <- 1
    A1 <- A1 + t(A1)
    A1[which(A1>1)] <- 1
  }
  #
  A2 <- kronecker(diag(nhh), matrix(1,nadult,nadult))
  A3 <- kronecker(diag(nhh), matrix(1,nchild,nadult))
  Ahh <- cbind(rbind(A1,t(A3)),rbind(A3,A2))
  Ahh <- Ahh- diag(diag(Ahh))
}
# Set parameters
n <- 200
h <- 5
nchild <- 2

# scenario 1
# Choose a range of m values
m_values <- seq(1, 15, 1)

# Initialize empty matrix to store results
results_case2 <- matrix(0, length(m_values), 3)

# Loop through m values
for (i in 1:length(m_values)) {
  m <- m_values[i]
  
  # Generate network for mixing scenario 1
  network <- netSummer(2, m, n, h, nchild)
  
  # Run simulations using multiSim function
  simulation_results_case2 <- multiSim(network, rep(1.5, length(m_values)))
  
  # Store the results
  results_case2[i, ] <- simulation_results_case2[i, ]
}

results_case_2 <- matrix(0, length(m_values), 3)

# Loop through m values
for (i in 1:length(m_values)) {
  m <- m_values[i]
  
  # Generate network for mixing scenario 1
  network <- netSummer(2, m, n, h, nchild)
  
  # R0 = 2.5
  simulation_results_case_2<- multiSim(network, rep(2.5, length(m_values)))
  
  # Store the results
  results_case_2[i, ] <- simulation_results_case_2[i, ]
}

# Plot the results
plot(m_values, results_case2[, 1], type = "l", ylim=c(0,1),xlab = "m", ylab = "Proportion of Households Infected",
     cex.axis = 1.5, cex.lab = 1.5,cex.main = 1.5,main = "Proportion of Households Infected (Mixing Scenario 2) with m(1-15),R0=1.5&2.5")
lines(m_values, results_case_2[, 1], col = "red")
legend("topleft", legend = c("R0= 1.5", "R0=2.5"), col = c("black", "red"),lty = 1,cex= 1.2)

# Add a line for the fixed R0 value
abline(h = R0, col = "green", lty = "dashed")


# scenario 2 
m_values <- seq(1, 15, 1)

# Initialize empty matrix to store results
results <- matrix(0, length(m_values), 3)

# Loop through m values
for (i in 1:length(m_values)) {
  m <- m_values[i]
  
  # Generate network for mixing scenario 1
  network <- netSummer(1, m, n, h, nchild)
  
  # Run simulations using multiSim function
  simulation_results <- multiSim(network, rep(1.5, length(m_values)))
  
  # Store the results
  results[i, ] <- simulation_results[i, ]
}

results_1 <- matrix(0, length(m_values), 3)

# Loop through m values
for (i in 1:length(m_values)) {
  m <- m_values[i]
  
  # Generate network for mixing scenario 1
  network <- netSummer(1, m, n, h, nchild)
  
  # R0 = 2.5
  simulation_results_1 <- multiSim(network, rep(2.5, length(m_values)))
  
  # Store the results
  results_1[i, ] <- simulation_results_1[i, ]
}

# Plot the results
plot(m_values, results[, 1], type = "l", ylim=c(0,1),xlab = "m", ylab = "Proportion of Households Infected",
     cex.axis = 1.5, cex.lab = 1.5,cex.main = 1.5,main = "Proportion of Households Infected (Mixing Scenario 1&2) with m(1-15),R0=1.5&2.5")
lines(m_values, results_1[, 1], col = "red")
lines(m_values,results_case2[,1],col= "black",lty=4)
lines(m_values,results_case_2[,1],col= "red",lty=4)
legend("bottomright", legend = c("R0= 1.5_S1", "R0=2.5_S1","R0= 1.5_S2","R0= 2.5_S2"), 
       col = c("black", "red","black","red"),lty = c(1, 1, 2, 2),,cex= 0.7,bty = "n",
       text.width = 0.5)
# Add a line for the fixed R0 value
abline(h = R0, col = "green", lty = "dashed")

# scenario 1
R0_values <- seq(1, 5, 0.1)
m_values <- 1:15
infected_households_1 <- matrix(0, nrow = length(R0_values), ncol = length(m_values))
for (i in seq_along(R0_values)) {
  for (j in seq_along(m_values)) {
    R0 <- R0_values[i]
    m <- m_values[j]
    net <- netSummer(case = 1, m = m)
    
    sim_result <- epiSIR(net, R0)
    infected_households_prop <- propHHinfected(sim_result)
    
    infected_households_1[i, j] <- infected_households_prop
  }
}
infected_households_df1 <- expand.grid(R0 = R0_values, m = m_values)
infected_households_df1$prop_infected_1 <- as.vector(infected_households_1)

ggplot(data = infected_households_df1, aes(x = R0, y = m, fill = prop_infected_1)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(title = "Proportion of Infected Households in scenario 1",
       x = "R0",
       y = "m",
       fill = "Proportion",cex.main = 1.5) +
  theme_minimal()
#new，for scenario 2
R0_values <- seq(1, 5, 0.1)
m_values <- 1:15
infected_households <- matrix(0, nrow = length(R0_values), ncol = length(m_values))
for (i in seq_along(R0_values)) {
  for (j in seq_along(m_values)) {
    R0 <- R0_values[i]
    m <- m_values[j]
    net <- netSummer(case = 2, m = m)
    
    sim_result <- epiSIR(net, R0)
    infected_households_prop <- propHHinfected(sim_result)
    
    infected_households[i, j] <- infected_households_prop
  }
}
infected_households_df <- expand.grid(R0 = R0_values, m = m_values)
infected_households_df$prop_infected <- as.vector(infected_households)

ggplot(data = infected_households_df, aes(x = R0, y = m, fill = prop_infected)) +
  geom_tile() +
  scale_fill_gradient(low = "white",  high = "red") +
  labs(title = "Proportion of Infected Households in scenario 2",
       x = "R0",
       y = "m",
       fill = "Proportion",cex.main = 1.5) +
  theme_minimal()






