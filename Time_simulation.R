library(dplyr)
library(readr)
library(softImpute)
library(tidyverse)
library(tictoc)
library(RcppArmadillo)
library(Rcpp)
library(ggplot2)
library(openxlsx)

sourceCpp("~/Dropbox/Uni/Master_Econometrie/Thesis/Code/R/gammaui.cpp")
source("Thesis_methods.R")

#Timing experiments-------

#' Get data function to simulate a dataset
#'
#' @param n various input parameters for the model
#' @param i 
#' @param dx 
#' @param dg 
#' @param factors 
#' @param sparsity 
#' @param unknown 
#' @param lambda 
#' @param f 
#' @param epsilon 
#'
#' @return List with datasets
#' @export
#'
#' @examples
Get_data <- function(n, i, dx, dg, factors, sparsity, unknown, lambda,f, epsilon){
  alpha <- runif(n, min=-5, max=0)
  beta <- runif(i, min=-5, max=0)
  X <- matrix(rnorm(n * dx, 0, 1), n, dx)
  G <-  matrix(rnorm(i * dg, 0, 1), i, dg)
  C <- matrix(rnorm(dx * factors, 0, 0.1), dx, factors)
  D <- matrix(rnorm(dg * factors, 0, 0.1), dg, factors)
  low_rankC <- cbind(X%*%C, alpha, rep(1, n))
  low_rankD <- cbind(G%*%D, rep(1,i), beta)
  gamma <- low_rankC %*% t(low_rankD) + matrix(rnorm(n*i, 0, 1), n, i)
  probability <- exp(gamma) / (1 + exp(gamma))
  
  # Create a train subset with a certain sparsity level
  if(sparsity != 1){
    USERID <- sample(1:n, sparsity*n*i, replace = TRUE)
    OFFERID <- sample(1:i, sparsity*n*i, replace = TRUE)
  }
  else{
    USERID <- rep((1:n), i)
    OFFERID <- rep((1:i), each = n)
  }
  df <- data.frame("USERID_ind_new" = USERID, "OFFERID_ind_new" = OFFERID)
  df <- unique(df)
  df$CLICK <- probability[as.matrix(df[ , c("USERID_ind_new", "OFFERID_ind_new")])]
  df$CLICK <- as.numeric(df$CLICK > 0.5)
  
  #Now we gather a subset of the rows in X and G
  #First we add indices to X and G
  X <- as.data.frame(X)
  G <- as.data.frame(G)
  X["USERID_ind_new"] <- 1:n
  G["OFFERID_ind_new"] <- 1:i
  X <- X %>% select("USERID_ind_new", everything())
  G <- G %>% select("OFFERID_ind_new", everything())
  
  #Draw a subset of random indices
  X <- X[1:(unknown*n),]
  G <- G[1:(unknown*i),]
  
  #Prepare the datasets
  prelim <- formatdfXG(df, X, G)
  df <- prelim$df
  X <- prelim$X
  G <- prelim$G
  
  return(list(df = df, X=X, G=G))
}

#' Simulates and returns full sets of data
#'
#' @param n 
#' @param i 
#' @param dx 
#' @param dg 
#' @param factors 
#' @param sparsity 
#' @param unknown 
#' @param lambda 
#' @param f 
#' @param epsilon 
#'
#' @return
#' @export
#'
#' @examples
Simulation <- function(n, i, dx, dg, factors, sparsity, unknown, lambda,f, epsilon){
  data <- Get_data(n, i, dx, dg, factors, sparsity, unknown, lambda,f, epsilon)
  df <- data$df
  X <- data$X
  G <- data$G
  
  #Estimate all methods on the dataset
  iter <- 1000
  iter_CD <- 50
  cat("--------------", "Running old model", "--------------", '\n')
  results_old <- parEst(df, factors, lambda, iter, 2, llh = T, T, epsilon = epsilon)
  cat("--------------", "Running SVD model", "--------------", '\n')
  results_SVD <- ParEst_rest_parrest_SVD(df, X, G, factors, lambda, iter, iter_CD, epsilon = epsilon, T, T )
  cat("--------------", "Running Eigen model", "--------------", '\n')
  results_Eigen <- ParEst_rest_parrest_Eigen(df, X, G, factors, lambda, iter, iter_CD, epsilon = epsilon, T, T )
  cat("--------------", "Running Deriv model", "--------------", '\n')
  results_deriv <- ParEst_rest_parrest_deriv(df, X, G, factors, lambda, iter, iter_CD, epsilon = epsilon, T, T)
  
  parameters <- list( "n"= n, "i"= i, "dx" = dx, "dg"=dg, "f"=f, "sparsity"=sparsity,
                      "unknown"=unknown, "factors" = factors, "lambda"=lambda, "iter" = iter, "iter_CD"=iter_CD, epsilon = epsilon)
  return(list("results_old" = results_old, "results_SVD" = results_SVD, 
              "results_Eigen" = results_Eigen, "results_deriv" = results_deriv, "parameters" = parameters))
}


#' Returns parameter estimates for full extra info (unknown = 1)
#'
#' @param n 
#' @param i 
#' @param dx 
#' @param dg 
#' @param factors 
#' @param sparsity 
#' @param unknown 
#' @param lambda 
#' @param f 
#' @param epsilon 
#'
#' @return
#' @export
#'
#' @examples
Simulation_full_extra_inf <- function(n, i, dx, dg, factors, sparsity, unknown, lambda,f, epsilon){
  alpha <- runif(n, min=-5, max=0)
  beta <- runif(i, min=-5, max=0)
  X <- matrix(rnorm(n * dx, 0, 1), n, dx)
  G <-  matrix(rnorm(i * dg, 0, 1), i, dg)
  C <- matrix(rnorm(dx * factors, 0, 0.1), dx, factors)
  D <- matrix(rnorm(dg * factors, 0, 0.1), dg, factors)
  varepsilon <- matrix(rnorm(n*i, 0, 1), n, i)
  gamma <- alpha%*%t(rep(1,i)) + rep(1,n)%*%t(beta) + X%*%C %*% t(D)%*%t(G) + varepsilon
  probability <- exp(gamma) / (1 + exp(gamma))
  
  # Create a train subset with a certain sparsity level
  if(sparsity != 1){
    USERID <- sample(1:n, sparsity*n*i, replace = TRUE)
    OFFERID <- sample(1:i, sparsity*n*i, replace = TRUE)
  }
  else{
    USERID <- rep((1:n), i)
    OFFERID <- rep((1:i), each = n)
  }
  df <- data.frame("USERID_ind" = USERID, "OFFERID_ind" = OFFERID)
  df <- unique(df)
  df$CLICK <- probability[as.matrix(df[ , c("USERID_ind", "OFFERID_ind")])]
  df$CLICK <- as.numeric(df$CLICK > 0.5)
  
  #Now we gather a subset of the rows in X and G
  #First we add indices to X and G
  X <- as.data.frame(X)
  G <- as.data.frame(G)
  X["USERID_ind"] <- 1:n
  G["OFFERID_ind"] <- 1:i
  X <- X %>% select("USERID_ind", everything())
  G <- G %>% select("OFFERID_ind", everything())
  
  X <- X[1:(unknown*n),]
  G <- G[1:(unknown*i),]
  
  #Prepare the datasets
  prelim <- formatdfXG(df, X, G)
  df <- prelim$df
  X <- prelim$X
  G <- prelim$G
  
  #Estimate all methods on the dataset
  factors <- f
  lambda <- lambda
  iter <- 10000
  iter_CD <- 10000
  cat("--------------", "Running old model", "--------------", '\n')
  results_old <- parEst(df, factors, lambda, iter, 2, llh = T, T, epsilon = epsilon)
  cat("--------------", "Running SVD model", "--------------", '\n')
  results_full <- ParEst_rest_fullrest(df, X, G, factors, lambda, iter, iter_CD, epsilon = epsilon, T, T )
  
  parameters <- list( "n"= n, "i"= i, "dx" = dx, "dg"=dg, "f"=f, "sparsity"=sparsity,
                      "unknown"=unknown, "factors" = factors, "lambda"=lambda, "iter" = iter, "iter_CD"=iter_CD, epsilon = epsilon)
  return(list("results_old" = results_old, "results_full" = results_full, "parameters" = parameters))
}

#' Timing simulation
#'
#' @param n 
#' @param i 
#' @param dx 
#' @param dg 
#' @param factors 
#' @param sparsity 
#' @param unknown 
#' @param lambda 
#' @param f 
#' @param epsilon 
#'
#' @return Data frame with different convergence times and set of parameters
#' @export
#'
#' @examples
Simulation_time <- function(n, i, dx, dg, factors, sparsity, unknown, lambda,f, epsilon){
  
  #Define additional variables
  iter <- 1000
  iter_CD <- 50
  
  output <- data.frame(n=n, Old= NA, SVD =NA, Eigen = NA, Deriv = rep(NA,5))
  for(run in 1:5){
    data <- Get_data(n, i, dx, dg, factors, sparsity, unknown, lambda,f, epsilon)
    df <- data$df
    X <- data$X
    G <- data$G
    
    #Run the models
    cat("--------------", "Running OLD model", "--------------", '\n')
    results_old <- parEst(df, factors, lambda, iter, 2, llh = T, T, epsilon = epsilon)
    output[run,2] <- sum(results_old$time)
    cat("--------------", "Running SVD model", "--------------", '\n')
    results_SVD <- ParEst_rest_parrest_SVD(df, X, G, factors, lambda, iter, iter_CD, epsilon = epsilon, T, T )
    output[run,3] <- sum(results_SVD$time)
    cat("--------------", "Running Eigen model", "--------------", '\n')
    results_Eigen <- ParEst_rest_parrest_Eigen(df, X, G, factors, lambda, iter, iter_CD, epsilon = epsilon, T, T )
    output[run,4] <- sum(results_Eigen$time)
    cat("--------------", "Running Deriv model", "--------------", '\n')
    results_deriv <- ParEst_rest_parrest_deriv(df, X, G, factors, lambda, iter, iter_CD, epsilon = epsilon, T, T)
    output[run,5] <- sum(results_deriv$time)
  }
  
  parameters <- list( "n"= n, "i"= i, "dx" = dx, "dg"=dg, "f"=f, "sparsity"=sparsity,
                      "unknown"=unknown, "factors" = factors, "lambda"=lambda, "iter" = iter, "iter_CD"=iter_CD, epsilon = epsilon)
  return(list("output"= output, "parameters" = parameters))
}

#' Same as Simulation_time but now for full unknown unknown = 1
#'
#' @param n 
#' @param i 
#' @param dx 
#' @param dg 
#' @param factors 
#' @param sparsity 
#' @param unknown 
#' @param lambda 
#' @param f 
#' @param epsilon 
#'
#' @return
#' @export
#'
#' @examples
Simulation_time_full <- function(n, i, dx, dg, factors, sparsity, unknown, lambda,f, epsilon){
  
  #Define additional variables
  iter <- 1000
  iter_CD <- 50
  
  output <- data.frame(n=n, Old= NA, full = rep(NA,5))
  for(run in 1:5){
    data <- Get_data(n, i, dx, dg, factors, sparsity, unknown, lambda, f, epsilon)
    df <- data$df
    X <- data$X
    G <- data$G
    
    #Run the models
    cat("--------------", "Running old model", "--------------", '\n')
    results_old <- parEst(df, factors, lambda, iter, 2, llh = T, T, epsilon = epsilon)
    output[run,2] <- sum(results_old$time)
    cat("--------------", "Running SVD model", "--------------", '\n')
    results_full <- ParEst_rest_fullrest(df, X, G, factors, lambda, iter, iter_CD, epsilon = epsilon, T, T )
    output[run,3] <- sum(results_full$time)
  }
  
  parameters <- list( "n"= n, "i"= i, "dx" = dx, "dg"=dg, "f"=f, "sparsity"=sparsity,
                      "unknown"=unknown, "factors" = factors, "lambda"=lambda, "iter" = iter, "iter_CD"=iter_CD, epsilon = epsilon)
  return(list("output"= output, "parameters" = parameters))
}


#' This function plots the results from the Simulation function
#'
#' @param results 
#'
#' @return
#' @export
#'
#' @examples
plot_results_full <- function(results){
  #Define variables used in the simulation
  n <- results$parameters$n
  i <- results$parameters$i
  dx <- results$parameters$dx
  dg <- results$parameters$dg
  f <- results$parameters$f
  sparsity <- results$parameters$sparsity
  unknown <- results$parameters$unknown
  lambda <- results$parameters$lambda
  
  #Retrieve all the results
  results_old <- results$results_old
  results_full <- results$results_full
  
  #Quick function which transforms data to wanted format
  getdf <- function(results_method, name){
    time <- results_method$time
    cumtime <- as.data.frame(cumsum(time))
    cumtime <- rbind(0, cumtime)
    logllh <- results_method$logllh
    relativellh <- logllh/logllh[1]
    df <- data.frame("logllh" = relativellh, "time" = as.numeric(unlist(cumtime)), model = name)
    return(df)
  }
  
  #transform all data to wanted format
  df_old <- getdf(results_old, "Unrestricted")
  df_full <- getdf(results_full, "Fully restricted")
  
  #Bind columns together
  df <- rbind(df_old, df_full)
  
  #Plot the graph
  chart_title <- paste(paste("n =",n),paste("i =",i),paste("delta","=",sparsity), paste("factors =", factors), 
                       paste("f =",f),paste("zeta","=",unknown),paste("lambda","=",lambda), sep = ", ")
  ggplot(data = df, aes(x = time, y = logllh, color = model)) + geom_line(size = 1, alpha = 0.5) +
    theme_bw()+ xlab("Time in seconds") + ylab("Relative objective function")+ 
    theme(legend.position = c(0.81, 0.8))
}

#' This function plots the results from the simulation_full_extra_info function
#'
#' @param results 
#'
#' @return
#' @export
#'
#' @examples
plot_results <- function(results){
  #Define variables used in the simulation
  n <- results$parameters$n
  i <- results$parameters$i
  dx <- results$parameters$dx
  dg <- results$parameters$dg
  f <- results$parameters$f
  sparsity <- results$parameters$sparsity
  unknown <- results$parameters$unknown
  lambda <- results$parameters$lambda
  
  #Retrieve all the results
  results_old <- results$results_old
  results_SVD <- results$results_SVD
  results_Eigen <- results$results_Eigen
  results_deriv <- results$results_deriv
  
  #Quick function which transforms data to wanted format
  getdf <- function(results_method, name){
    time <- results_method$time
    cumtime <- as.data.frame(cumsum(time))
    cumtime <- rbind(0, cumtime)
    logllh <- results_method$logllh
    relativellh <- logllh/logllh[1]
    df <- data.frame("logllh" = relativellh, "time" = as.numeric(unlist(cumtime)), model = name)
    return(df)
  }
  
  #transform all data to wanted format
  df_old <- getdf(results_old, "Unrestricted")
  df_SVD <- getdf(results_SVD, "SVD")
  df_Eigen <- getdf(results_Eigen, "Eigen")
  df_Deriv <- getdf(results_deriv, "Derivative")
  
  #Bind columns together
  df <- rbind(df_old, df_SVD, df_Eigen, df_Deriv)
  
  #Plot the graph
  chart_title <- paste(paste("n =",n),paste("i =",i),paste("delta","=",sparsity), paste("factors =", factors), 
                       paste("f =",f),paste("zeta","=",unknown),paste("lambda","=",lambda), sep = ", ")
  ggplot(data = df, aes(x = time, y = logllh, color = model)) + geom_line(size = 1, alpha = 0.5) +
    theme_bw()+ xlab("Time in seconds") + ylab("Relative objective function")+ 
    theme(legend.position = c(0.83, 0.67))
}

#Visualizing the probability distribution
#df_prop <- data.frame("Probability" = as.vector(probability))
#Probabilities <- ggplot(data = df_prop, aes(Probability)) + geom_histogram(bins = 15, color = "black", fill = "grey") + 
#  geom_vline(aes(xintercept=mean(Probability)), linetype = "dashed") + scale_y_log10(breaks = c(0,10,100,1000,10000,100000)) + xlab("Probability") + ylab("Frequency")+
#  theme_bw()
#Probabilities

#Timing simulation setup----------

#Different simulations for nonsparse
i <- 100
dx <- 100
dg <- 100
sparsity <- 0.05
unknown <- 0.5
lambda <- 0.5
factors <- 5
f <- 5
epsilon <- 0.01

#Vector with options for N
N <- c(10^2, round(10^(2.5)), 10^3, round(10^(3.5)), 10^4, round(10^(4.5)), 10^5 )

#Define iteration as zero
iteration <- 0

#formulate results list
results_sparse <- list()

for(n in N){
  name <- paste("results_", as.character(iteration), sep = "")
  results_sparse[[name]] <- Simulation_time(n=n,i=i,dx=dx,dg=dg,factors=factors,sparsity=sparsity,unknown=unknown,lambda=lambda
                          , f=f, epsilon = epsilon)
  iteration = iteration + 1
}

saveRDS(results_sparse, "Results_sparse_5.RDS")

#Second simulation, now set sparsity to 1
sparsity <- 1

#Define iteration as zero
iteration <- 0
results_dense <- list()

#formulate results list
for(n in N){
  name <- paste("results_", as.character(iteration), sep = "")
  results_dense[[name]] <- Simulation_time(n=n,i=i,dx=dx,dg=dg,factors=factors,sparsity=sparsity,unknown=unknown,lambda=lambda
                                       , f=f, epsilon = epsilon)
  iteration = iteration + 1
}

saveRDS(results_dense, "Results_dense_5.RDS")

#Same data for full additional data
sparsity <- 0.05
unknown <- 1
epsilon <- 0.005

#Define iteration as zero
iteration <- 0

#formulate results list
results_sparse_full_small <- list()

for(n in N){
  name <- paste("results_", as.character(iteration), sep = "")
  results_sparse_full_small[[name]] <- Simulation_time_full(n=n,i=i,dx=dx,dg=dg,factors=factors,sparsity=sparsity,unknown=unknown,lambda=lambda
                                            , f=f, epsilon = epsilon)
  iteration = iteration + 1
}

saveRDS(results_sparse_full_small, "Results_sparse_full_5_small.RDS")

#Second simulation, now set sparsity to 1
sparsity <- 1

#Define iteration as zero
iteration <- 0
results_dense_full_small <- list()

#formulate results list
for(n in N){
  name <- paste("results_", as.character(iteration), sep = "")
  results_dense_full_small[[name]] <- Simulation_time_full(n=n,i=i,dx=dx,dg=dg,factors=factors,sparsity=sparsity,unknown=unknown,lambda=lambda
                                           , f=f, epsilon = epsilon)
  iteration = iteration + 1
}

saveRDS(results_dense_full_small, "Results_dense_full_5_small.RDS")

#Earlier singular simulations----
n <- 10000
i <- 200
dx <- 100
dg <- 100
sparsity <- 0.05
unknown <- 1
lambda <- 0.5
factors <- 5
f <- 5
epsilon <- 0.005

results_3 <- Simulation_full_extra_inf(n=n,i=i,dx=dx,dg=dg,factors=factors,sparsity=sparsity,unknown=unknown,lambda=lambda
                                       , f=f, epsilon = epsilon)
plot_results_full(results_3)

#plot for full additional information
n <- 10000
i <- 200
dx <- 100
dg <- 100
sparsity <- 1
unknown <- 1
lambda <- 0.5
factors <- 5
f <- 5
epsilon <- 0.005

results_4 <- Simulation_full_extra_inf(n=n,i=i,dx=dx,dg=dg,factors=factors,sparsity=sparsity,unknown=unknown,lambda=lambda
                                       , f=f, epsilon = epsilon)
plot_results_full(results_4)


#Output plots-----
#' Retrieve difference and absolute stats from raw inputs of Simulation_time
#'
#' @param input1 "Dense" matrix
#' @param input2 "Sparse" matrix
#'
#' @return list of differences and absolute values
#' @export
#'
#' @examples
stats <- function(input1, input2){
  difference <- data.frame(N = c(10^2, round(10^(2.5)), 10^3, round(10^(3.5)), 10^4, round(10^(4.5)), 10^5),
                           old = NA, SVD = NA, Eigen = NA, Deriv = NA)
  absolute_sparse <- data.frame(N = c(10^2, round(10^(2.5)), 10^3, round(10^(3.5)), 10^4, round(10^(4.5)), 10^5),
                                old = NA, SVD = NA, Eigen = NA, Deriv = NA)
  for(i in 0:(length(input1)-1)){
    result_name <- paste("results_",as.character(i),sep = "")
    
    #Extract difference various models
    difference[(i+1),2] <- mean(input1[[result_name]]$output$Old) - mean(input2[[result_name]]$output$Old)
    difference[(i+1),3] <- mean(input1[[result_name]]$output$SVD) - mean(input2[[result_name]]$output$SVD)
    difference[(i+1),4] <- mean(input1[[result_name]]$output$Eigen) - mean(input2[[result_name]]$output$Eigen)
    difference[(i+1),5] <- mean(input1[[result_name]]$output$Deriv) - mean(input2[[result_name]]$output$Deriv)
    
    #Extract absolute time of various models
    absolute_sparse[(i+1),2] <- mean(input2[[result_name]]$output$Old)
    absolute_sparse[(i+1),3] <- mean(input2[[result_name]]$output$SVD)
    absolute_sparse[(i+1),4] <- mean(input2[[result_name]]$output$Eigen)
    absolute_sparse[(i+1),5] <- mean(input2[[result_name]]$output$Deriv)
  }
  return(list(difference = difference, absolute_sparse = absolute_sparse))
}

stats_full <- function(input1, input2){
  difference <- data.frame(N = c(10^2, round(10^(2.5)), 10^3, round(10^(3.5)), 10^4, round(10^(4.5)), 10^5),
                           old = NA, full = NA)
  absolute_sparse <- data.frame(N = c(10^2, round(10^(2.5)), 10^3, round(10^(3.5)), 10^4, round(10^(4.5)), 10^5),
                                old = NA, full = NA)
  for(i in 0:(length(input1)-1)){
    result_name <- paste("results_",as.character(i),sep = "")
    
    #Extract difference various models
    difference[(i+1),2] <- mean(input1[[result_name]]$output$Old) - mean(input2[[result_name]]$output$Old)
    difference[(i+1),3] <- mean(input1[[result_name]]$output$full) - mean(input2[[result_name]]$output$full)

    #Extract absolute time of various models
    absolute_sparse[(i+1),2] <- mean(input2[[result_name]]$output$Old)
    absolute_sparse[(i+1),3] <- mean(input2[[result_name]]$output$full)
  }
  return(list(difference = difference, absolute_sparse = absolute_sparse))
}


#Load the results of limited extra information
Results_sparse_5 <- readRDS("~/Dropbox/Uni/Master_Econometrie/Thesis/Code/R/Results_sparse_5.RDS")
Results_dense_5 <- readRDS("~/Dropbox/Uni/Master_Econometrie/Thesis/Code/R/Results_dense_5.RDS")

restricted <- stats(Results_dense_5, Results_sparse_5)

difference <- restricted$difference %>% gather(Model,Difference, -N)
absolute <- restricted$absolute_sparse %>% gather(Model,Absolute, -N)

#Graph the difference
ggplot(data = difference, aes(x = N, y = Difference, color = Model)) + geom_line(size = 1, alpha = 0.5) +
  theme_bw()+ xlab("Number of users") + ylab("Difference in convergence time in seconds")+
  scale_x_log10(breaks = c(0,10,100,1000,10000,100000))+ theme(legend.position = c(0.83, 0.67)) + scale_y_log10()

#Graph the absolute
ggplot(data = absolute, aes(x = N, y = Absolute, color = Model)) + geom_line(size = 1, alpha = 0.5) +
  theme_bw()+ xlab("Number of users") + ylab("Difference in convergence time in seconds")+scale_x_log10(breaks = c(0,10,100,1000,10000,100000))


#Load results of the full additional information plots

Results_sparse_full_5 <- readRDS("~/Dropbox/Uni/Master_Econometrie/Thesis/Code/R/Results_sparse_full_5.RDS")
Results_dense_full_5 <- readRDS("~/Dropbox/Uni/Master_Econometrie/Thesis/Code/R/Results_dense_full_5.RDS")

full <- stats_full(Results_dense_full_5, Results_sparse_full_5)

difference_full <- full$difference %>% gather(Model,Difference, -N)
absolute_full <- full$absolute_sparse %>% gather(Model,Absolute, -N)

ggplot(data = difference_full, aes(x = N, y = Difference, color = Model)) + geom_line(size = 1, alpha = 0.5) +
  theme_bw()+ xlab("Number of users") + ylab("Difference in convergence time in seconds")+
  scale_x_log10(breaks = c(0,10,100,1000,10000,100000))+ theme(legend.position = c(0.83, 0.67)) + scale_y_log10()

#Graph the absolute
ggplot(data = absolute_full, aes(x = N, y = Absolute, color = Model)) + geom_line(size = 1, alpha = 0.5) +
  theme_bw()+ xlab("Number of users") + ylab("Difference in convergence time in seconds")+scale_x_log10(breaks = c(0,10,100,1000,10000,100000))

