library(dplyr)
library(readr)
library(softImpute)
# library(spam)
library(tidyverse)
library(tictoc)
library(RcppArmadillo)
library(Rcpp)
library(ggplot2)
library(WeightIt)

sourceCpp("~/Dropbox/Uni/Master_Econometrie/Thesis/Code/R/gammaui.cpp")
sourceCpp("C:/Users/Thijs/Dropbox/Uni/Master_Econometrie/Thesis/Code/R/gammaui.cpp")

# Tools ------------------------------------------------------------------------

# Partial loglikelihood when observed value is 1
logllh1 <- function(x){
  return(log(1 + exp(-x)))
}

# Partial loglikelihood when observed value is 0
logllh0 <- function(x){
  return(x + log(1 + exp(-x)))
}

# Derivative of logllh1
derf1 <- function(x){
  return(-1 / (1 + exp(x)))
}

# Derivative of logllh0
derf2 <- function(x){
  return(1 / (1 + exp(-x)))
}

# Sigmoid function
mu <- function(x){
  return(1 / (1 + exp(-x)))
}






#Basemodel functions -------
#' Create a train/test split, given a training set
#'
#' @param df df containing indices, click, and "ratio" columns (see "1.1 Data preparation")
#' @param onlyVar logical variable specifying whether rows/columns without variation (only)
#' @param cv indicates whether train test split is made for CV
#' @param ind vector with fold indices in case of CV
#' @param fold fold that should currently be the test set
#' zeroes or ones are omitted
#'
#' @return returns a training set and the test set
trainTest <- function(df, onlyVar, cv=FALSE, ind=NULL, fold=NULL){
  # Formatting
  names(df) <- c("USERID_ind", "OFFERID_ind", "CLICK", "ratioU", "ratioO")
  
  # Make the test train split (test is 1)
  if (cv){
    # In case of cross validation (recode to zeroes and ones for ease)
    df$train_test <- 0
    df$train_test[ind == fold] <- 1
  }
  else {
    # In case of a random draws
    df$train_test <- rbinom(n = nrow(df), size = 1, prob = 0.2)
  }
  
  # Pre-allocate column for predictions
  df$prediction <- NA
  
  # Deleting the rows/columns without variation
  if (onlyVar) {
    
    # Split them (temporarily)
    df_test <- df[as.logical(df$train_test), ]
    df_train <- df[!(as.logical(df$train_test)), ]
    
    # Assign the 0 or 1 to test set obs where a ratio is 0 or 1 (prediction in advance)
    df_test$prediction[(df_test$ratioU == 0 | df_test$ratioO == 0)] <- 0
    df_test$prediction[(df_test$ratioU == 1 | df_test$ratioO == 1)] <- 1
    
    # Drop the train obs where a ratio is 0 or 1
    df_train <- df_train[!(df_train$ratioU == 0 | df_train$ratioO == 0 | 
                             df_train$ratioU == 1 | df_train$ratioO == 1), ]
    
    # Merge the two to make indices
    df <- dplyr::bind_rows(df_train, df_test)
  }
  
  # Create new indices. Make sure test is at bottom
  df <- df[order(df$train_test), ]
  df <- df %>% 
    mutate(USERID_ind_new = group_indices(., factor(USERID_ind, levels = unique(USERID_ind))))
  df <- df %>% 
    mutate(OFFERID_ind_new = group_indices(., factor(OFFERID_ind, levels = unique(OFFERID_ind))))
  
  # Split sets
  df_test <- df[as.logical(df$train_test), ]
  df_train <- df[!(as.logical(df$train_test)), c("USERID_ind_new", "OFFERID_ind_new", "CLICK", 
                                                 "ratioU", "ratioO")]
  
  # Return
  output <- list("df_train" = df_train, "df_test" = df_test)
  return(output)
}

#' Gives initial estimates for alpha, beta, C and D
#'
#' @param df training df containing ONLY indices and click
#' @param factors depth of C and D
#' @param priorsdu variance of normal distr for C
#' @param priorsdi variance of normal distr for D
#' @param initType method used for initialization (integer value)
#'
#' @return initial estimates of alpha, beta, C and D
#'
initChoose <- function(df, factors, lambda, initType, a_in = NULL, b_in = NULL,
                       C_in = NULL, D_in = NULL){
  # Formatting
  names(df) <- c("USERID_ind", "OFFERID_ind", "CLICK")
  
  nu <- max(df[ ,"USERID_ind"])
  ni <- max(df[ ,"OFFERID_ind"])
  
  # 0 for alpha and beta. Normal for C and D with mean 0
  if (initType == 1){ 
    #Alpha and beta's initialized with a zero, C and D with normal priors
    alpha <- rep(0, nu)
    beta <- rep(0, ni)
    C <- matrix(rnorm(nu * factors, 0, 1/lambda), nu, factors)
    D <- matrix(rnorm(ni * factors, 0, 1/lambda), ni, factors)
  
    # Take alpha's that create avg click rate per user
  } else if (initType == 2){
    
    # Make user click averages
    temp <- df %>%
      group_by(USERID_ind) %>%
      summarize(meanCLICK = mean(CLICK)) %>%
      select(meanCLICK)
    
    # Give some value when this is 0 or 1 (otherwise gamma -> inf)
    temp[temp == 0] <- 0.01 # THINK ABOUT THIS
    temp[temp == 1] <- 0.99 
    
    # Calculate the gamma's that produce these click rates
    alpha <- as.matrix(-1 * log(1/temp - 1))
    
    # Simple zero means for the other parameters
    beta <- rep(0, ni)
    C <- matrix(rnorm(nu * factors, 0, 1/lambda), nu, factors)
    D <- matrix(rnorm(ni * factors, 0, 1/lambda), ni, factors)
  }
  
  # If a specific input for alpha, beta C or D is given then overwrite
  if (!is.null(a_in)){
    alpha <- a_in
  }
  if (!is.null(b_in)){
    beta <- b_in
  }
  if (!is.null(C_in)){
    C <- C_in
  }
  if (!is.null(D_in)){
    D <- D_in
  }
  
  output <- list("alpha" = alpha, "beta" = beta, "C" = C, "D" = D)
  return(output)
}

#' Main algorithm for attaining alpha, beta, C and D
#'
#' @param df Dataframe consisting of userid, orderid and click. Id's should run continuously
#' from 1 to end.
#' @param factors "Width" of C and D
#' @param lambda lambda
#' @param iter Iterlation limit
#' @param epsilon Convergence criteria
#'
#' @return returns parameters alpha, beta, C and D
parEst <- function(df, factors, lambda, iter, initType, llh = TRUE, time, df_test=NULL, 
                   epsilon=NULL, a_in = NULL, b_in = NULL, C_in = NULL, D_in = NULL) {
  names(df) <- c("USERID_ind", "OFFERID_ind", "CLICK")
  
  # Initialization
  initPars <- initChoose(df, factors, lambda, initType, a_in, b_in, C_in, 
                         D_in)
  
  #Center required parameters for identification
  alpha <- initPars$alpha
  beta <- scale(initPars$beta, scale = FALSE)
  C <- scale(initPars$C, scale = FALSE)
  D <- scale(initPars$D, scale = FALSE)
  
  #Define parameters
  df <- as.matrix(df)
  logllh_old <- NULL
  
  nu <- max(df[,"USERID_ind"])
  ni <- max(df[,"OFFERID_ind"])
  
  #Retrieve indices for y=1 and y=0 from the input data
  y1 <- df[which(df[ ,"CLICK"] == 1), c("USERID_ind", "OFFERID_ind")]
  y0 <- df[which(df[ ,"CLICK"] == 0), c("USERID_ind", "OFFERID_ind")]
  
  df1 <- cbind(y1, "deriv" = NA)
  df0 <- cbind(y0, "deriv" = NA)
  
  gamma_y1 <- get_gamma0(y1[,1], y1[,2], alpha, beta, C, D)
  gamma_y0 <- get_gamma0(y0[,1], y0[,2], alpha, beta, C, D)
  
  run <- 1
  
  if (llh) {
    # Keeping track of likelihoods
    logllh_all <- rep(NA, (iter+1))
  } else {
    logllh_all <- NA
  }
  
  if(time){
    # Keeping track of time
    tic.clearlog()
  } else {
    timings <- NA
  }
  
  while (run <= iter) {
    tic("time of an iteration")
    
    if (!is.null(epsilon) || llh) {
      logllh <- sum(logllh1(gamma_y1)) + sum(logllh0(gamma_y0))
      logllh_new <- logllh +
        lambda / 2 * norm(C, type = "F") ^ 2 + lambda / 2 * norm(D, type = "F") ^ 2
    }
    
    if (llh){
      # Log Likelihood of current iteration
      logllh_all[run] <- logllh_new
    }
    
    if(run>1){
      if (!is.null(epsilon)) {
        print(logllh_new)
        print((logllh_new-logllh_old)/logllh_old)
        if (abs((logllh_new-logllh_old)/logllh_old) < epsilon) break
      }
    }
    # tic(paste("Complete iteration", run, sep = " "))
    #Define low rank representation of gamma0
    low_rankC <- cbind(C, alpha, rep(1, nu))
    low_rankD <- cbind(D, rep(1,ni), beta)
    
    #Calculate gamma0
    # gamma0 <- low_rankC %*% t(low_rankD)
    
    #Calculate respective first derivatives for both y=1 and y=0
    df1[,"deriv"] <- -4 * derf1(gamma_y1)
    df0[,"deriv"] <- -4 * derf2(gamma_y0)
    
    #Combine the results in one matrix
    df01 <- rbind(df1, df0)
    
    #Turn this matrix to sparse, notice that the dims had to be manually set
    # (for missing items probably)
    sparse <- sparseMatrix(i = (df01[ ,"USERID_ind"]), j = (df01[ ,"OFFERID_ind"]),
                           x = df01[ ,"deriv"], dims = c(nu, ni))
    
    #Calculating the H matrix for alpha update
    H_slr <- splr(sparse, low_rankC, low_rankD)
    
    #Updating alpha and beta
    newalpha <- as.matrix((1/ni) * H_slr %*% rep(1, ni))
    
    #Subtract the rowmean from H for the update for beta
    low_rankC <- cbind(C, (alpha - newalpha), rep(1, nu))
    H_slr_rowmean <- splr(sparse, low_rankC, low_rankD)
    newbeta <- as.matrix((1/nu) * t(t(rep(1, nu)) %*% H_slr_rowmean))
    
    #Updating the C and D
    #Remove row and column mean from H
    low_rankD <- cbind(D, rep(1, ni), (beta - colMeans(H_slr_rowmean)))
    H_slr_rowandcolmean <-splr(sparse, low_rankC, low_rankD)
    
    #Retrieve C and D from the svd.als function
    results <- svd.als(H_slr_rowandcolmean, rank.max = factors, lambda = lambda / 2)
    
    # Updates
    alpha <- newalpha
    beta <- newbeta
    
    # With one factor the ohter code gives an error due to d being scalar
    if (factors == 1) {
      C <- results$u %*% sqrt(results$d)
      D <- results$v %*% sqrt(results$d)
    }
    else {
      C <- results$u %*% diag(sqrt(results$d))
      D <- results$v %*% diag(sqrt(results$d))
    }
    
    # Updating gamma
    gamma_y1 <- get_gamma0(y1[,1], y1[,2], alpha, beta, C, D)
    gamma_y0 <- get_gamma0(y0[,1], y0[,2], alpha, beta, C, D)
    
    #Save old log likelihood value
    if (!is.null(epsilon)) {
      logllh_old <- logllh_new
    }
    
    toc(log = TRUE)
    run <- run + 1
  }
  
  if(time){
    log.lst <- tic.log(format = FALSE)
    timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
  }
  
  output <- list("alpha" = alpha, "beta" = beta, "C" = C, "D" = D, "logllh" = logllh_all[1:run], 
                 "run" = run, time = timings)
  return(output)
}

#' Get predictions for a test set
#'
#' @param df Test set consisting of user and offer id's and CLICK (NA or value)
#' @param alpha parameter estimate alpha
#' @param beta parameter estimate beta
#' @param C parameter estimate C
#' @param D parameter estimate D
#' @param uniqueU basically a vector used as a dictionary
#' @param uniqueI basically a vector used as a dictionary
#'
#' @return dataframe including predictions, NA for unknown user/item
getPredict <- function(df, alpha, beta, C, D){
  names(df) <- c("USERID_ind", "OFFERID_ind", "CLICK", "ratioU", "ratioO", "prediction")
  
  # By using the size of C and D, we can infer which obs are missing in training
  maxU <- nrow(C)
  maxI <- nrow(D)
  
  # Marking offer/items that are non existent in the training set
  df$nonMiss <- ((df[ ,"USERID_ind"]  <= maxU) & (df[ ,"OFFERID_ind"] <= maxI))
  
  # Predciting for the non missing obs
  # Get the non missing indices
  nonMiss <- as.matrix(df[df$nonMiss, c("USERID_ind", "OFFERID_ind")])
  
  # Calculating gamma
  gamma <- get_gamma0(nonMiss[,1], nonMiss[,2], alpha, beta, C, D)
  
  # And predict useing the llh function and gamma
  df$prediction[df$nonMiss] <- mu(gamma)
  
  return(df)
}

#' Run the full algorithm
#'
#' @param df_train
#' @param df_test 
#' @param factors 
#' @param lambda 
#' @param iter 
#' @param initType 
#' @param onlyVar 
#'
#' @return
#' @export
#'
#' @examples
fullAlg <- function(df_train, df_test, factors, lambda, iter, initType, llh=FALSE, 
                    rmse=FALSE, epsilon=NULL, a_in = NULL, b_in = NULL, C_in = NULL, D_in = NULL){
  # Estimating parameters
  #tic("2. Estimating parameters")
  pars <- parEst(df_train, factors, lambda, iter, initType, llh, rmse, df_test, 
                 epsilon, a_in, b_in, C_in, D_in)
  toc()
  
  # Getting predictions
  tic("3. Getting predictions")
  results <- getPredict(df_test[ ,c("USERID_ind_new", "OFFERID_ind_new", "CLICK",
                                    "ratioU", "ratioO", "prediction")], 
                        pars$alpha, pars$beta, pars$C, pars$D)
  toc()
  
  # RMSE
  # What to do with the NANs (some majority rule)
  results$prediction[is.na(results$prediction)] <- 0
  RMSE <- sqrt(mean((results$prediction - results$CLICK)^2))
  
  #Calculate MPR
  results <- results %>% group_by(USERID_ind) %>% mutate("perc_rank" = rank(1-prediction)/n())
  MPR <- sum(results$CLICK*results$perc_rank)/sum(results$CLICK)
  
  
  # Calculate confusion matrix
  threshold <- 0.02191325 # average click rate
  results$predictionBin <- rep(0, length(results$prediction))
  results$predictionBin[results$prediction > threshold] <- 1
  
  # True positives:
  TP <- sum(results$predictionBin == 1 & results$CLICK == 1)
  TN <- sum(results$predictionBin == 0 & results$CLICK == 0)
  FP <- sum(results$predictionBin == 1 & results$CLICK == 0)
  FN <- sum(results$predictionBin == 0 & results$CLICK == 1)
  
  confusion <- list("TP" = TP, "TN" = TN, "FP" = FP, "FN" = FN)
  
  # Output
  output <- list("parameters" = pars, "prediction" = results, "RMSE" = RMSE, "MPR" = MPR,
                 "confusion" = confusion)
  return(output)
}

crossValidate <- function(df, FACTORS, LAMBDA, INITTYPE, ONLYVAR, folds, iter, 
                          epsilon, warm){
  
  # Initialize a multidimensional output array
  # Rows are all the possible permutations of the huperparameters * folds
  rows <- (length(ONLYVAR) * length(FACTORS) * length(LAMBDA) * length(INITTYPE)) * folds
  
  # Columns for the hyperparameters, plus a name variable, and then all the results you want
  # these are: rmse, TP (true positive(1)), TN, FP, FN, number of iterations, best baseline, epsilon
  columns <- 15
  
  # Initialize the df (depth is the number of folds)
  CVoutput <- data.frame(matrix(NA, nrow = rows, ncol = columns))
  names(CVoutput) <- c("Factor", "Lambda", "InitType", "OnlyVar", "Epsilon", "Specification",
                       "RMSE", "MPR", "TP", "TN", "FP", "FN", "Iter", "rmseUser", "DifferenceRMSE")
  
  # Now we loop
  # First we make the folds
  # Randomly shuffle the data
  set.seed(123)
  df <- df[sample(nrow(df)), ]
  
  # Then assign 1-5 fold indices
  foldInd <- cut(seq(1, nrow(df)), breaks = folds, labels = FALSE)
  
  row <- 1
  # Looping over your folds
  for (z in 1:folds){
    # Do onlyvar first because the train test split depends on it
    for (a in 1:length(ONLYVAR)){
      # Do onlyvar first because the train test split depends on it
      onlyVar <- ONLYVAR[a]
      
      # Make the train test split by using the foldInd and fold as input (see trainTest)
      split <- trainTest(df, onlyVar, cv = TRUE, ind = foldInd, fold = z)
      df_train <-split$df_train[ ,c("USERID_ind_new", "OFFERID_ind_new", "CLICK")]
      df_test <- split$df_test
      
      # Loop the other hyperparameters
      
      for (b in 1:length(FACTORS)){
        # Initialize the warm start objects (can only be used within a certain factor size)
        a_in <- NULL
        b_in <- NULL
        C_in <- NULL
        D_in <- NULL
        for (c in 1:length(LAMBDA)){
          for (d in 1:length(INITTYPE)){
            tic(paste("Run", row, "out of", rows, "in fold", z, "out of ", folds))
            factors <- FACTORS[b]
            lambda <- LAMBDA[c]
            initType <- INITTYPE[d]
            
            # Run the algorithm
            output <- fullAlg(df_train, df_test, factors, lambda, iter, initType, 
                              epsilon = epsilon, a_in=a_in, b_in=b_in, C_in=C_in, D_in=D_in)
            
            # Fill the array with output
            CVoutput$Factor[row] <- factors
            CVoutput$Lambda[row] <- lambda
            CVoutput$InitType[row] <- initType
            CVoutput$OnlyVar[row] <- onlyVar
            CVoutput$Epsilon[row] <- epsilon
            
            # The name
            CVoutput$Specification[row] <- paste("Factor = ", factors, ", Lambda = ", lambda, 
                                                 ", initType = ", initType, ", onlyVar = ", onlyVar,
                                                 ", warm = ", warm, sep = "")
            
            # Performance variables
            CVoutput$RMSE[row] <- output$RMSE
            CVoutput$MPR[row] <- output$MPR
            CVoutput$TP[row] <- output$confusion$TP
            CVoutput$TN[row] <- output$confusion$TN
            CVoutput$FP[row] <- output$confusion$FP
            CVoutput$FN[row] <- output$confusion$FN
            CVoutput$Iter[row] <- output$parameters$run - 1
            CVoutput$rmseUser[row] <- baselinePred(df_train, df_test)$rmseUser
            CVoutput$DifferenceRMSE[row] <- CVoutput$RMSE[row]-CVoutput$rmseUser[row]
            
            row <- row+1
            
            # In case of warm starts, keep track of the last parameters
            if (warm){
              a_in <- output$parameters$alpha
              b_in <- output$parameters$beta
              C_in <- output$parameters$C
              D_in <- output$parameters$D
            }
            
            toc()
            
            
          }
        }
      }
    }
  }
  
  # Create a mean table
  CVmean <- CVoutput %>% 
    group_by(Epsilon, Specification) %>%
    summarise_all(mean)
  
  return(list("CVoutput" = CVoutput, "CVmean" = CVmean))
}

baselinePred <- function(df_train, df_test){
  # initialize column with majority
  df_test$predUser <- 0
  df_test$predOffer <- 0
  df_test$predOverall <- mean(df_train$CLICK)
  df_test$predMajority <- 0
  
  # Fill in predictions where available
  df_test$predUser <- df_test$ratioU[!is.na(df_test$ratioU)]
  df_test$predOffer <- df_test$ratioO[!is.na(df_test$ratioO)]
  
  rmseUser <- sqrt(mean((df_test$predUser - df_test$CLICK)^2))
  rmseOffer <- sqrt(mean((df_test$predOffer - df_test$CLICK)^2))
  rmseOverall <- sqrt(mean((df_test$predOverall - df_test$CLICK)^2))
  rmseMajority <- sqrt(mean((df_test$predMajority - df_test$CLICK)^2))
  
  output <- list("df_test" = df_test, "rmseUser" = rmseUser, "rmseOffer" = rmseOffer, 
                 "rmseOverall" = rmseOverall, "rmseMajority" = rmseMajority)
  
  return(output)
}


#Thesis functions--------
#' Create a train/test split, given a training set
#'
#' @param df df containing indices, click, and "ratio" columns (see "1.1 Data preparation")
#' @param onlyVar logical variable specifying whether rows/columns without variation (only)
#' @param cv indicates whether train test split is made for CV
#' @param ind vector with fold indices in case of CV
#' @param fold fold that should currently be the test set
#' zeroes or ones are omitted
#' @param X Data frame containing the additional information for users
#' @param G Data frame containing the additional information for items
#'
#' @return returns a training set and the test set
trainTest_thes <- function(df, X, G, onlyVar, cv=FALSE, ind=NULL, fold=NULL){
  
  # Formatting
  names(df) <- c("USERID_ind", "OFFERID_ind", "CLICK", "ratioU", "ratioO")
  
  # Make the test train split (test is 1)
  if (cv){
    # In case of cross validation (recode to zeroes and ones for ease)
    df$train_test <- 0
    df$train_test[ind == fold] <- 1
  }
  else {
    # In case of a random draws
    df$train_test <- rbinom(n = nrow(df), size = 1, prob = 0.2)
  }
  
  # Pre-allocate column for predictions
  df$prediction <- NA
  
  # Deleting the rows/columns without variation
  if (onlyVar) {
    
    # Split them (temporarily)
    df_test <- df[as.logical(df$train_test), ]
    df_train <- df[!(as.logical(df$train_test)), ]
    
    # Assign the 0 or 1 to test set obs where a ratio is 0 or 1 (prediction in advance)
    df_test$prediction[(df_test$ratioU == 0 | df_test$ratioO == 0)] <- 0
    df_test$prediction[(df_test$ratioU == 1 | df_test$ratioO == 1)] <- 1
    
    # Drop the train obs where a ratio is 0 or 1
    df_train <- df_train[!(df_train$ratioU == 0 | df_train$ratioO == 0 | 
                             df_train$ratioU == 1 | df_train$ratioO == 1), ]
    
    # Merge the two to make indices
    df <- dplyr::bind_rows(df_train, df_test)
  }
  
  # Create new keys. Make sure test is at bottom
  df <- df[order(df$train_test), ]
  
  key_user <- data.frame("USERID_ind" = unique(df$USERID_ind)) %>%
    mutate(USERID_ind_new = group_indices(., factor(USERID_ind, levels = unique(USERID_ind))))
  key_item <- data.frame("OFFERID_ind" = unique(df$OFFERID_ind)) %>% 
    mutate(OFFERID_ind_new = group_indices(., factor(OFFERID_ind, levels = unique(OFFERID_ind))))
  
  #Left join keys with matrices
  df <- df %>% left_join(key_user) %>% left_join(key_item)
  X <- X %>% left_join(key_user) %>% select(-c("USERID_ind"))
  G <- G %>% left_join(key_item) %>% select(-c("OFFERID_ind"))
  
  # Split sets
  #Keep old and new keys for test
  df_test <- df[as.logical(df$train_test), ]
  #keep only the new keys for the train, to be in lign with X and G
  df_train <- df[!(as.logical(df$train_test)), c("USERID_ind_new", "OFFERID_ind_new", "CLICK", 
                                                 "ratioU", "ratioO")]
  
  # Return
  output <- list("df_train" = df_train, "df_test" = df_test, "X" = X, "G" = G, key_user, key_item)
  return(output)
}

#' Reformat both df, X and G
#'
#' @param df 
#' @param X 
#' @param G 
#'
#' @return df_res, the restricted observations, df, the remaining variables, Coorectly partioned X and G
#' @export
#'
#' @examples
formatdfXG <- function(df_train, X, G){
  names(df_train) <- c("USERID_ind_new", "OFFERID_ind_new", "CLICK")
  
  #select only users and items which are in the training set
  X <- X[which(X$USERID_ind_new %in% df_train$USERID_ind_new),]
  G <- G[which(G$OFFERID_ind_new %in% df_train$OFFERID_ind_new),]
  
  #define matrix characteristics
  nu <- max(df_train$USERID_ind_new)
  ni <- max(df_train$OFFERID_ind_new)
  nx <- nrow(X)
  ng <- nrow(G)
  
  #create empty list for the keys
  keys <- list()
  
  #Create two keys, one with user and one with item indexes
  key_user <-  data.frame(USERID_ind_new = unique(df_train$USERID_ind_new))
  key_item <-  data.frame(OFFERID_ind_new = unique(df_train$OFFERID_ind_new))
  
  #recode the X and the G so that the indices of the users and items in the X and G are first
  X_newind <- X %>% 
    mutate(USERID_recoded = group_indices(., factor(USERID_ind_new, levels = unique(USERID_ind_new))))
  G_newind <- G %>% 
    mutate(OFFERID_recoded = group_indices(., factor(OFFERID_ind_new, levels = unique(OFFERID_ind_new))))
  
  #Merge with the corresponding keys and sort on the new keys 
  key_user <- left_join(key_user, select(X_newind, c("USERID_ind_new", "USERID_recoded"))) %>% arrange(USERID_recoded)
  key_item <- left_join(key_item, select(G_newind, c("OFFERID_ind_new", "OFFERID_recoded"))) %>% arrange(OFFERID_recoded)
  
  #fill the remaining keys
  key_user[(nx+1):nu,"USERID_recoded"]<- (nx+1):nu
  key_item[(ng+1):ni,"OFFERID_recoded"]<- (ng+1):ni
  
  #Recode and rename X, G and df_train with the new keys
  df_train <- df_train %>% left_join(key_user)%>%left_join(key_item) %>% select(-c("USERID_ind_new", "OFFERID_ind_new")) %>%
    rename(USERID_ind_new = USERID_recoded) %>% rename(OFFERID_ind_new = OFFERID_recoded) %>% 
    select(c("USERID_ind_new", "OFFERID_ind_new", "CLICK"))
  X <- X %>% left_join(key_user) %>%select(-c("USERID_ind_new")) %>% rename(USERID_ind_new = USERID_recoded)
  G <- G %>% left_join(key_item) %>% select(-c("OFFERID_ind_new")) %>% rename(OFFERID_ind_new = OFFERID_recoded)
  
  #remove the user and item numbers from X as they are not considered variables, additionally, make full rank
  X <- select(X, -c("USERID_ind_new")) %>% make_full_rank()
  G <- select(G, -c("OFFERID_ind_new")) %>% make_full_rank()
  
  keys$key_user <- key_user
  keys$key_item <- key_item
  
  return(list("df_train" = df_train, "X" = X, "G" = G, "keys" = keys))
}

recode_results <- function(pars, model, keys){
  #Extract the keys
  key_user <- keys$key_user
  key_item <- keys$key_item
  
  #Resort the keys, to reshufle the results to the old keys
  key_user <- key_user[order(key_user$USERID_ind_new), ]
  key_item <- key_item[order(key_item$OFFERID_ind_new), ]
  
  #Two different recoding criteria, one for full restriction and one for partial restriction
  if(model == "Full"){
    alpha <- pars$alpha
    beta <- pars$beta
    C_pr <- pars$XC_r
    D_pr <- pars$GD_r
  }
  else{
    alpha <- pars$alpha
    beta <- pars$beta
    C_pr <- pars$C_pr
    D_pr <- pars$D_pr
  }
  
  #Return correct order of variables
  return(list("alpha" = alpha[key_user$USERID_recoded], "beta" = beta[key_item$OFFERID_recoded], 
         "C_pr" = C_pr[key_user$USERID_recoded,], "D_pr" = D_pr[key_item$OFFERID_recoded,]))
}

#' Gives initial estimates for alpha, beta, C and D
#'
#' @param df training df containing ONLY indices and click
#' @param factors depth of C and D
#' @param priorsdu variance of normal distr for C
#' @param priorsdi variance of normal distr for D
#' @param initType method used for initialization (integer value)
#'
#' @return initial estimates of alpha, beta, C_r and D_r
#'
initChoose_rest <- function(df, X, G, factors, lambda, a_in = NULL, b_in = NULL,
                            C_r_in = NULL, D_r_in = NULL, C_in = NULL, D_in =NULL){
  nu <- max(df[ ,"USERID_ind"])
  ni <- max(df[ ,"OFFERID_ind"])
  nx <- nrow(X)
  ng <- nrow(G)
  d_x <- ncol(X)
  d_g <- ncol(G)
  
  #click averages
  temp <- df %>%
    group_by(USERID_ind) %>%
    summarize(meanCLICK = mean(CLICK)) %>%
    select(meanCLICK)
  
  # Give some value when this is 0 or 1 (otherwise gamma -> inf)
  temp[temp == 0] <- 0.01 # THINK ABOUT THIS
  temp[temp == 1] <- 0.99 
  
  # Calculate the gamma's that produce these click rates
  alpha <- as.matrix(-1 * log(1/temp - 1))
  
  # Simple zero means for the other parameters
  beta <- rep(0, ni)
  
  XC_r <- matrix(rnorm(nx * factors, 0, 1/lambda), nx, factors)
  GD_r <- matrix(rnorm(ng * factors, 0, 1/lambda), ng, factors)
  C <- matrix(rnorm((nu-nx) * factors, 0, 1/lambda), nu-nx, factors)
  D <- matrix(rnorm((ni-ng) * factors, 0, 1/lambda), ni-ng, factors)
  C_r <- NULL
  D_r <- NULL
  
  # If a specific input for alpha, beta C_r or D_r is given then overwrite
  if (!is.null(a_in)){
    alpha <- a_in
  }
  if (!is.null(b_in)){
    beta <- b_in
  }
  if (!is.null(C_r_in)){
    C_r <- C_r_in
  }
  if (!is.null(D_r_in)){
    D_r <- D_r_in
  }
  if (!is.null(C_in)){
    C <- C_in
  }
  if (!is.null(D_in)){
    D <- D_in
  }
  
  output <- list("alpha" = alpha, "beta" = beta, "C_r" = C_r, "D_r" = D_r, "C" = C, "D" = D, "XC_r" = XC_r, "GD_r" = GD_r)
  return(output)
}

#' Parameter estimation with a restricted solution
#'
#' @param df Clicks of the user
#' @param user Additional user data
#' @param item Additional item data
#' @param factors Amount of factors in the SVD
#' @param lambda Lambda used
#' @param iter Amount of iterations until convergence
#' @param iter_CD Max iterations for finding C and D
#' @param epsilon stop algorithm when log likelihood converged to epsilon (or max iter)
#'
#' @return
#' @export
#'
#' @examples
ParEst_rest_fullrest <- function(df, X, G, factors, lambda, iter, iter_CD, epsilon, llh, time, a_in = NULL, b_in = NULL
                                 , C_r_in = NULL, D_r_in = NULL){
  names(df) <- c("USERID_ind", "OFFERID_ind", "CLICK")
  
  # Initialization of parameters
  initPars <- initChoose_rest(df, X, G, factors, lambda, a_in, b_in, C_r_in, D_r_in)
  
  #Center required parameters for identification
  alpha <- initPars$alpha
  beta <- scale(initPars$beta, scale = FALSE)
  X <- scale(X, scale = FALSE)
  G <- scale(G, scale = FALSE)
  if(is.null(C_r_in)){
    XC_r <- scale(initPars$XC_r, scale = F)
    GD_r <- scale(initPars$GD_r, scale = F)
    C_r <- initPars$C_r
    D_r <- initPars$D_r
  }
  else{
    C_r <- initPars$C_r
    D_r <- initPars$D_r
    XC_r <- X%*%C_r
    GD_r <- G%*%D_r
  }
  
  df <- as.matrix(df)
  
  #Initialization of parameters
  nu <- max(df[,"USERID_ind"])
  ni <- max(df[,"OFFERID_ind"])
  nx <- nrow(X)
  ng <- nrow(G)
  dx <- ncol(X)
  dg <- ncol(G)
  logllh_old <- NULL
  
  #some memmory stored usefull matrices
  GprimeG <- t(G)%*%G
  XprimeX <- t(X)%*%X
  XprimXminhalf <- expm::sqrtm(Matrix::solve(XprimeX))
  GprimGminhalf <- expm::sqrtm(Matrix::solve(GprimeG))
  XprimXhalf <- expm::sqrtm(XprimeX)
  GprimGhalf <- expm::sqrtm(GprimeG)
  
  #Calculate Z, the sparse matrix
  #Retrieve indices for y=1 and y=0 from the input data
  y1 <- df[which(df[ ,"CLICK"] == 1), c("USERID_ind", "OFFERID_ind")]
  y0 <- df[which(df[ ,"CLICK"] == 0), c("USERID_ind", "OFFERID_ind")]
  
  df1 <- cbind(y1, "deriv" = NA)
  df0 <- cbind(y0, "deriv" = NA)
  
  if (llh) {
    # Keeping track of likelihoods
    logllh_all <- rep(NA, (iter+1))
  } else {
    logllh_all <- NA
  }
  
  if(time){
    # Keeping track of time
    tic.clearlog()
  } else {
    timings <- NA
  }
  
  run <- 1
  while(run<=iter){
    tic("time of an iteration")
    #define usefull matrixes
    if(run>1){
      XC_r <- X%*%C_r
      GD_r <- G%*%D_r
    }
    
    #Retrieve the gamma only for y1 and y0
    gamma_y1 <- get_gamma0(y1[,1], y1[,2], alpha, beta, XC_r, GD_r)
    gamma_y0 <- get_gamma0(y0[,1], y0[,2], alpha, beta, XC_r, GD_r)
    
    #Print log likelihood
    logllh <- sum(logllh1(gamma_y1)) + sum(logllh0(gamma_y0))
    logllh_new <- logllh +
      (lambda / 2) * norm(XC_r, type = "F") ^ 2 + lambda / 2 * norm(GD_r, type = "F") ^ 2
    print(logllh_new)
    #check if the algorithm has reached the convergence criteria epsilon
    if(run>1){
      if(!is.null(epsilon)){
        decrease <- abs((logllh_new-logllh_old)/logllh_old)
        cat("Decrease in logllh: ", decrease, "\n")
        if( decrease < epsilon){
          break
        }
      }
    }
    
    if (llh) {
      # Calculate log likelihood
      logllh_all[run] <- logllh
    }
    
    #Calculate respective first derivatives for both y=1 and y=0
    df1[,"deriv"] <- -4 * derf1(gamma_y1)
    df0[,"deriv"] <- -4 * derf2(gamma_y0)
    
    #Combine the results in one matrix
    df01 <- rbind(df1, df0)
    
    #Turn this matrix to sparse, notice that the dims had to be manually set
    # (for missing items probably)
    Z <- sparseMatrix(i = (df01[ ,"USERID_ind"]), j = (df01[ ,"OFFERID_ind"]),
                      x = df01[ ,"deriv"], dims = c(nu, ni))
    
    #Lowrank part of sparse + lowrank
    low_rankC <- cbind(XC_r, alpha, rep(1, nu))
    low_rankD <- cbind(GD_r, rep(1,ni), beta)
    
    #Calculating the H matrix for alpha update
    H_slr <- splr(Z, low_rankC, low_rankD)
    
    #Updating alpha and beta
    newalpha <- as.matrix((1/ni) * H_slr %*% rep(1, ni))
    
    #Subtract the rowmean from H for the update for beta
    low_rankC <- cbind(XC_r, (alpha - rowMeans(H_slr)), rep(1, nu))
    H_slr_rowmean <- splr(Z, low_rankC, low_rankD)
    newbeta <- as.matrix((1/nu) * t(t(rep(1, nu)) %*% H_slr_rowmean))
    
    #Updating C_r and D_r
    #Define variables of first iteration
    V <- matrix(rnorm(dg*factors), dg, factors)
    V <- svd(V)$u
    U <- matrix(0L, nrow = dx, ncol = factors)
    Phi <- diag(factors, nrow = factors, ncol = factors)
    
    V_prev <- NULL
    U_prev <- NULL
    Phi_prev <- NULL
    
    #Calculate common matrices for efficiency
    GZX <- GprimGminhalf%*%t(G)%*%t(Z)%*%X%*%XprimXminhalf
    if(run>1){
      GGDCXX <- GprimGhalf%*%D_r%*%t(C_r)%*%XprimXhalf
    }
    else{
      GGDCXX <- GprimGminhalf%*%t(G)%*%GD_r%*%t(XC_r)%*%X%*%XprimXminhalf
    }
    
    #test if objective function decreases
    #low_rankD <- cbind(GD_r, rep(1,ni), beta - colMeans(H_slr_rowmean))
    #cat("regularized majorized value equation (47) + regularization: ", (1/8)*norm(Z, type = "F")^2 +(lambda/2)*(norm(XC_r, type = "F")^2+norm(GD_r, type = "F")^2), "\n")
    #H_test <- Z + low_rankC%*%t(low_rankD)
    
    #Initialize convergence criterion matrices:
    iteration_CD <- 0
    while(iteration_CD<iter_CD){
      #cat("Iter: ", iteration_CD, "\n")
      #Update U
      Uhatsvd <- svd(Phi%*%t(V)%*%GZX + Phi%*%t(V)%*%GGDCXX)
      U <- Uhatsvd$v%*%t(Uhatsvd$u)
      
      #Update V
      Vhatsvd <- svd(GZX%*%U%*%Phi + GGDCXX%*%U%*%Phi)
      V <- Vhatsvd$u%*%t(Vhatsvd$v)
      
      #Update Phi
      Ftilde <- t(V)%*%GZX%*%U + t(V)%*%GGDCXX%*%U
      for (i in 1:nrow(Phi)){
        Phi[i,i] <- max((Ftilde[i,i]-4*lambda),0)
      }
      
      #Check the exit criteria
      if(iteration_CD>1){
        criteria <- abs(((norm(V,type="F")+norm(U,type="F")+norm(Phi,type="F"))-(norm(V_prev,type="F")+norm(U_prev,type="F")+norm(Phi_prev,type="F")))/(norm(V_prev,type="F")+norm(U_prev,type="F")+norm(Phi_prev,type="F")))
        if(criteria<0.001){
          #cat("Stopped at iteration: ", iteration_CD, "\n")
          break
        }
      }
      
      iteration_CD = iteration_CD+1
      
      V_prev <- V
      U_prev <- U
      Phi_prev <- Phi
      
    }
    
    #save old likelihood value
    if(!is.null(epsilon)){
      logllh_old <- logllh_new
    }
    
    alpha <- newalpha
    beta <- newbeta
    C_r <- XprimXminhalf%*%U%*%sqrt(Phi)
    D_r <- GprimGminhalf%*%V%*%sqrt(Phi)
    
    run <- run+1
    toc(log = TRUE)
  }
  
  if(time){
    log.lst <- tic.log(format = FALSE)
    timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
  }
  
  output <- list("alpha" = alpha, "beta" = beta,"XC_r" = XC_r, "GD_r" = GD_r, "C_r" = C_r, "D_r" = D_r, 
                 "logllh" = logllh_all[1:run], "time" = timings, "run" = run)
  
  return(output)
}


ParEst_rest_parrest_SVD <- function(df, X, G, factors, lambda, iter, iter_CD, epsilon, llh, time, a_in = NULL, b_in= NULL, 
                                    C_r_in= NULL, D_r_in= NULL, D_in= NULL, C_in= NULL){
  names(df) <- c("USERID_ind", "OFFERID_ind", "CLICK")
  
  # Initialization of parameters
  initPars <- initChoose_rest(df, X, G, factors, lambda, a_in, b_in, C_r_in, D_r_in, C_in, D_in)
  
  #Center required parameters for identification
  alpha <- initPars$alpha
  beta <- scale(initPars$beta, scale = FALSE)
  X <- scale(X, scale = FALSE)
  G <- scale(G, scale = FALSE)
  
  #If no warm start, extract variables
  if(is.null(C_r_in)){
    XC_r <- scale(initPars$XC_r, scale = F)
    GD_r <- scale(initPars$GD_r, scale = F)
    C <- scale(initPars$C, scale = F)
    D <- scale(initPars$D, scale = F)
    C_r <- initPars$C_r
    D_r <- initPars$D_r
  }
  #if warm start, extract old variables
  else{
    C_r <- initPars$C_r
    D_r <- initPars$D_r
    XC_r <- X%*%C_r
    GD_r <- G%*%D_r
    C <- initPars$C
    D <- initPars$D
  }
  
  df <- as.matrix(df)
  #Initialization of parameters
  nu <- max(df[,"USERID_ind"])
  ni <- max(df[,"OFFERID_ind"])
  nx <- nrow(X)
  ng <- nrow(G)
  dx <- ncol(X)
  dg <- ncol(G)
  logllh_old <- NULL
  
  #some memmory stored usefull matrices
  GprimeG <- t(G)%*%G
  XprimeX <- t(X)%*%X
  
  GprimeGhalf <- expm::sqrtm(t(G)%*%G)
  XprimeXhalf <- expm::sqrtm(t(X)%*%X)
  
  GprimeGinverse <- solve(GprimeG)
  XprimeXinverse <- solve(XprimeX)
  
  GprimeGminhalf <- expm::sqrtm(GprimeGinverse)
  XprimeXminhalf <- expm::sqrtm(XprimeXinverse)
  
  GprimGminhalfGprim <- GprimeGminhalf%*%t(G)
  XXprimXminhalf <- X%*%XprimeXminhalf
  
  #Retrieve indices for y=1 and y=0 from the input data
  y1 <- df[which(df[ ,"CLICK"] == 1), c("USERID_ind", "OFFERID_ind")]
  y0 <- df[which(df[ ,"CLICK"] == 0), c("USERID_ind", "OFFERID_ind")]
  
  #Split y1 and y0 into the four different blocks, also correct the indices to be in line with a separate block structure
  y1_11 <- y1[which(y1[,"USERID_ind"] %in% 1:nx & y1[,"OFFERID_ind"] %in% 1:ng),]
  y1_12 <- y1[which(y1[,"USERID_ind"] %in% 1:nx & !(y1[,"OFFERID_ind"] %in% 1:ng)),]
  y1_12[,2] <- y1_12[,2]-ng
  y1_21 <- y1[which(!(y1[,"USERID_ind"] %in% 1:nx) & y1[,"OFFERID_ind"] %in% 1:ng),]
  y1_21[,1] <- y1_21[,1]-nx
  y1_22 <- y1[which(!(y1[,"USERID_ind"] %in% 1:nx) & !(y1[,"OFFERID_ind"] %in% 1:ng)),]
  y1_22[,2] <- y1_22[,2]-ng
  y1_22[,1] <- y1_22[,1]-nx
  
  y0_11 <- y0[which(y0[,"USERID_ind"] %in% 1:nx & y0[,"OFFERID_ind"] %in% 1:ng),]
  y0_12 <- y0[which(y0[,"USERID_ind"] %in% 1:nx & !(y0[,"OFFERID_ind"] %in% 1:ng)),]
  y0_12[,2] <- y0_12[,2]-ng
  y0_21 <- y0[which(!(y0[,"USERID_ind"] %in% 1:nx) & y0[,"OFFERID_ind"] %in% 1:ng),]
  y0_21[,1] <- y0_21[,1]-nx
  y0_22 <- y0[which(!(y0[,"USERID_ind"] %in% 1:nx) & !(y0[,"OFFERID_ind"] %in% 1:ng)),]
  y0_22[,2] <- y0_22[,2]-ng
  y0_22[,1] <- y0_22[,1]-nx
  
  df1_11 <- cbind(y1_11, "deriv" = NA)
  df1_12 <- cbind(y1_12, "deriv" = NA)
  df1_21 <- cbind(y1_21, "deriv" = NA)
  df1_22 <- cbind(y1_22, "deriv" = NA)
  
  df0_11 <- cbind(y0_11, "deriv" = NA)
  df0_12 <- cbind(y0_12, "deriv" = NA)
  df0_21 <- cbind(y0_21, "deriv" = NA)
  df0_22 <- cbind(y0_22, "deriv" = NA)
  
  if (llh) {
    # Keeping track of likelihoods
    logllh_all <- rep(NA, (iter+1))
  } else {
    logllh_all <- NA
  }
  
  if(time){
    # Keeping track of time
    tic.clearlog()
  } else {
    timings <- NA
  }
  
  
  run <- 1
  while(run<=iter){
    #define usefull matrixes
    tic("time of an iteration")
    if(run == 1){
      C_r <- NULL
      D_r <- NULL
    }
    else{
      XC_r <- X%*%C_r
      GD_r <- G%*%D_r}
    
    C_pr <- rbind(XC_r, C)
    D_pr <- rbind(GD_r, D)
    
    #define subvectors of alpha and beta
    alpha_1 <- alpha[1:nx]
    alpha_2 <- alpha[(nx+1):nu]
    beta_1 <- beta[1:ng]
    beta_2 <- beta[(ng+1):ni]
    
    #Retrieve the gamma only for y1 and y0
    
    gamma_y1_11 <- get_gamma0(y1_11[,1], y1_11[,2], alpha_1, beta_1, XC_r, GD_r)
    gamma_y1_12 <- get_gamma0(y1_12[,1], y1_12[,2], alpha_1, beta_2, XC_r, D)
    gamma_y1_21 <- get_gamma0(y1_21[,1], y1_21[,2], alpha_2, beta_1, C, GD_r)
    gamma_y1_22 <- get_gamma0(y1_22[,1], y1_22[,2], alpha_2, beta_2, C, D)
    
    gamma_y0_11 <- get_gamma0(y0_11[,1], y0_11[,2], alpha_1, beta_1, XC_r, GD_r)
    gamma_y0_12 <- get_gamma0(y0_12[,1], y0_12[,2], alpha_1, beta_2, XC_r, D)
    gamma_y0_21 <- get_gamma0(y0_21[,1], y0_21[,2], alpha_2, beta_1, C, GD_r)
    gamma_y0_22 <- get_gamma0(y0_22[,1], y0_22[,2], alpha_2, beta_2, C, D)
    
    #Print log likelihood
    logllh <- sum(logllh1(gamma_y1_11)) + sum(logllh1(gamma_y1_12))+ sum(logllh1(gamma_y1_21))+ sum(logllh1(gamma_y1_22)) +
      sum(logllh0(gamma_y0_11))+ sum(logllh0(gamma_y0_21))+ sum(logllh0(gamma_y0_12))+ sum(logllh0(gamma_y0_22))
    logllh_new <- logllh +
      (lambda / 2) * norm(C_pr, type = "F") ^ 2 + lambda / 2 * norm(D_pr, type = "F") ^ 2
    cat("Log likelihood value iteration ", run,": ", logllh_new, "\n")
    
    #check if the algorithm has reached the convergence criteria epsilon
    if(run>1){
      if(!is.null(epsilon)){
        change <- abs((logllh_new-logllh_old)/logllh_old)
        cat("*****Change in logllh: ", change,"*****","\n")
        if( abs(change) < epsilon){
          break
        }
      }
    }
    
    if (llh) {
      # Calculate log likelihood
      logllh_all[run] <- logllh
    }
    
    #Calculate respective first derivatives for both y=1 and y=0
    df1_11[,"deriv"] <- -4 * derf1(gamma_y1_11)
    df1_12[,"deriv"] <- -4 * derf1(gamma_y1_12)
    df1_21[,"deriv"] <- -4 * derf1(gamma_y1_21)
    df1_22[,"deriv"] <- -4 * derf1(gamma_y1_22)
    
    df0_11[,"deriv"] <- -4 * derf2(gamma_y0_11)
    df0_12[,"deriv"] <- -4 * derf2(gamma_y0_12)
    df0_21[,"deriv"] <- -4 * derf2(gamma_y0_21)
    df0_22[,"deriv"] <- -4 * derf2(gamma_y0_22)
    
    #Combine the results in one matrix
    df01_11 <- rbind(df1_11, df0_11)
    df01_12 <- rbind(df1_12, df0_12)
    df01_21 <- rbind(df1_21, df0_21)
    df01_22 <- rbind(df1_22, df0_22)
    
    #Define the four different sparse matrices
    Z_11 <- sparseMatrix(i = (df01_11[ ,"USERID_ind"]), j = (df01_11[ ,"OFFERID_ind"]),
                         x = df01_11[ ,"deriv"], dims = c(nx, ng))
    Z_12 <- sparseMatrix(i = (df01_12[ ,"USERID_ind"]), j = (df01_12[ ,"OFFERID_ind"]),
                         x = df01_12[ ,"deriv"], dims = c(nx, (ni-ng)))
    Z_21 <- sparseMatrix(i = (df01_21[ ,"USERID_ind"]), j = (df01_21[ ,"OFFERID_ind"]),
                         x = df01_21[ ,"deriv"], dims = c((nu-nx), ng))
    Z_22 <- sparseMatrix(i = (df01_22[ ,"USERID_ind"]), j = (df01_22[ ,"OFFERID_ind"]),
                         x = df01_22[ ,"deriv"], dims = c((nu-nx), (ni-ng)))
    
    #Define the four different low rank matrices
    low_rankC_11 <- cbind(XC_r, alpha_1, rep(1, nx))
    low_rankD_11 <- cbind(GD_r, rep(1,ng), beta_1)
    low_rankC_12 <- low_rankC_11
    low_rankD_12 <- cbind(D, rep(1,(ni-ng)), beta_2)
    low_rankC_21 <- cbind(C, alpha_2, rep(1, (nu-nx)))
    low_rankD_21 <- low_rankD_11
    low_rankC_22 <- low_rankC_21
    low_rankD_22 <- low_rankD_12
    
    #Combine the four sparse + low rank data sets
    H_slr_11 <- splr(Z_11, low_rankC_11, low_rankD_11)
    H_slr_12 <- splr(Z_12, low_rankC_12, low_rankD_12)
    H_slr_21 <- splr(Z_21, low_rankC_21, low_rankD_21)
    H_slr_22 <- splr(Z_22, low_rankC_22, low_rankD_22)
    
    #Updating alpha
    newalpha_1 <- (ng*rowMeans(H_slr_11)+(ni-ng)*rowMeans(H_slr_12))/ni
    newalpha_2 <- (ng*rowMeans(H_slr_21)+(ni-ng)*rowMeans(H_slr_22))/ni
    
    #subtract rowmeans from H and formulate beta updates
    low_rankC_rowmean_11 <- cbind(XC_r, alpha_1-rowMeans(H_slr_11), rep(1, nx))
    low_rankC_rowmean_12 <- cbind(XC_r, alpha_1-rowMeans(H_slr_12), rep(1, nx))
    low_rankC_rowmean_21 <- cbind(C, alpha_2-rowMeans(H_slr_21), rep(1, (nu-nx)))
    low_rankC_rowmean_22 <- cbind(C, alpha_2-rowMeans(H_slr_22), rep(1, (nu-nx)))
    
    H_slr_rowmean_11 <- splr(Z_11, low_rankC_rowmean_11, low_rankD_11)
    H_slr_rowmean_12 <- splr(Z_12, low_rankC_rowmean_12, low_rankD_12)
    H_slr_rowmean_21 <- splr(Z_21, low_rankC_rowmean_21, low_rankD_21)
    H_slr_rowmean_22 <- splr(Z_22, low_rankC_rowmean_22, low_rankD_22)
    
    newbeta_1 <- (nx*colMeans(H_slr_rowmean_11)+(nu-nx)*colMeans(H_slr_rowmean_21))/nu
    newbeta_2 <- (nx*colMeans(H_slr_rowmean_12)+(nu-nx)*colMeans(H_slr_rowmean_22))/nu
    
    #Start updating for C_pr and D_pr
    #Define variables of first iteration
    V_pr <- matrix(rnorm((dg+(ni-ng))*factors), dg+(ni-ng), factors)
    V_pr <- svd(V_pr)$u
    U_pr <- matrix(0L, nrow = dx+(nu-nx), ncol = factors)
    Phi <- diag(factors, nrow = factors, ncol = factors)
    
    V_pr_prev <- NULL
    U_pr_prev <- NULL
    Phi_prev <- NULL
    #Split V into V_r and V
    V_r = V_pr[(1:dg),(1:ncol(V_pr))]
    V = V_pr[((dg+1):nrow(V_pr)),(1:ncol(V_pr))]
    
    #Split U into U_r and U
    U_r = U_pr[(1:dx),(1:ncol(U_pr))]
    U = U_pr[((dx+1):nrow(U_pr)),(1:ncol(U_pr))]
    
    #Calculate matrices invariant in the iterations 
    Xi_11 <- GprimGminhalfGprim%*%t(Z_11)%*%XXprimXminhalf + GprimGminhalfGprim%*%GD_r%*%t(XC_r)%*%XXprimXminhalf
    Xi_12 <- t(scale(t(GprimGminhalfGprim%*%t(Z_21)), scale = FALSE))+GprimGminhalfGprim%*%GD_r%*%t(C)
    Xi_21 <- scale(t(Z_12)%*%XXprimXminhalf, scale = FALSE)+D%*%t(XC_r)%*%XXprimXminhalf
    Xi_22_withoutZ <- D%*%t(C)
    
    #Initialize convergence criterion matrices:
    iteration_CD <- 0
    while(iteration_CD<iter_CD){
      #Update U
      if(iteration_CD == 0){
        VJZJ <- t(scale(t(t(scale(V, scale = F))%*%t(Z_22)), scale = F))
        Uhatsvd <- svd(Phi%*%cbind(t(V_r)%*%Xi_11+t(V)%*%Xi_21, t(V_r)%*%Xi_12 + VJZJ + t(V)%*%Xi_22_withoutZ))
        fixed1 <- 0
        fixed2 <- 0
      }
      else{
        Uhatsvd <- svd(Phi%*%cbind(fixed1, fixed2))
      }
      U_pr <- Uhatsvd$v%*%t(Uhatsvd$u)
      
      #Split U into U_r and U
      U_r = U_pr[(1:dx),(1:ncol(U_pr))]
      U = U_pr[((dx+1):nrow(U_pr)),(1:ncol(U_pr))]
      
      #Update V
      JZJU <- scale(t(Z_22)%*%scale(U, scale = F),scale = F)
      Vhatsvd <- svd(rbind(Xi_11%*%U_r + Xi_12%*%U, Xi_21%*%U_r + JZJU + Xi_22_withoutZ%*%U))
      V_pr <- Vhatsvd$u%*%t(Vhatsvd$v)
      
      #Split V into V_r and V
      V_r = V_pr[(1:dg),(1:ncol(V_pr))]
      V = V_pr[((dg+1):nrow(V_pr)),(1:ncol(V_pr))]
      
      #Update Phi
      VJZJ <- t(scale(t(t(scale(V, scale = F))%*%t(Z_22)), scale = F))
      fixed1 <- (t(V_r)%*%Xi_11+t(V)%*%Xi_21)
      fixed2 <- (t(V_r)%*%Xi_12 + VJZJ + t(V)%*%Xi_22_withoutZ)
      F_tilde <- fixed1%*%U_r + fixed2%*%U
      
      for (i in 1:nrow(Phi)){
        Phi[i,i] <- max((F_tilde[i,i]-4*lambda),0)
      }
      
      #Check the exit criteria
      if(iteration_CD>1){
        criteria <- abs(((norm(V_pr,type="F")+norm(U_pr,type="F")+norm(Phi,type="F"))-(norm(V_pr_prev,type="F")+norm(U_pr_prev,type="F")+norm(Phi_prev,type="F")))/(norm(V_pr_prev,type="F")+norm(U_pr_prev,type="F")+norm(Phi_prev,type="F")))
        if(criteria<0.001){
          #cat("Stopped at iteration: ", iteration_CD, "\n")
          break
        }
      }
      
      iteration_CD = iteration_CD+1
      
      V_pr_prev <- V_pr
      U_pr_prev <- U_pr
      Phi_prev <- Phi
    }
    
    #save old likelihood value
    if(!is.null(epsilon)){
      logllh_old <- logllh_new
    }
    
    alpha <- rbind(as.matrix(newalpha_1), as.matrix(newalpha_2))
    beta <- rbind(as.matrix(newbeta_1), as.matrix(newbeta_2))
    
    C_r <- as.matrix(XprimeXminhalf%*%U_r%*%sqrt(Phi))
    D_r <- as.matrix(GprimeGminhalf%*%V_r%*%sqrt(Phi))
    C <- as.matrix(U%*%sqrt(Phi))
    D <- as.matrix(V%*%sqrt(Phi))
    
    run <- run+1
    toc(log = TRUE)
  }
  
  if(time){
    log.lst <- tic.log(format = FALSE)
    timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
  }
  output <- list("alpha" = alpha, "beta" = beta, "C_pr" = C_pr, "D_pr" = D_pr,"C_r" = C_r, "D_r" = D_r,
                 "C" = C, "D" = D, "logllh" = logllh_all[1:run], "time" = timings, "run" = run)
  
  return(output)
}


#' Parameter estimation with a restricted solution
#'
#' @param df Clicks of the user
#' @param user Additional user data
#' @param item Additional item data
#' @param factors Amount of factors in the SVD
#' @param lambda Lambda used
#' @param iter Amount of iterations until convergence
#' @param iter_CD Max iterations for finding C and D
#' @param epsilon stop algorithm when log likelihood converged to epsilon (or max iter)
#'
#' @return
#' @export
#'
#' @examples
ParEst_rest_parrest_Eigen <- function(df, X, G, factors, lambda, iter, iter_CD, epsilon, llh, time, a_in = NULL, b_in= NULL, 
                                      C_r_in= NULL, D_r_in= NULL, D_in= NULL, C_in= NULL){
  names(df) <- c("USERID_ind", "OFFERID_ind", "CLICK")
  
  # Initialization of parameters
  initPars <- initChoose_rest(df, X, G, factors, lambda, a_in, b_in, C_r_in, D_r_in, C_in, D_in)
  
  #Center required parameters for identification
  alpha <- initPars$alpha
  beta <- scale(initPars$beta, scale = FALSE)
  X <- scale(X, scale = FALSE)
  G <- scale(G, scale = FALSE)
  
  #If no warm start, extract variables
  if(is.null(C_r_in)){
    XC_r <- scale(initPars$XC_r, scale = F)
    GD_r <- scale(initPars$GD_r, scale = F)
    C <- scale(initPars$C, scale = F)
    D <- scale(initPars$D, scale = F)
    C_r <- initPars$C_r
    D_r <- initPars$D_r
  }
  #if warm start, extract old variables
  else{
    C_r <- initPars$C_r
    D_r <- initPars$D_r
    XC_r <- X%*%C_r
    GD_r <- G%*%D_r
    C <- initPars$C
    D <- initPars$D
  }
  
  df <- as.matrix(df)
  #Initialization of parameters
  nu <- max(df[,"USERID_ind"])
  ni <- max(df[,"OFFERID_ind"])
  nx <- nrow(X)
  ng <- nrow(G)
  dx <- ncol(X)
  dg <- ncol(G)
  logllh_old <- NULL
  
  #Initialize both C_r and D_r seperately (needed in Algorithm 5) as empty matrix
  C_r = matrix(0, nrow = dx, ncol = factors)
  D_r = matrix(0, nrow = dg, ncol = factors)
  
  #some memmory stored usefull matrices
  GprimeG <- t(G)%*%G
  XprimeX <- t(X)%*%X
  
  GprimeGhalf <- expm::sqrtm(t(G)%*%G)
  XprimeXhalf <- expm::sqrtm(t(X)%*%X)
  
  GprimeGinverse <- solve(GprimeG)
  XprimeXinverse <- solve(XprimeX)
  
  GprimeGminhalf <- expm::sqrtm(GprimeGinverse)
  XprimeXminhalf <- expm::sqrtm(XprimeXinverse)
  
  GprimGminhalfGprim <- GprimeGminhalf%*%t(G)
  XXprimXminhalf <- X%*%XprimeXminhalf
  
  #Retrieve indices for y=1 and y=0 from the input data
  y1 <- df[which(df[ ,"CLICK"] == 1), c("USERID_ind", "OFFERID_ind")]
  y0 <- df[which(df[ ,"CLICK"] == 0), c("USERID_ind", "OFFERID_ind")]
  
  #Split y1 and y0 into the four different blocks, also correct the indices to be in line with a separate block structure
  y1_11 <- y1[which(y1[,"USERID_ind"] %in% 1:nx & y1[,"OFFERID_ind"] %in% 1:ng),]
  y1_12 <- y1[which(y1[,"USERID_ind"] %in% 1:nx & !(y1[,"OFFERID_ind"] %in% 1:ng)),]
  y1_12[,2] <- y1_12[,2]-ng
  y1_21 <- y1[which(!(y1[,"USERID_ind"] %in% 1:nx) & y1[,"OFFERID_ind"] %in% 1:ng),]
  y1_21[,1] <- y1_21[,1]-nx
  y1_22 <- y1[which(!(y1[,"USERID_ind"] %in% 1:nx) & !(y1[,"OFFERID_ind"] %in% 1:ng)),]
  y1_22[,2] <- y1_22[,2]-ng
  y1_22[,1] <- y1_22[,1]-nx
  
  y0_11 <- y0[which(y0[,"USERID_ind"] %in% 1:nx & y0[,"OFFERID_ind"] %in% 1:ng),]
  y0_12 <- y0[which(y0[,"USERID_ind"] %in% 1:nx & !(y0[,"OFFERID_ind"] %in% 1:ng)),]
  y0_12[,2] <- y0_12[,2]-ng
  y0_21 <- y0[which(!(y0[,"USERID_ind"] %in% 1:nx) & y0[,"OFFERID_ind"] %in% 1:ng),]
  y0_21[,1] <- y0_21[,1]-nx
  y0_22 <- y0[which(!(y0[,"USERID_ind"] %in% 1:nx) & !(y0[,"OFFERID_ind"] %in% 1:ng)),]
  y0_22[,2] <- y0_22[,2]-ng
  y0_22[,1] <- y0_22[,1]-nx
  
  df1_11 <- cbind(y1_11, "deriv" = NA)
  df1_12 <- cbind(y1_12, "deriv" = NA)
  df1_21 <- cbind(y1_21, "deriv" = NA)
  df1_22 <- cbind(y1_22, "deriv" = NA)
  
  df0_11 <- cbind(y0_11, "deriv" = NA)
  df0_12 <- cbind(y0_12, "deriv" = NA)
  df0_21 <- cbind(y0_21, "deriv" = NA)
  df0_22 <- cbind(y0_22, "deriv" = NA)
  
  if (llh) {
    # Keeping track of likelihoods
    logllh_all <- rep(NA, (iter+1))
  } else {
    logllh_all <- NA
  }
  
  if(time){
    # Keeping track of time
    tic.clearlog()
  } else {
    timings <- NA
  }
  
  
  run <- 1
  while(run<=iter){
    #define usefull matrixes
    tic("time of an iteration")
    if(run == 1){
      C_r <- matrix(0, dx, factors)
      D_r <- matrix(0, dg, factors)
    }
    else{
      XC_r <- X%*%C_r
      GD_r <- G%*%D_r}
    
    C_pr <- rbind(XC_r, C)
    D_pr <- rbind(GD_r, D)
    
    #define subvectors of alpha and beta
    alpha_1 <- alpha[1:nx]
    alpha_2 <- alpha[(nx+1):nu]
    beta_1 <- beta[1:ng]
    beta_2 <- beta[(ng+1):ni]
    
    #Retrieve the gamma only for y1 and y0
    
    gamma_y1_11 <- get_gamma0(y1_11[,1], y1_11[,2], alpha_1, beta_1, XC_r, GD_r)
    gamma_y1_12 <- get_gamma0(y1_12[,1], y1_12[,2], alpha_1, beta_2, XC_r, D)
    gamma_y1_21 <- get_gamma0(y1_21[,1], y1_21[,2], alpha_2, beta_1, C, GD_r)
    gamma_y1_22 <- get_gamma0(y1_22[,1], y1_22[,2], alpha_2, beta_2, C, D)
    
    gamma_y0_11 <- get_gamma0(y0_11[,1], y0_11[,2], alpha_1, beta_1, XC_r, GD_r)
    gamma_y0_12 <- get_gamma0(y0_12[,1], y0_12[,2], alpha_1, beta_2, XC_r, D)
    gamma_y0_21 <- get_gamma0(y0_21[,1], y0_21[,2], alpha_2, beta_1, C, GD_r)
    gamma_y0_22 <- get_gamma0(y0_22[,1], y0_22[,2], alpha_2, beta_2, C, D)
    
    #Print log likelihood
    logllh <- sum(logllh1(gamma_y1_11)) + sum(logllh1(gamma_y1_12))+ sum(logllh1(gamma_y1_21))+ sum(logllh1(gamma_y1_22)) +
      sum(logllh0(gamma_y0_11))+ sum(logllh0(gamma_y0_21))+ sum(logllh0(gamma_y0_12))+ sum(logllh0(gamma_y0_22))
    logllh_new <- logllh +
      (lambda / 2) * norm(C_pr, type = "F") ^ 2 + lambda / 2 * norm(D_pr, type = "F") ^ 2
    cat("Log likelihood value iteration ", run,": ", logllh_new, "\n")
    #check if the algorithm has reached the convergence criteria epsilon
    if(run>1){
      if(!is.null(epsilon)){
        change <- abs((logllh_new-logllh_old)/logllh_old)
        cat("*****Change in logllh: ", change,"*****","\n")
        if( abs(change) < epsilon){
          break
        }
      }
    }
    
    if (llh) {
      # Calculate log likelihood
      logllh_all[run] <- logllh
    }
    
    #Calculate respective first derivatives for both y=1 and y=0
    df1_11[,"deriv"] <- -4 * derf1(gamma_y1_11)
    df1_12[,"deriv"] <- -4 * derf1(gamma_y1_12)
    df1_21[,"deriv"] <- -4 * derf1(gamma_y1_21)
    df1_22[,"deriv"] <- -4 * derf1(gamma_y1_22)
    
    df0_11[,"deriv"] <- -4 * derf2(gamma_y0_11)
    df0_12[,"deriv"] <- -4 * derf2(gamma_y0_12)
    df0_21[,"deriv"] <- -4 * derf2(gamma_y0_21)
    df0_22[,"deriv"] <- -4 * derf2(gamma_y0_22)
    
    #Combine the results in one matrix
    df01_11 <- rbind(df1_11, df0_11)
    df01_12 <- rbind(df1_12, df0_12)
    df01_21 <- rbind(df1_21, df0_21)
    df01_22 <- rbind(df1_22, df0_22)
    
    #Define the four different sparse matrices
    Z_11 <- sparseMatrix(i = (df01_11[ ,"USERID_ind"]), j = (df01_11[ ,"OFFERID_ind"]),
                         x = df01_11[ ,"deriv"], dims = c(nx, ng))
    Z_12 <- sparseMatrix(i = (df01_12[ ,"USERID_ind"]), j = (df01_12[ ,"OFFERID_ind"]),
                         x = df01_12[ ,"deriv"], dims = c(nx, (ni-ng)))
    Z_21 <- sparseMatrix(i = (df01_21[ ,"USERID_ind"]), j = (df01_21[ ,"OFFERID_ind"]),
                         x = df01_21[ ,"deriv"], dims = c((nu-nx), ng))
    Z_22 <- sparseMatrix(i = (df01_22[ ,"USERID_ind"]), j = (df01_22[ ,"OFFERID_ind"]),
                         x = df01_22[ ,"deriv"], dims = c((nu-nx), (ni-ng)))
    
    #Define the four different low rank matrices
    low_rankC_11 <- cbind(XC_r, alpha_1, rep(1, nx))
    low_rankD_11 <- cbind(GD_r, rep(1,ng), beta_1)
    low_rankC_12 <- low_rankC_11
    low_rankD_12 <- cbind(D, rep(1,(ni-ng)), beta_2)
    low_rankC_21 <- cbind(C, alpha_2, rep(1, (nu-nx)))
    low_rankD_21 <- low_rankD_11
    low_rankC_22 <- low_rankC_21
    low_rankD_22 <- low_rankD_12
    
    #Combine the four sparse + low rank data sets
    H_slr_11 <- splr(Z_11, low_rankC_11, low_rankD_11)
    H_slr_12 <- splr(Z_12, low_rankC_12, low_rankD_12)
    H_slr_21 <- splr(Z_21, low_rankC_21, low_rankD_21)
    H_slr_22 <- splr(Z_22, low_rankC_22, low_rankD_22)
    
    #Updating alpha
    newalpha_1 <- (ng*rowMeans(H_slr_11)+(ni-ng)*rowMeans(H_slr_12))/ni
    newalpha_2 <- (ng*rowMeans(H_slr_21)+(ni-ng)*rowMeans(H_slr_22))/ni
    
    #subtract rowmeans from H and formulate beta updates
    low_rankC_rowmean_11 <- cbind(XC_r, alpha_1-rowMeans(H_slr_11), rep(1, nx))
    low_rankC_rowmean_12 <- cbind(XC_r, alpha_1-rowMeans(H_slr_12), rep(1, nx))
    low_rankC_rowmean_21 <- cbind(C, alpha_2-rowMeans(H_slr_21), rep(1, (nu-nx)))
    low_rankC_rowmean_22 <- cbind(C, alpha_2-rowMeans(H_slr_22), rep(1, (nu-nx)))
    
    H_slr_rowmean_11 <- splr(Z_11, low_rankC_rowmean_11, low_rankD_11)
    H_slr_rowmean_12 <- splr(Z_12, low_rankC_rowmean_12, low_rankD_12)
    H_slr_rowmean_21 <- splr(Z_21, low_rankC_rowmean_21, low_rankD_21)
    H_slr_rowmean_22 <- splr(Z_22, low_rankC_rowmean_22, low_rankD_22)
    
    newbeta_1 <- (nx*colMeans(H_slr_rowmean_11)+(nu-nx)*colMeans(H_slr_rowmean_21))/nu
    newbeta_2 <- (nx*colMeans(H_slr_rowmean_12)+(nu-nx)*colMeans(H_slr_rowmean_22))/nu
    
    #Update C_pr and D_pr
    
    #Define variables of first iteration
    V_pr <- matrix(rnorm((dg+(ni-ng))*factors), dg+(ni-ng), factors)
    V_pr <- svd(V_pr)$u
    U_pr <- matrix(0L, nrow = dx+(nu-nx), ncol = factors)
    Phi <- diag(factors, nrow = factors, ncol = factors)
    
    V_pr_prev <- NULL
    U_pr_prev <- NULL
    Phi_prev <- NULL
    #Split V into V_r and V
    V_r = V_pr[(1:dg),(1:ncol(V_pr))]
    V = V_pr[((dg+1):nrow(V_pr)),(1:ncol(V_pr))]
    
    #Split U into U_r and U
    U_r = U_pr[(1:dx),(1:ncol(U_pr))]
    U = U_pr[((dx+1):nrow(U_pr)),(1:ncol(U_pr))]
    
    #Calculate matrices invariant in the iterations 
    Xi_11 <- GprimGminhalfGprim%*%t(Z_11)%*%XXprimXminhalf + GprimGminhalfGprim%*%GD_r%*%t(XC_r)%*%XXprimXminhalf
    Xi_12 <- t(scale(t(GprimGminhalfGprim%*%t(Z_21)), scale = FALSE))+GprimGminhalfGprim%*%GD_r%*%t(C)
    Xi_21 <- scale(t(Z_12)%*%XXprimXminhalf, scale = FALSE)+D%*%t(XC_r)%*%XXprimXminhalf
    Xi_22_withoutZ <- D%*%t(C)
    
    #Initialize convergence criterion matrices:
    iteration_CD <- 0
    while(iteration_CD<iter_CD){
      #Update U
      if(iteration_CD == 0){
        VJZJ <- t(scale(t(t(scale(V, scale = F))%*%t(Z_22)), scale = F))
        F_1 <- Phi%*%cbind(t(V_r)%*%Xi_11+t(V)%*%Xi_21, t(V_r)%*%Xi_12 + VJZJ + t(V)%*%Xi_22_withoutZ)
        fixed1 <- 0
        fixed2 <- 0
      }
      else{
        F_1 <- Phi%*%cbind(fixed1, fixed2)
      }
      F_1F_1trans <- F_1%*%t(F_1)
      Uhateigen <- eigen(F_1F_1trans)
      Q <- solve(sqrt(diag(Uhateigen$values)))%*%t(Uhateigen$vectors)%*%F_1
      U_pr <- t(Q)%*%t(Uhateigen$vectors)
      
      #Split U into U_r and U
      U_r = U_pr[(1:dx),(1:ncol(U_pr))]
      U = U_pr[((dx+1):nrow(U_pr)),(1:ncol(U_pr))]
      
      #Update V
      JZJU <- scale(t(Z_22)%*%scale(U, scale = F), scale = F)
      F_2 <- rbind(Xi_11%*%U_r + Xi_12%*%U, Xi_21%*%U_r + JZJU + Xi_22_withoutZ%*%U)
      F_2transF_2 <- t(F_2)%*%F_2
      Vhateigen <- eigen(F_2transF_2)
      P <- F_2%*%Vhateigen$vectors%*%solve(sqrt(diag(Vhateigen$values)))
      V_pr <- P%*%t(Vhateigen$vectors)
      
      #Split V into V_r and V
      V_r = V_pr[(1:dg),(1:ncol(V_pr))]
      V = V_pr[((dg+1):nrow(V_pr)),(1:ncol(V_pr))]
      
      #Update Phi
      VJZJ <- t(scale(t(t(scale(V, scale = F))%*%t(Z_22)), scale = F))
      fixed1 <- (t(V_r)%*%Xi_11+t(V)%*%Xi_21)
      fixed2 <- (t(V_r)%*%Xi_12 + VJZJ + t(V)%*%Xi_22_withoutZ)
      F_tilde <- fixed1%*%U_r + fixed2%*%U
      
      for (i in 1:nrow(Phi)){
        Phi[i,i] <- max((F_tilde[i,i]-4*lambda),0)
      }
      
      #Check the exit criteria
      if(iteration_CD>1){
        criteria <- abs(((norm(V_pr,type="F")+norm(U_pr,type="F")+norm(Phi,type="F"))-(norm(V_pr_prev,type="F")+norm(U_pr_prev,type="F")+norm(Phi_prev,type="F")))/(norm(V_pr_prev,type="F")+norm(U_pr_prev,type="F")+norm(Phi_prev,type="F")))
        if(criteria<0.001){
          #cat("Stopped at iteration: ", iteration_CD, "\n")
          break
        }
      }
      
      iteration_CD = iteration_CD+1
      
      V_pr_prev <- V_pr
      U_pr_prev <- U_pr
      Phi_prev <- Phi
    }
    
    #save old likelihood value
    if(!is.null(epsilon)){
      logllh_old <- logllh_new
    }
    
    alpha <- rbind(as.matrix(newalpha_1), as.matrix(newalpha_2))
    beta <- rbind(as.matrix(newbeta_1), as.matrix(newbeta_2))
    
    C_r <- as.matrix(XprimeXminhalf%*%U_r%*%sqrt(Phi))
    D_r <- as.matrix(GprimeGminhalf%*%V_r%*%sqrt(Phi))
    C <- as.matrix(U%*%sqrt(Phi))
    D <- as.matrix(V%*%sqrt(Phi))
    
    run <- run+1
    toc(log = TRUE)
  }
  
  if(time){
    log.lst <- tic.log(format = FALSE)
    timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
  }
  output <- list("alpha" = alpha, "beta" = beta, "C_pr" = C_pr, "D_pr" = D_pr,"C_r" = C_r, "D_r" = D_r,
                 "C" = C, "D" = D, "logllh" = logllh_all[1:run], "time" = timings, "run" = run)
  return(output)
}


#' Parameter estimation with a restricted solution
#'
#' @param df Clicks of the user
#' @param user Additional user data
#' @param item Additional item data
#' @param factors Amount of factors in the SVD
#' @param lambda Lambda used
#' @param iter Amount of iterations until convergence
#' @param iter_CD Max iterations for finding C and D
#' @param epsilon stop algorithm when log likelihood converged to epsilon (or max iter)
#'
#' @return
#' @export
#'
#' @examples
ParEst_rest_parrest_deriv <- function(df, X, G, factors, lambda, iter, iter_CD, epsilon, llh, time, a_in = NULL, b_in= NULL, 
                                      C_r_in= NULL, D_r_in= NULL, D_in= NULL, C_in= NULL){
  names(df) <- c("USERID_ind", "OFFERID_ind", "CLICK")
  
  # Initialization of parameters
  initPars <- initChoose_rest(df, X, G, factors, lambda, a_in, b_in, C_r_in, D_r_in, C_in, D_in)
  
  #Center required parameters for identification
  alpha <- initPars$alpha
  beta <- scale(initPars$beta, scale = FALSE)
  X <- scale(X, scale = FALSE)
  G <- scale(G, scale = FALSE)
  
  #If no warm start, extract variables
  if(is.null(C_r_in)){
    XC_r <- scale(initPars$XC_r, scale = F)
    GD_r <- scale(initPars$GD_r, scale = F)
    C <- scale(initPars$C, scale = F)
    D <- scale(initPars$D, scale = F)
    C_r <- initPars$C_r
    D_r <- initPars$D_r
  }
  #if warm start, extract old variables
  else{
    C_r <- initPars$C_r
    D_r <- initPars$D_r
    XC_r <- X%*%C_r
    GD_r <- G%*%D_r
    C <- initPars$C
    D <- initPars$D
  }
  
  df <- as.matrix(df)
  #Initialization of parameters
  nu <- max(df[,"USERID_ind"])
  ni <- max(df[,"OFFERID_ind"])
  nx <- nrow(X)
  ng <- nrow(G)
  dx <- ncol(X)
  dg <- ncol(G)
  logllh_old <- NULL
  
  #some memmory stored usefull matrices
  GprimeG <- t(G)%*%G
  XprimeX <- t(X)%*%X
  
  GprimeGinverse <- solve(GprimeG)
  XprimeXinverse <- solve(XprimeX)
  #Retrieve indices for y=1 and y=0 from the input data
  y1 <- df[which(df[ ,"CLICK"] == 1), c("USERID_ind", "OFFERID_ind")]
  y0 <- df[which(df[ ,"CLICK"] == 0), c("USERID_ind", "OFFERID_ind")]
  
  #Split y1 and y0 into the four different blocks, also correct the indices to be in line with a separate block structure
  y1_11 <- y1[which(y1[,"USERID_ind"] %in% 1:nx & y1[,"OFFERID_ind"] %in% 1:ng),]
  y1_12 <- y1[which(y1[,"USERID_ind"] %in% 1:nx & !(y1[,"OFFERID_ind"] %in% 1:ng)),]
  y1_12[,2] <- y1_12[,2]-ng
  y1_21 <- y1[which(!(y1[,"USERID_ind"] %in% 1:nx) & y1[,"OFFERID_ind"] %in% 1:ng),]
  y1_21[,1] <- y1_21[,1]-nx
  y1_22 <- y1[which(!(y1[,"USERID_ind"] %in% 1:nx) & !(y1[,"OFFERID_ind"] %in% 1:ng)),]
  y1_22[,2] <- y1_22[,2]-ng
  y1_22[,1] <- y1_22[,1]-nx
  
  y0_11 <- y0[which(y0[,"USERID_ind"] %in% 1:nx & y0[,"OFFERID_ind"] %in% 1:ng),]
  y0_12 <- y0[which(y0[,"USERID_ind"] %in% 1:nx & !(y0[,"OFFERID_ind"] %in% 1:ng)),]
  y0_12[,2] <- y0_12[,2]-ng
  y0_21 <- y0[which(!(y0[,"USERID_ind"] %in% 1:nx) & y0[,"OFFERID_ind"] %in% 1:ng),]
  y0_21[,1] <- y0_21[,1]-nx
  y0_22 <- y0[which(!(y0[,"USERID_ind"] %in% 1:nx) & !(y0[,"OFFERID_ind"] %in% 1:ng)),]
  y0_22[,2] <- y0_22[,2]-ng
  y0_22[,1] <- y0_22[,1]-nx
  
  df1_11 <- cbind(y1_11, "deriv" = NA)
  df1_12 <- cbind(y1_12, "deriv" = NA)
  df1_21 <- cbind(y1_21, "deriv" = NA)
  df1_22 <- cbind(y1_22, "deriv" = NA)
  
  df0_11 <- cbind(y0_11, "deriv" = NA)
  df0_12 <- cbind(y0_12, "deriv" = NA)
  df0_21 <- cbind(y0_21, "deriv" = NA)
  df0_22 <- cbind(y0_22, "deriv" = NA)
  
  if (llh) {
    # Keeping track of likelihoods
    logllh_all <- rep(NA, (iter+1))
  } else {
    logllh_all <- NA
  }
  
  if(time){
    # Keeping track of time
    tic.clearlog()
  } else {
    timings <- NA
  }
  
  run <- 1
  while(run<=iter){
    #define usefull matrixes
    tic("time of an iteration")
    if(run == 1){
      C_r <- matrix(0, dx, factors)
      D_r <- matrix(0, dg, factors)
    }
    else{
      XC_r <- X%*%C_r
      GD_r <- G%*%D_r}
    
    C_pr <- rbind(XC_r, C)
    D_pr <- rbind(GD_r, D)
    
    #define subvectors of alpha and beta
    alpha_1 <- alpha[1:nx]
    alpha_2 <- alpha[(nx+1):nu]
    beta_1 <- beta[1:ng]
    beta_2 <- beta[(ng+1):ni]
    
    #Retrieve the gamma only for y1 and y0
    
    gamma_y1_11 <- get_gamma0(y1_11[,1], y1_11[,2], alpha_1, beta_1, XC_r, GD_r)
    gamma_y1_12 <- get_gamma0(y1_12[,1], y1_12[,2], alpha_1, beta_2, XC_r, D)
    gamma_y1_21 <- get_gamma0(y1_21[,1], y1_21[,2], alpha_2, beta_1, C, GD_r)
    gamma_y1_22 <- get_gamma0(y1_22[,1], y1_22[,2], alpha_2, beta_2, C, D)
    
    gamma_y0_11 <- get_gamma0(y0_11[,1], y0_11[,2], alpha_1, beta_1, XC_r, GD_r)
    gamma_y0_12 <- get_gamma0(y0_12[,1], y0_12[,2], alpha_1, beta_2, XC_r, D)
    gamma_y0_21 <- get_gamma0(y0_21[,1], y0_21[,2], alpha_2, beta_1, C, GD_r)
    gamma_y0_22 <- get_gamma0(y0_22[,1], y0_22[,2], alpha_2, beta_2, C, D)
    
    #Print log likelihood
    logllh <- sum(logllh1(gamma_y1_11)) + sum(logllh1(gamma_y1_12))+ sum(logllh1(gamma_y1_21))+ sum(logllh1(gamma_y1_22)) +
      sum(logllh0(gamma_y0_11))+ sum(logllh0(gamma_y0_21))+ sum(logllh0(gamma_y0_12))+ sum(logllh0(gamma_y0_22))
    logllh_new <- logllh +
      (lambda / 2) * norm(C_pr, type = "F") ^ 2 + lambda / 2 * norm(D_pr, type = "F") ^ 2
    cat("Log likelihood value iteration ", run,": ", logllh_new, "\n")
    #check if the algorithm has reached the convergence criteria epsilon
    if(run>1){
      if(!is.null(epsilon)){
        change <- abs((logllh_new-logllh_old)/logllh_old)
        cat("*****Change in logllh: ", change,"*****","\n")
        if( abs(change) < epsilon){
          break
        }
      }
    }
    
    if (llh) {
      # Calculate log likelihood
      logllh_all[run] <- logllh
    }
    
    #Calculate respective first derivatives for both y=1 and y=0
    df1_11[,"deriv"] <- -4 * derf1(gamma_y1_11)
    df1_12[,"deriv"] <- -4 * derf1(gamma_y1_12)
    df1_21[,"deriv"] <- -4 * derf1(gamma_y1_21)
    df1_22[,"deriv"] <- -4 * derf1(gamma_y1_22)
    
    df0_11[,"deriv"] <- -4 * derf2(gamma_y0_11)
    df0_12[,"deriv"] <- -4 * derf2(gamma_y0_12)
    df0_21[,"deriv"] <- -4 * derf2(gamma_y0_21)
    df0_22[,"deriv"] <- -4 * derf2(gamma_y0_22)
    
    #Combine the results in one matrix
    df01_11 <- rbind(df1_11, df0_11)
    df01_12 <- rbind(df1_12, df0_12)
    df01_21 <- rbind(df1_21, df0_21)
    df01_22 <- rbind(df1_22, df0_22)
    
    #Define the four different sparse matrices
    Z_11 <- sparseMatrix(i = (df01_11[ ,"USERID_ind"]), j = (df01_11[ ,"OFFERID_ind"]),
                         x = df01_11[ ,"deriv"], dims = c(nx, ng))
    Z_12 <- sparseMatrix(i = (df01_12[ ,"USERID_ind"]), j = (df01_12[ ,"OFFERID_ind"]),
                         x = df01_12[ ,"deriv"], dims = c(nx, (ni-ng)))
    Z_21 <- sparseMatrix(i = (df01_21[ ,"USERID_ind"]), j = (df01_21[ ,"OFFERID_ind"]),
                         x = df01_21[ ,"deriv"], dims = c((nu-nx), ng))
    Z_22 <- sparseMatrix(i = (df01_22[ ,"USERID_ind"]), j = (df01_22[ ,"OFFERID_ind"]),
                         x = df01_22[ ,"deriv"], dims = c((nu-nx), (ni-ng)))
    
    #Define the four different low rank matrices
    low_rankC_11 <- cbind(XC_r, alpha_1, rep(1, nx))
    low_rankD_11 <- cbind(GD_r, rep(1,ng), beta_1)
    low_rankC_12 <- low_rankC_11
    low_rankD_12 <- cbind(D, rep(1,(ni-ng)), beta_2)
    low_rankC_21 <- cbind(C, alpha_2, rep(1, (nu-nx)))
    low_rankD_21 <- low_rankD_11
    low_rankC_22 <- low_rankC_21
    low_rankD_22 <- low_rankD_12
    
    #Combine the four sparse + low rank data sets
    H_slr_11 <- splr(Z_11, low_rankC_11, low_rankD_11)
    H_slr_12 <- splr(Z_12, low_rankC_12, low_rankD_12)
    H_slr_21 <- splr(Z_21, low_rankC_21, low_rankD_21)
    H_slr_22 <- splr(Z_22, low_rankC_22, low_rankD_22)
    
    #Updating alpha
    newalpha_1 <- (ng*rowMeans(H_slr_11)+(ni-ng)*rowMeans(H_slr_12))/ni
    newalpha_2 <- (ng*rowMeans(H_slr_21)+(ni-ng)*rowMeans(H_slr_22))/ni
    
    #subtract rowmeans from H and formulate beta updates
    low_rankC_rowmean_11 <- cbind(XC_r, alpha_1-rowMeans(H_slr_11), rep(1, nx))
    low_rankC_rowmean_12 <- cbind(XC_r, alpha_1-rowMeans(H_slr_12), rep(1, nx))
    low_rankC_rowmean_21 <- cbind(C, alpha_2-rowMeans(H_slr_21), rep(1, (nu-nx)))
    low_rankC_rowmean_22 <- cbind(C, alpha_2-rowMeans(H_slr_22), rep(1, (nu-nx)))
    
    H_slr_rowmean_11 <- splr(Z_11, low_rankC_rowmean_11, low_rankD_11)
    H_slr_rowmean_12 <- splr(Z_12, low_rankC_rowmean_12, low_rankD_12)
    H_slr_rowmean_21 <- splr(Z_21, low_rankC_rowmean_21, low_rankD_21)
    H_slr_rowmean_22 <- splr(Z_22, low_rankC_rowmean_22, low_rankD_22)
    
    newbeta_1 <- (nx*colMeans(H_slr_rowmean_11)+(nu-nx)*colMeans(H_slr_rowmean_21))/nu
    newbeta_2 <- (nx*colMeans(H_slr_rowmean_12)+(nu-nx)*colMeans(H_slr_rowmean_22))/nu
    
    #Update C_r, D_r, C and D
    #Revert the indices of df01 back to the original indices
    df01_12[,2] <- df01_12[,2]+ng
    df01_21[,1] <- df01_21[,1]+nx
    df01_22[,1] <- df01_22[,1]+nx
    df01_22[,2] <- df01_22[,2]+ng
    
    #Define the Z_sparse matrices needed for the updates
    df01_C_r <- rbind(df01_11, df01_12)
    df01_D_r <- rbind(df01_11, df01_21)
    df01_C <- rbind(df01_21, df01_22)
    df01_D <- rbind(df01_12, df01_22)
    
    Z_C_r <- sparseMatrix(i = (df01_C_r[ ,"USERID_ind"]), j = (df01_C_r[ ,"OFFERID_ind"]),
                          x = df01_C_r[ ,"deriv"], dims = c(nx, ni))
    Z_D_r <- sparseMatrix(i = (df01_D_r[ ,"USERID_ind"]), j = (df01_D_r[ ,"OFFERID_ind"]),
                          x = df01_D_r[ ,"deriv"], dims = c(nu, ng))
    Z_C <- sparseMatrix(i = (df01_C[ ,"USERID_ind"])-nx, j = (df01_C[ ,"OFFERID_ind"]),
                        x = df01_C[ ,"deriv"], dims = c((nu-nx), ni))
    Z_D <- sparseMatrix(i = (df01_D[ ,"USERID_ind"]), j = (df01_D[ ,"OFFERID_ind"])-ng,
                        x = df01_D[ ,"deriv"], dims = c(nu, (ni-ng)))
    
    
    #Calculate matrices invariant in the iterations 
    fixedC_r <- XprimeXinverse%*%t(X)%*%Z_C_r  #XprimeXinverse%*%(t(X)%*%XC_r)%*%(t(D_pr)%*%D_pr
    fixedD_r <- GprimeGinverse%*%t(G)%*%t(Z_D_r)
    
    #Define matrices from the current run
    C_pr_old <- C_pr
    D_pr_old <- D_pr
    C_old <- C
    D_old <- D
    XC_r_old <- XC_r
    GD_r_old <- GD_r
    
    #Initialize convergence criterion matrices:
    iteration_CD <- 0
    while(iteration_CD<iter_CD){
      #Define parameters for convergence
      C_pr_prev <- C_pr
      D_pr_prev <- D_pr
      
      #Update C_r
      C_r <- (fixedC_r%*%D_pr + XprimeXinverse%*%t(X)%*%XC_r_old%*%(t(D_pr_old)%*%D_pr))%*%solve(t(D_pr)%*%D_pr+4*lambda*diag(factors))
      #update params with C_r
      XC_r <- X%*%C_r
      C_pr <- rbind(XC_r, C)
      
      #similar for D_r
      D_r <- (fixedD_r%*%C_pr + GprimeGinverse%*%t(G)%*%GD_r_old%*%(t(C_pr_old)%*%C_pr))%*%solve(t(C_pr)%*%C_pr+4*lambda*diag(factors))
      GD_r <- G%*%D_r
      D_pr <- rbind(GD_r, D)
      
      #similar for C
      C <- (scale(Z_C%*%D_pr, scale = F) +scale(alpha_2 %*% t(colSums(D_pr)), scale = F) + C_old%*%(t(D_pr_old)%*%D_pr))%*%solve(t(D_pr)%*%D_pr+4*lambda*diag(factors))
      C_pr <- rbind(XC_r, C)
      
      #silmilar for D
      D <- (scale(t(Z_D)%*%C_pr, scale = F) + scale(beta_2 %*% t(colSums(C_pr)), scale = F)+ D_old%*%(t(C_pr_old)%*%C_pr))%*%solve(t(C_pr)%*%C_pr+4*lambda*diag(factors))
      D_pr <- rbind(GD_r, D)
      
      #Check the exit criteria
      if(iteration_CD>1){
        criteria <- abs(((norm(C_pr,type="F")+norm(D_pr,type="F"))-(norm(C_pr_prev,type="F")+norm(D_pr_prev,type="F")))/(norm(C_pr_prev,type="F")+norm(D_pr_prev,type="F")))
        if(criteria<0.01){
          #cat("Stopped at iteration: ", iteration_CD, "\n")
          break
        }
      }
      iteration_CD = iteration_CD+1
      
    }
    
    #save old likelihood value
    if(!is.null(epsilon)){
      logllh_old <- logllh_new
    }
    
    #Update the parameters
    alpha <- rbind(as.matrix(newalpha_1), as.matrix(newalpha_2))
    beta <- rbind(as.matrix(newbeta_1), as.matrix(newbeta_2))
    
    C_r <- as.matrix(C_r)
    D_r <- as.matrix(D_r)
    C <- as.matrix(C)
    D <- as.matrix(D)
    
    run <- run+1
    toc(log = TRUE)
  }
  if(time){
    log.lst <- tic.log(format = FALSE)
    timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
  }
  output <- list("alpha" = alpha, "beta" = beta, "C_pr" = C_pr, "D_pr" = D_pr,"C_r" = C_r, "D_r" = D_r,
                 "C" = C, "D" = D, "logllh" = logllh_all[1:run], "time" = timings, "run" = run)
  return(output)
}

getPredict_thes <- function(df_test, alpha, beta, C_pr, D_pr){
  names(df_test) <- c("USERID_ind", "OFFERID_ind", "CLICK", "ratioU", "ratioO", "prediction")
  
  # By using the size of C and D, we can infer which obs are missing in training
  maxU <- nrow(C_pr)
  maxI <- nrow(D_pr)
  
  # Marking offer/items that are non existent in the training set
  df_test$nonMiss <- ((df_test[ ,"USERID_ind"]  <= maxU) & (df_test[ ,"OFFERID_ind"] <= maxI))
  
  # Predciting for the non missing obs
  # Get the non missing indices
  nonMiss <- as.matrix(df_test[df_test$nonMiss, c("USERID_ind", "OFFERID_ind")])
  
  # Calculating gamma
  gamma <- get_gamma0(nonMiss[,1], nonMiss[,2], alpha, beta, C_pr, D_pr)
  
  # And predict using the llh function and gamma
  df_test$prediction[df_test$nonMiss] <- mu(gamma)
  
  return(df_test)
}

fullAlg_rest <- function(df_train, X, G, df_test, factors, lambda, iter, iter_CD, initType, method_name, llh = F, 
                   time = F, epsilon=NULL, a_in = NULL, b_in = NULL, C_in = NULL, D_in = NULL, 
                   C_r_in = NULL, D_r_in = NULL){
  # Format the train, X and G matrices
  tic("1. Format the train, X and G matrices")
  prelim <- formatdfXG(df_train, X, G)
  df_train <- prelim$df_train[ ,c("USERID_ind_new", "OFFERID_ind_new", "CLICK")]
  X <- prelim$X
  G <- prelim$G
  
  # Estimating parameters
  if(method_name == "Full"){
    tic("2. Estimating parameters")
    pars <- ParEst_rest_fullrest(df_train, X, G, factors, lambda, iter, iter_CD, epsilon, llh, time, a_in, 
                                  b_in , C_r_in, D_r_in)
    toc()
  }
  else if(method_name == "SVD"){
    tic("2. Estimating parameters")
    pars <- ParEst_rest_parrest_SVD(df_train, X, G, factors, lambda, iter, iter_CD, epsilon, llh, time, a_in, b_in, 
                         C_r_in, D_r_in, D_in, C_in)
    toc()
  }
  else if(method_name == "Eigen"){
    tic("2. Estimating parameters")
    pars <- ParEst_rest_parrest_Eigen(df_train, X, G, factors, lambda, iter, iter_CD, epsilon, llh, time, a_in, b_in, 
                                 C_r_in, D_r_in, D_in, C_in)
    toc()
  }
  else if(method_name == "Deriv"){
    tic("2. Estimating parameters")
    pars <- ParEst_rest_parrest_deriv(df_train, X, G, factors, lambda, iter, iter_CD, epsilon, llh, time, a_in, b_in, 
                                 C_r_in, D_r_in, D_in, C_in)
    toc()
  }
  
  # Getting predictions
  tic("3. Getting predictions")
  #Recode parameters
  newpars <- recode_results(pars, method_name, prelim$keys)
  results <- getPredict_thes(df_test[ ,c("USERID_ind_new", "OFFERID_ind_new", "CLICK",
                                    "ratioU", "ratioO", "prediction")], 
                        newpars$alpha, newpars$beta, newpars$C_pr, newpars$D_pr)
  toc()
  
  # Performance metrics
  # What to do with the NANs (some majority rule)
  results$prediction[is.na(results$prediction)] <- 0
  
  # RMSE
  RMSE <- sqrt(mean((results$prediction - results$CLICK)^2))
  
  #Calculate MPR
  results <- results %>% group_by(USERID_ind) %>% mutate("perc_rank" = rank(1-prediction)/n())
  MPR <- sum(results$CLICK*results$perc_rank)/sum(results$CLICK)
  
  # Calculate confusion matrix
  threshold <- 0.02191325 # average click rate
  results$predictionBin <- rep(0, length(results$prediction))
  results$predictionBin[results$prediction > threshold] <- 1
  
  # True positives:
  TP <- sum(results$predictionBin == 1 & results$CLICK == 1)
  TN <- sum(results$predictionBin == 0 & results$CLICK == 0)
  FP <- sum(results$predictionBin == 1 & results$CLICK == 0)
  FN <- sum(results$predictionBin == 0 & results$CLICK == 1)
  
  confusion <- list("TP" = TP, "TN" = TN, "FP" = FP, "FN" = FN)
  
  # Output
  output <- list("parameters" = pars, "prediction" = results, "MPR" = MPR, "RMSE" = RMSE,
                 "confusion" = confusion)
  return(output)
}

baselinePred_mpr <- function(df_train, df_test){
  
  # initialize column with majority
  df_test$predUser <- 0
  df_test$predOffer <- 0
  df_test$predOverall <- mean(df_train$CLICK)
  df_test$predMajority <- 0
  
  # Fill in predictions where available
  df_test$predUser <- df_test$ratioU[!is.na(df_test$ratioU)]
  df_test$predOffer <- df_test$ratioO[!is.na(df_test$ratioO)]
  
  df_test <- df_test %>% group_by(USERID_ind_new) %>% mutate("perc_rank_predUser" = rank(1-predUser)/n(), 
                                                             "perc_rank_predOffer" = rank(1-predOffer)/n(),
                                                             "perc_rank_predOverall" = rank(1-predOverall)/n(),
                                                             "perc_rank_predMajority" = rank(1-predMajority)/n())
  
  MPRUser <- sum(df_test$CLICK*df_test$perc_rank_predUser)/sum(df_test$CLICK)
  MPROffer <- sum(df_test$CLICK*df_test$perc_rank_predOffer)/sum(df_test$CLICK)
  MPROverall <- sum(df_test$CLICK*df_test$perc_rank_predOverall)/sum(df_test$CLICK)
  MPRMajority <- sum(df_test$CLICK*df_test$perc_rank_predMajority)/sum(df_test$CLICK)
  
  output <- list("df_test" = df_test, "MPRUser" = MPRUser, "MPROffer" = MPROffer, 
                 "MPROverall" = MPROverall, "MPRMajority" = MPRMajority)
  
  return(output)
}

crossValidate_thes <- function(df, X, G, FACTORS, LAMBDA, INITTYPE, ONLYVAR, folds, iter, iter_CD, method_name,
                          epsilon, warm){
  
  # Initialize a multidimensional output array
  # Rows are all the possible permutations of the huperparameters * folds
  rows <- (length(ONLYVAR) * length(FACTORS) * length(LAMBDA) * length(INITTYPE)) * folds
  
  # Columns for the hyperparameters, plus a name variable, and then all the results you want
  # these are: rmse, TP (true positive(1)), TN, FP, FN, number of iterations, best baseline, epsilon
  columns <- 11
  
  # Initialize the df (depth is the number of folds)
  CVoutput <- data.frame(matrix(NA, nrow = rows, ncol = columns))
  names(CVoutput) <- c("Epsilon", "Specification",
                       "MPR","RMSE", "TP", "TN", "FP", "FN", "Iter", "rmseUser", "DifferenceRMSE")
  
  # Now we loop
  # First we make the folds
  # Randomly shuffle the data
  df <- df[sample(nrow(df)), ]
  
  # Then assign fold indices
  foldInd <- cut(seq(1, nrow(df)), breaks = folds, labels = FALSE)
  
  row <- 1
  # Looping over your folds
  for (z in 1:folds){
    # Do onlyvar first because the train test split depends on it
    for (a in 1:length(ONLYVAR)){
      # Do onlyvar first because the train test split depends on it
      onlyVar <- ONLYVAR[a]
      
      # Make the train test split by using the foldInd and fold as input (see trainTest)
      split <- trainTest_thes(df,X, G, onlyVar, cv = TRUE, ind = foldInd, fold = z)
      df_train <-split$df_train[ ,c("USERID_ind_new", "OFFERID_ind_new", "CLICK")]
      df_test <- split$df_test
      
      # Loop the other hyperparameters
      
      for (b in 1:length(FACTORS)){
        # Initialize the warm start objects (can only be used within a certain factor size)
        a_in <- NULL
        b_in <- NULL
        C_in <- NULL
        D_in <- NULL
        C_r_in <- NULL
        D_r_in <- NULL
        
        for (c in 1:length(LAMBDA)){
          for (d in 1:length(INITTYPE)){
            tic(paste("Run", row, "out of", rows, "in fold", z, "out of ", folds))
            factors <- FACTORS[b]
            lambda <- LAMBDA[c]
            initType <- INITTYPE[d]
            
            # Run the algorithm
            if(method_name == "Eigen"){
              output <- try(fullAlg_rest(df_train, X, G, df_test, factors, lambda, iter, iter_CD, initType, method_name, 
                                         epsilon = epsilon, a_in=a_in, b_in=b_in, C_in=C_in, D_in=D_in, C_r_in = C_r_in,
                                         D_r_in = D_r_in))
              if (class(output) == 'try-error') {
                cat('Eigen method failed')
                output <- list()
                output$MPR <- NA
                output$RMSE <- NA
                output$TP <- NA
                output$TN <- NA
                output$FP <- NA
                output$FN <- NA
                output$parameters$run <- NA
              }
            }
            else{
              output <- fullAlg_rest(df_train, X, G, df_test, factors, lambda, iter, iter_CD, initType, method_name, 
                                     epsilon = epsilon, a_in=a_in, b_in=b_in, C_in=C_in, D_in=D_in, C_r_in = C_r_in,
                                     D_r_in = D_r_in)
            }
            
            
            
            # Fill the array with output
            CVoutput$Epsilon[row] <- epsilon
            
            # The name
            CVoutput$Specification[row] <- paste("Method name = ", method_name, ", Factor = ", factors, ", Lambda = ", lambda, 
                                                 ", initType = ", initType, ", onlyVar = ", onlyVar,
                                                 ", warm = ", warm, sep = "")
            
            # Performance variables
            CVoutput$MPR[row] <- output$MPR
            CVoutput$RMSE[row] <- output$RMSE
            CVoutput$TP[row] <- output$confusion$TP
            CVoutput$TN[row] <- output$confusion$TN
            CVoutput$FP[row] <- output$confusion$FP
            CVoutput$FN[row] <- output$confusion$FN
            CVoutput$Iter[row] <- output$parameters$run
            CVoutput$rmseUser[row] <- baselinePred(df_train, df_test)$rmseUser
            CVoutput$DifferenceRMSE[row] <- CVoutput$RMSE[row]-CVoutput$rmseUser[row]
            
            row <- row+1
            
            #If the output was an error, skip paramamter updates (as no updates available)
            if (class(output) == 'try-error') {
              next
              }
            
            # In case of warm starts, keep track of the last parameters
            if (warm){
              a_in <- output$parameters$alpha
              b_in <- output$parameters$beta
              C_r_in <- output$parameters$C_r
              D_r_in <- output$parameters$D_r
              if (method_name != "Old"){
                C_in <- output$parameters$C
                D_in <- output$parameters$D
              }
            }
            
            toc()
            
            
          }
        }
      }
    }
  }
  
  # Create a mean table
  CVmean <- CVoutput %>% 
    group_by(Epsilon, Specification) %>%
    summarise_all(mean)
  
  return(list("CVoutput" = CVoutput, "CVmean" = CVmean))
}



#1 Seminar --------------------------------------------------------------------------------
# 1.1 Data preparation ---------
# This is how you should import the data.
# The sequence here is important. We want to have a continuous sequence, starting at 1
# for the indices for user and order in our training set.

# Import train and game (test) set from whereever you store them
df_train <- read_delim("Observations_Report.csv",
                       ";", escape_double = FALSE, trim_ws = TRUE)
df_test <- read_delim("Observations_Game.csv",
                      ";", escape_double = FALSE, trim_ws = TRUE)

# df_train <- read_delim("~/Google Drive/Seminar 2020/Data/Observations_Report.csv",
#                        ";", escape_double = FALSE, trim_ws = TRUE)
# df_test <- read_delim("~/Google Drive/Seminar 2020/Data/Observations_Game.csv",
#                       ";", escape_double = FALSE, trim_ws = TRUE)


# Merge to create indices
df_test$CLICK <- NA
df <- rbind(df_train, df_test)

#Combine Mail and Offer id's
df[,"MailOffer"] <- df %>%
  unite("MailOffer", c("MAILID" ,"OFFERID"), sep = "_") %>%
  select("MailOffer")

# Order on click first such that NA are at bottom (no missings indices in training data)
df <- df[order(df$CLICK), ]
df <- df %>% 
  mutate(USERID_ind = group_indices(., factor(USERID, levels = unique(USERID))))
df <- df %>% 
  mutate(OFFERID_ind = group_indices(., factor(MailOffer, levels = unique(MailOffer))))

# Make it neat
df <- df[order(df$USERID_ind), c("USERID_ind", "OFFERID_ind", "CLICK")]

# Create ratios of CLICK per offer or user (== 1 or == 0 indicates no variation)
df <- df %>%
  group_by(USERID_ind) %>%
  mutate(ratioU = mean(CLICK, na.rm = TRUE)) %>%
  ungroup()

df <- df %>%
  group_by(OFFERID_ind) %>%
  mutate(ratioO = mean(CLICK, na.rm = TRUE)) %>%
  ungroup()

# Split
df_test <- df[is.na(df$CLICK), ]
df_train <- df[!(is.na(df$CLICK)), ]

# Save. Use the df_train.RDS file in CV
saveRDS(df_train, )
saveRDS(df_test, )

# saveRDS(df_train, "~/Google Drive/Seminar 2020/Data/df_train")
# saveRDS(df_test, "~/Google Drive/Seminar 2020/Data/df_test")


# 1.2 Train/test pred. -----------------------------------------------------------------
# Makes predictions for a train/test split for the FULL training set
# Also, includes columns/rows with only 0 or 1

# Use "Preparing data" first to get the df_train object
df <- readRDS(df_train)
df <- df[ ,c("USERID_ind", "OFFERID_ind", "CLICK", "ratioU", "ratioO")]

# Setting parameters
factors <- 2
lambda <- 1
iter <- 100
initType <- 4
onlyVar <- TRUE
llh <- TRUE
rmse <- TRUE
epsilon <- 0.01

set.seed(50)
split <- trainTest(df, onlyVar)
df_train <- split$df_train[ ,c("USERID_ind_new", "OFFERID_ind_new", "CLICK")]
df_test <- split$df_test[ ,c("USERID_ind_new", "OFFERID_ind_new", "CLICK", "ratioU", "ratioO", "prediction")]
rm("split")

output <- fullAlg(df_train, df_test, factors, lambda, iter, initType, llh, rmse, 
                  epsilon)

baseline <- baselinePred(df_train, df_test)

# Visualization
hist(output$prediction$prediction)
plot(output$parameters$logllh)

# 1.3 Train/test pred. SUB --------------------------------------------------------------
# Makes predictions for a train/test split for A SUBSET of the training set
# Also, includes columns/rows with only 0 or 1

# Use "Preparing data" first to get the df_train object
df <- readRDS("df_train")
df <- df[df$USERID_ind < 10000, c("USERID_ind", "OFFERID_ind", "CLICK", "ratioU", "ratioO")]

# Setting parameters
factors <- 4
lambda <- 1
iter <- 200
initType <- 4
onlyVar <- TRUE
llh <- TRUE
rmse <- TRUE
epsilon <- 0.01

set.seed(50)
split <- trainTest(df, onlyVar)
df_train <- split$df_train[ ,c("USERID_ind_new", "OFFERID_ind_new", "CLICK")]
df_test <- split$df_test[ ,c("USERID_ind_new", "OFFERID_ind_new", "CLICK", "ratioU", "ratioO", "prediction")]
rm("split")

output <- fullAlg(df_train, df_test, factors, lambda, iter, initType, llh, rmse, epsilon)

baseline <- baselinePred(df_train2, df_test2)


# Visualization
hist(output2$prediction$prediction)
xdata <- seq(1, iter+1)
plot(xdata, output$parameters$logllh, col="blue")
plot(xdata, output$parameters$rmse_it, col="red")

# 1.4 Cross validation -----------------------------------------------------------------------
# Import train set
# Make sure the names are correct
df <- readRDS("df_train")
df <- df[df$USERID_ind < 10000, c("USERID_ind", "OFFERID_ind", "CLICK", "ratioU", "ratioO")]


# Input whichever hyperparameters you want to test
FACTORS <- c(1)
LAMBDA <- c(1, 2)
INITTYPE <- c(2)
ONLYVAR <- c(TRUE)
folds <- 2
iter <- 10
epsilon <- 0.01
warm <- FALSE

CVoutput <- crossValidate(df, FACTORS, LAMBDA, INITTYPE, ONLYVAR, folds, iter, 
                          epsilon, warm)

CVoutput_mean <- CVoutput %>% group_by(epsilon, Specification) %>% summarise_all(mean)

# Visualizing output
CVoutput$Specification <- as.factor(CVoutput$Specification)

p <- ggplot(CVoutput, aes(x=Specification, y=RMSE)) + 
  geom_boxplot()
p
  
CVoutput$RMSE

# 1.5 Final predictions ----------------------------------------------------------------------
# If you want predictions for the final set

# Import test and train set
df_train <- readRDS("df_train")
df_test <- readRDS("df_test")

#Caclulcating parameters
#Hyperparameters
factors <- 2
priorsdu <- 1
priorsdi <- 1
priorlambdau <- 1/priorsdu
priorlambdai <- 1/priorsdi

pars <- getPars(df_train[ ,c("USERID_ind", "OFFERID_ind", "CLICK")], 
                factors, priorsdu, priorsdi, priorlambdau, priorlambdai)

gameResults <- getPredict(df_test, pars$alpha, pars$beta, pars$C, pars$D)

#2 Thesis----------------------------
#2.1 Data preparation-------
#data cleaning done in data_cleaning.R

#Load the cleaned data
user_data <- readRDS("~/Dropbox/Uni/Master_Econometrie/Thesis/Code/R/user_final.RDS")
item_data <- readRDS("~/Dropbox/Uni/Master_Econometrie/Thesis/Code/R/item_final.Rds")
df <- read_delim("~/Google Drive/Data scriptie/Thijs_Observations.csv", 
                           ";", escape_double = FALSE, trim_ws = TRUE)
df[,"MailOffer"] <- df %>% unite("MailOffer", c("MAILID" ,"OFFERID"), sep = "_") %>% select("MailOffer")
df <- df %>% select(-c("MAILID", "OFFERID"))%>% select("USERID", "MailOffer", everything())

#We need to recode both the offer and user indices to correspond to sparse indexes
df <- df %>% 
  mutate(USERID_ind = group_indices(., factor(USERID, levels = unique(USERID))))
df <- df %>% 
  mutate(OFFERID_ind = group_indices(., factor(MailOffer, levels = unique(MailOffer))))

#Extract a "key" to recode the users and items in user_data and item_data
users <- df %>% group_by(USERID) %>% summarize(USERID_ind = mean(USERID_ind))
items <- df %>% group_by(MailOffer) %>% summarize(OFFERID_ind = mean(OFFERID_ind))

user_data <- left_join(user_data, users)%>%select("USERID", "USERID_ind", everything()) %>% select(-c("USERID"))
item_data <- left_join(item_data, items) %>% select("MailOffer", "OFFERID_ind", everything()) %>% select(-c("MailOffer"))
df <- select(df, c("USERID_ind", "OFFERID_ind", "CLICK"))

#save the prepared datasets
saveRDS(user_data, "user_data_prep.RDS")
saveRDS(item_data, "item_data_prep.RDS")
saveRDS(df, "df_prep.RDS")


#2.1.1 Create a subset used for testing -----
set.seed(1)
#We take a subset the observations to test with. We choose 5% of the data.
sub_df <- sample_n(df, nrow(df)*0.05)
#Take the subset of the item and user data corresponding to the items in sub_df
sub_item <- item_data[item_data$OFFERID_ind %in% unique(sub_df$OFFERID_ind),]
sub_user <- user_data[user_data$USERID_ind %in% unique(sub_df$USERID_ind),]

#Now for the subset do the same recoding as before.
sub_df <- rename(sub_df, "USERID_ind_old"="USERID_ind")
sub_df <- rename(sub_df, "OFFERID_ind_old"="OFFERID_ind")
sub_item <- rename(sub_item, "OFFERID_ind_old"="OFFERID_ind")
sub_user <- rename(sub_user, "USERID_ind_old"="USERID_ind")

sub_df <- sub_df %>% 
  mutate(USERID_ind = group_indices(., factor(USERID_ind_old, levels = unique(USERID_ind_old))))
sub_df <- sub_df %>% 
  mutate(OFFERID_ind = group_indices(., factor(OFFERID_ind_old, levels = unique(OFFERID_ind_old))))

users <- sub_df %>% group_by(USERID_ind_old) %>% summarize(USERID_ind = mean(USERID_ind))
items <- sub_df %>% group_by(OFFERID_ind_old) %>% summarize(OFFERID_ind = mean(OFFERID_ind))

sub_user <- left_join(sub_user, users)%>%select("USERID_ind_old", "USERID_ind", everything()) %>% select(-c("USERID_ind_old"))
sub_item <- left_join(sub_item, items) %>% select("OFFERID_ind_old", "OFFERID_ind",everything()) %>% select(-c("OFFERID_ind_old"))
sub_df <- select(sub_df, c("USERID_ind", "OFFERID_ind", "CLICK"))

#Save the subsets
saveRDS(sub_df, "sub_df.RDS")
saveRDS(sub_item, "sub_item.RDS")
saveRDS(sub_user, "sub_user.RDS")

#Starting variables for current coding---------
sub_df <- readRDS("~/Dropbox/Uni/Master_Econometrie/Thesis/Code/R/Data/sub_df.RDS")
sub_item <- readRDS("~/Dropbox/Uni/Master_Econometrie/Thesis/Code/R/Data/sub_item.RDS")
sub_user <- readRDS("~/Dropbox/Uni/Master_Econometrie/Thesis/Code/R/Data/sub_user.RDS")

prelim <- formatdfXG(sub_df, sub_user, sub_item)
df <- prelim$df
X <- prelim$X
G <- prelim$G

factors <- 5
lambda <- 5
iter <- 5
iter_CD <- 100
epsilon <- 0.001

sourceCpp("~/Dropbox/Uni/Master_Econometrie/Thesis/Code/R/gammaui.cpp")
#initialize all tools
test <- ParEst_rest_parrest(df, X, G, factors, lambda, iter, iter_CD,  NULL, T, T)
test <- ParEst_rest_parrest_SVD(df, X, G, factors, lambda, iter, iter_CD, NULL, T, T)
test <- ParEst_rest_parrest_Eigen(df, X, G, factors, lambda, iter, iter_CD, NULL, T, T)


#Model with full dataset
df_prep <-readRDS("~/Dropbox/Uni/Master_Econometrie/Thesis/Code/R/df_prep.RDS")
user_data_prep <- readRDS("~/Dropbox/Uni/Master_Econometrie/Thesis/Code/R/user_data_prep.RDS")
item_data_prep <- readRDS("~/Dropbox/Uni/Master_Econometrie/Thesis/Code/R/item_data_prep.RDS")

prelim_full <- formatdfXG(df_prep, user_data_prep, item_data_prep)
df_full <- prelim$df
X_full <- prelim$X
G_full <- prelim$G

factors <- 5
lambda <- 0.1
iter <- 100
iter_CD <- 100
epsilon <- 0.0001

test_full <- ParEst_rest_parrest(df_full, X_full, G_full, factors, lambda, iter, iter_CD, epsilon)
test_full <- ParEst_rest_parrest_SVD(df_full, X_full, G_full, factors, lambda, iter, iter_CD, epsilon)

#Old model
factors <- 5
lambda <- 100
iter <- 100
initType <- 2
onlyVar <- TRUE
llh <- F
rmse <- FALSE
epsilon <- 0.01

output <- fullAlg(df_prep, NULL, factors, lambda, iter, initType, llh, rmse, 
                  epsilon)

# 3.1 Data preparation ---------
# This is how the full dataset and cleaned user and item data is imported and split into an train and test split

# Import full train set and cleaned dataset
df <- read_delim("~/Google Drive/Data scriptie/Thijs_Observations.csv", 
                 ";", escape_double = FALSE, trim_ws = TRUE)
user_data <- readRDS("~/Dropbox/Uni/Master_Econometrie/Thesis/Code/R/Data/user_final.RDS")
item_data <- readRDS("~/Dropbox/Uni/Master_Econometrie/Thesis/Code/R/Data/item_final.Rds")

#Combine Mail and Offer id's
df[,"MailOffer"] <- df %>%
  unite("MailOffer", c("MAILID" ,"OFFERID"), sep = "_") %>%
  select("MailOffer")

#create new unique indices for both users and items
key_user <- data.frame("USERID" = unique(df$USERID))
key_item <- data.frame("MailOffer" = (unique(df$MailOffer)))

key_user <- key_user %>% 
  mutate(USERID_ind = group_indices(., factor(USERID, levels = unique(USERID))))
key_item <- key_item %>% 
  mutate(OFFERID_ind = group_indices(., factor(MailOffer, levels = unique(MailOffer))))

#Match new names with data
df <- left_join(df, key_user) %>% left_join(key_item) %>% select(c("USERID_ind", "OFFERID_ind", "CLICK"))
user_data <- left_join(user_data, key_user) %>% select(-c("USERID")) %>% select(c("USERID_ind", everything()))
item_data <- left_join(item_data, key_item) %>% select(-c("MailOffer")) %>% select(c("OFFERID_ind", everything()))

# Create ratios of CLICK per offer or user (== 1 or == 0 indicates no variation)
df <- df %>%
  group_by(USERID_ind) %>%
  mutate(ratioU = mean(CLICK, na.rm = TRUE)) %>%
  ungroup()

df <- df %>%
  group_by(OFFERID_ind) %>%
  mutate(ratioO = mean(CLICK, na.rm = TRUE)) %>%
  ungroup()

#Split df in train and test set
split <- trainTest_thes(df, user_data, item_data, onlyVar = FALSE)

df_train <- split$df_train[ ,c("USERID_ind_new", "OFFERID_ind_new", "CLICK")]
df_test <- split$df_test[ ,c("USERID_ind_new", "OFFERID_ind_new", "CLICK", "ratioU", "ratioO", "prediction")]
X <- split$X
G <- split$G
key_user <- split$key_user
key_item <- split$key_item


# Save. Use the df.RDS file in CV
saveRDS(df,  "~/Google Drive/Data scriptie/df.RDS")
saveRDS(df_train, "~/Google Drive/Data scriptie/df_train.RDS")
saveRDS(df_test, "~/Google Drive/Data scriptie/df_test.RDS")
saveRDS(X, "~/Google Drive/Data scriptie/X.RDS")
saveRDS(G, "~/Google Drive/Data scriptie/G.RDS")

#Also save the translation keys, to be used in the formatdfXG function
saveRDS(key_user, "~/Google Drive/Data scriptie/key_user.RDS")
saveRDS(key_item, "~/Google Drive/Data scriptie/key_item.RDS")


# 3.2 Train/test pred. -----------------------------------------------------------------
# Makes predictions for a train/test split from the data preperation section

# Use "Data preparation" first to get the df_train object
df_train <- readRDS("~/Google Drive/Data scriptie/df_train.RDS")
df_test <- readRDS("~/Google Drive/Data scriptie/df_test.RDS")
X <- readRDS("~/Google Drive/Data scriptie/X.RDS")
G <- readRDS("~/Google Drive/Data scriptie/G.RDS")

# Setting parameters
factors <- 5
lambda <- 1
iter <- 2
iter_CD <- 10000
initType <- 2
llh <- TRUE
time <- TRUE
epsilon <- 0.1
method_name <- "SVD"

output <- fullAlg_rest(df_train, X, G, df_test, factors, lambda, iter, iter_CD, initType, method_name, llh, 
                  time, epsilon= epsilon)

baseline <- baselinePred(df_train, df_test)

# Visualization
hist(output$prediction$prediction)
plot(output$parameters$logllh)

# 3.3 Train/test pred. SUB --------------------------------------------------------------
# Makes predictions for a train/test split for A SUBSET of the training set

# Use "Data preparation" first to get the df_train object
df_train <- readRDS("~/Google Drive/Data scriptie/df_train.RDS")
X <- readRDS("~/Google Drive/Data scriptie/X.RDS")
G <- readRDS("~/Google Drive/Data scriptie/G.RDS")

df_subset <- df_train[df_train$USERID_ind_new < 10000,]

# Setting parameters
factors <- 5
lambda <- 1
iter <- 1000
iter_CD <- 10000
initType <- 2
onlyVar <- TRUE
llh <- TRUE
time <- TRUE
rmse <- TRUE
epsilon <- 0.005
method_name <- "SVD"

#Create new test and training split
split <- trainTest_thes(df_subset, X, G, onlyVar = onlyVar)
df_train <- split$df_train
df_test <- split$df_test
X <- split$X
G <- split$G
df_train <- split$df_train[ ,c("USERID_ind_new", "OFFERID_ind_new", "CLICK")]
df_test <- split$df_test[ ,c("USERID_ind_new", "OFFERID_ind_new", "CLICK", "ratioU", "ratioO", "prediction")]
rm("split")

output <- fullAlg_rest(df_train, X, G, df_test, factors, lambda, iter, iter_CD, initType, method_name, llh, 
                       time, epsilon= epsilon)
baseline <- baselinePred(df_train, df_test)


# Visualization
hist(output$prediction$prediction)
xdata <- seq(1, iter+1)
plot(xdata, output$parameters$logllh, col="blue")
plot(xdata, output$parameters$rmse_it, col="red")

# 3.4 Cross validation -----------------------------------------------------------------------
# Import train set
df_train <- readRDS("~/Google Drive/Data scriptie/df_train.RDS")
X <- readRDS("~/Google Drive/Data scriptie/X.RDS")
G <- readRDS("~/Google Drive/Data scriptie/G.RDS")

#To test the method we again use a subset of observations
df_subset <- df_train[df_train$USERID_ind_new < 10000,]

split <- trainTest_thes(df_subset, X, G, onlyVar = T)
X <- split$X
G <- split$G
df_train <- split$df_train[ ,c("USERID_ind_new", "OFFERID_ind_new", "CLICK")]
df_test <- split$df_test[ ,c("USERID_ind_new", "OFFERID_ind_new", "CLICK", "ratioU", "ratioO", "prediction")]
rm("split")


# Input whichever hyperparameters you want to test
FACTORS <- c(2,5,8,10,15,20)
LAMBDA <- c(250, 100, 50, 25, 20, 15, 10, 8, 6, 5, 4, 3, 2, 1, 0.5)
INITTYPE <- c(2)
ONLYVAR <- c(TRUE)
folds <- 2
iter <- 1000
iter_CD <- 10000
method_name <- "SVD"
epsilon <- 0.001
warm <- TRUE

CVoutput <- crossValidate_thes(df_subset, X, G, FACTORS, LAMBDA, INITTYPE, ONLYVAR, folds, iter, iter_CD, method_name,
                               epsilon, warm)

CVoutput_mean <- CVoutput %>% group_by(epsilon, Specification) %>% summarise_all(mean)

# Visualizing output
CVoutput$Specification <- as.factor(CVoutput$Specification)

p <- ggplot(CVoutput, aes(x=Specification, y=RMSE)) + 
  geom_boxplot()
p

CVoutput$RMSE

# 3.5 Final predictions ----------------------------------------------------------------------
# If you want predictions for the final set

df_train <- readRDS("~/Google Drive/Data scriptie/df_train.RDS")
df_test <- readRDS("~/Google Drive/Data scriptie/df_test.RDS")

X <- readRDS("~/Google Drive/Data scriptie/X.RDS")
G <- readRDS("~/Google Drive/Data scriptie/G.RDS")

#Use optimal hyperparameters
factors <- 2
lambda <- 8
iter <- 1000
iter_CD <- 10000
initType <- 2
onlyVar <- TRUE
llh <- TRUE
time <- TRUE
rmse <- TRUE
epsilon <- 0.00001
method_name <- "SVD"

pars <- fullAlg_rest(df_train, X, G, df_test, factors, lambda, iter, iter_CD, initType, method_name, llh, 
                     time, epsilon= epsilon)

gameResults <- getPredict(df_test, pars$alpha, pars$beta, pars$C, pars$D)

#3.6 Get the PR plots -------
#Load predictions
predictions <- readRDS("~/Dropbox/Uni/Master_Econometrie/Thesis/Code/R/predictions.RDS")
predictions_old <- readRDS("~/Dropbox/Uni/Master_Econometrie/Thesis/Code/R/predictions_old.RDS")
predictions_baseline <- readRDS("~/Dropbox/Uni/Master_Econometrie/Thesis/Code/R/baseline.RDS")

#Retrieve only the predictions
predictions_svd <- predictions$prediction[,c("CLICK", "prediction")]
predictions_old <- predictions_old$prediction[,c("CLICK", "prediction")]
predictions_baseline_user <- predictions_baseline$df_test[,c("CLICK", "predUser")]
predictions_baseline_offer <- predictions_baseline$df_test[,c("CLICK", "predOffer")]

#Calculate the Precision Recall curves
pr_svd <- pr.curve(scores.class0 = predictions_svd$prediction, weights.class0 = predictions_svd$CLICK, curve=T)
pr_old <- pr.curve(scores.class0 = predictions_old$prediction, weights.class0 = predictions_old$CLICK, curve = T)
pr_user <- pr.curve(scores.class0 = predictions_baseline_user$predUser, weights.class0 = predictions_baseline_user$CLICK, curve = T)
pr_offer <- pr.curve(scores.class0 = predictions_baseline_offer$predOffer, weights.class0 = predictions_baseline_user$CLICK, curve = T)

#Draw a subset of results for plotting
random_ind_svd <- sample(nrow(pr_svd$curve),10000)
random_ind_old <- sample(nrow(pr_old$curve),10000)
random_ind_user <- sample(nrow(pr_user$curve),10000)
random_ind_offer <- sample(nrow(pr_offer$curve),50000)


#Create dataframes with results
pr_svd_sub <- cbind(pr_svd$curve[random_ind_svd,1:2],"SVD")
pr_old_sub <- cbind(pr_old$curve[random_ind_old,1:2],"Unrestricted")
pr_user_sub <- cbind(pr_user$curve[random_ind_user,1:2],"User")
pr_user_offer <- cbind(pr_offer$curve[random_ind_offer,1:2],"Offer")


#Combine all results
combined_data <- rbind(pr_svd_sub,pr_old_sub,pr_user_sub, pr_user_offer)
data_plot_pr <- data.frame("Recall" = as.numeric(combined_data[,1]), "Precision"= as.numeric(combined_data[,2]),
                                                                  "Method" = as.factor(combined_data[,3]))

#Plot the final result
ggplot(data = data_plot_pr, aes(x = Recall, y = Precision, color = Method)) + geom_line(size = 1, alpha = 0.5) +
  theme_bw()+ xlab("Recall") + ylab("Precision")+ theme(legend.position = c(0.84, 0.75))


library(PRROC)
library(PerfMeas)


