# Disclaimer: I was provided a sample code from Dr. Steingrimsson.
# FOR the Copenhagen Stroke Study 
set.seed(437)

library(pec) # Copenhagen Stroke Data
library(survival) # Surv, survreg
library(randomForest) # randomForest
library(randomForestSRC) # rfsrc


#######################################################################################################
# helper functions for data imputation (directly from Dr. Steingrimsson)
#######################################################################################################
# Calculating m1, with conditional survival function estimated using RSF procedure (p.378 from the article)
random.forest = function(obs,delta,xx) # mfunc just calls this
{
  nu = length(obs)
  
  
  # Fitting the Cox model
  # Creating the data frame
  data.used <- data.frame(obs, delta, xx)
  rand.for = rfsrc(Surv(obs,delta)~ ., data = data.used, importance = "none", ntime = NULL)
  
  # Getting the Survival Curves.
  pred.rf <- predict(rand.for, proximity = FALSE)
  
  # Calculation of the survivor function
  m1=matrix(0,nu, nu)
  
  # Finding unique event times
  time.used <- pred.rf$time.interest
  
  # Finding the jumps in the estimator \hat P(T >t|W_i)
  surv.diff <- matrix(0, ncol = sum(delta), nrow = nu) # sum(delta) = 702
  
  for(i in 1:nu){
    # Calculating the jumps in the random forest model survival curve estimator
    surv.diff[i, ] <- c(1, pred.rf$survival[i, ][-length(pred.rf$survival[i, ])]) - pred.rf$survival[i, ]
  }
  
  for(j in 1:nu)
  {
    if(delta[j]==FALSE){
      
      for(i in 1:nu)
      {
        if(obs[j]<=obs[i]){
          
          if(sum(surv.diff[i, ][time.used > obs[j]]) != 0){
            # Calculating the conditional expectation
            m1[j,i]=  sum(log(time.used[time.used > obs[j]]) * surv.diff[i, ][time.used > obs[j]])/sum(surv.diff[i, ][time.used > obs[j]])
          }
        }
      }
      if (sum(surv.diff[i, ][time.used > obs[j]]) == 0){
        m1[j,]=log(obs[j])}
    }
  }
  
  return(m1)
}

# Calculates The estimator for the censoring distribution using a KM estimator
# Input: obs = observed time, delta = failure indicator, dtype = truncation level and method.
# Output: List of bar G(tilde T|W), failure indicator changed to account for truncation
# and observed time also changed to account for truncation. 
GfuncKM=function(obs,delta,dtype)
{ 
  n <- length(obs)
  # Changing the dataset to account for truncation. 
  aa=datach(obs,delta,dtype)
  # Observed time after truncation
  obs=aa[[1]]
  # Failure indicator after truncation.
  delta=aa[[2]]
  
  #Calculating the hazard estimator. 
  hazC=mapply(function(xx,dd){dd/sum((obs>=xx))},xx=obs,dd=1-delta)
  surC_km=mapply(function(xx){prod(1-hazC[obs<=xx])},xx=obs)
  return(list(surC_km,obs,delta))	
}

# Changes the dataset to account for truncation
# Input: obs = observed time, delta = failure indicator, dtype = how truncation is done.
# b is method 2 and a is method 1 and the number afterwards indicates what the truncation level is.
# Output: obs = observed failure time adjusted for truncation, delta = failure indicator adjusted for truncaiton.
datach=function(obs,delta,dtype) # method: refer to 2016 paper
{
  nu = length(obs)
  if(dtype=="b0")
  {
    delta[obs==max(obs)]=TRUE
  }
  
  
  if(dtype=="b5")
  {
    delta[order(obs)][floor(nu-0.05*nu):nu]=TRUE
  }
  
  
  if(dtype=="b10")
  {
    delta[order(obs)][floor(nu-0.10*nu):nu]=TRUE
  }
  
  
  if(dtype=="b15")
  {
    delta[order(obs)][floor(nu-0.15*nu):nu]=TRUE
  }
  
  if(dtype=="a5")
  {    
    delta[order(obs)][floor(nu-0.05*nu):nu]=T
    obs[order(obs)][floor(nu-0.05*nu):nu]=obs[order(obs)][floor(nu-0.05*nu)]    
  }	
  if(dtype=="a10")
  {
    delta[order(obs)][floor(nu-0.10*nu):nu]=T
    obs[order(obs)][floor(nu-0.10*nu):nu]=obs[order(obs)][floor(nu-0.10*nu)]    
  }	
  
  
  if(dtype=="a15")
  {
    delta[order(obs)][floor(nu-0.15*nu):nu]=T
    obs[order(obs)][floor(nu-0.15*nu):nu]=obs[order(obs)][floor(nu-0.15*nu)]    
  }		
  
  if(dtype=="none")
  {
    obs=obs
    delta=delta
  }		
  return(list(obs,delta))	
}	

# Calculating the Tree Based conditional expectation
# Input: obs = observed time, delta = failure indicator, x1, ldots, x25 = vector of covariates
# xx = matrix of covariates, mtype = what type of cond exp is being calculated. 
# Output: Matrix of conditional expectations m_{1i}(tilde T_j) for j, i = 1, ldots ,n
mfunc=function(obs,delta,xx, mtype)
{
  
  if (mtype=="rand.for")
  {
    m1=random.forest(obs,delta,xx)
  }
  
  return(m1)
}

# external calculates the parameter vector using the L_2 loss and no simulation from martingale.
# Inputs: obs = observed time, delta = failure indicator, x1-x5 covariates,
# mtype = type of conditional expectation, dtype = truncation method and level.
externalRegularKM=function(obs,delta,xx, mtype,dtype) # for DoublyRobustL2
{
  n = length(obs)
  nu = n
  # number of observations
  
  # Calculating the conditional expectation    
  m1 = mfunc(obs,delta,xx, mtype)
  
  # Calculating the conditional censoring distribution.
  tem = GfuncKM(obs,delta,dtype)
  surC_km = tem[[1]]
  obs = tem[[2]]
  delta = tem[[3]]
  
  # Calculating a0, a1, b0, b1, c0, c1
  a0=delta/surC_km
  a1=a0*log(obs)
  
  b0=(1-delta)/surC_km
  b1=b0*diag(m1)
  
  kk=mapply(function(tt){sum((tt<=obs))},tt=obs)
  c0=mapply(function(tt){sum(b0*(obs<=tt)/kk)},tt=obs)
  c1=mapply(function(tt,i){sum(b0*(obs<=tt)*m1[,i]/kk)},tt=obs,i=1:nu)
  
  parms = c(a0,a1,b0,b1,c0,c1,obs,delta,n)
  
  return(parms)
}

# external calculates the parameter vector.
# Inputs: obs = observed time, delta = failure indicator, x1-x25 covariates,
# mtype = type of conditional expectation, dtype = truncation method and level.
BJexternal=function(obs,delta,xx, mtype,dtype, time.point) # for BuckleyJamesL2
{
  n = length(obs)
  nu = n
  # Creating the new T(t) dataset
  
  # Calculating the conditional expectation    
  m1 = mfunc(obs,delta,xx, mtype)
  
  a1 = delta *  log(obs) + (1 - delta) * diag(m1)
  
  parms = a1
  
  return(parms)
}


#######################################################################################################
# models (with modification)
#######################################################################################################

default.Wei <- function(data.train, data.test) {#
  model <- survreg(Surv(time, status) ~ . , data = data.train, dist = "wei")
  
  predict.test <- predict(model, newdata = data.test, type = "lp")
  
  # survival function from STAT 437 lecture 18 slide 
  # k shape, p scale
  S_t <- function(t, survreg.scale, survreg.lp){
    shape <- 1/survreg.scale 
    scale <- exp(survreg.lp)
    return(exp(-(t/scale)**shape))
  }
  
  surv.prob = matrix(1, nrow = nrow(data.test), ncol = 1)
  
  for(i in 1:nrow(data.test)){
    # Calculating P(T > tau1|W_i) with tau1 = 3
    surv.prob[i] <- S_t(3 * 365, model$scale, predict.test[i])
  }
  
  return(surv.prob)
}


default.RF <- function(data.train, n.tree, tau1, data.test){
  obs.train = data.train[, 1]
  delta.train = data.train[, 2]
  xx.train = data.train[, -c(1,2)]
  p = dim(xx.train)[2]
  n <- length(obs.train)
  nu <- n
  num = n
  names(data.train)[1:2] = c("time", "event")
  
  # Fitting the random survival forest using the log-rank splitting rule. 
  tree.default <- rfsrc(Surv(time, event)~ ., data = data.train, ntree = n.tree, mtry = ceiling(sqrt(p)), splitrule = "logrank")
  # Predicting from the RSF
  pred.surv.forest <- predict(tree.default, newdata = data.test)
  
  # Finding which column of pred.surv.forest$survival corresponds to
  # P(T > tau1|W)
  used.time <- sum(pred.surv.forest[['time.interest']] <= tau1) + 1
  surv.est.default <- rep(0, nrow(data.test))
  for(i in 1:nrow(data.test)){
    # Calculating P(T > tau1|W_i)
    surv.est.default[i] <- c(1, pred.surv.forest[['survival']][i, ])[used.time]
  }
  
  return(surv.est.default)
}


DoublyRobustL2 <- function(data.train, n.tree, tau1, data.test, mtype, dtype){
  obs.train = data.train[, 1] # time
  delta.train = data.train[, 2] # status
  xx.train = data.train[, -c(1,2)]
  p = dim(xx.train)[2]
  n <- length(obs.train)
  num = n
  
  # surv.est[i,j] = P(T >tau1|W_i, tree j)
  surv.est.sim <- matrix(1, nrow = nrow(data.test), ncol = n.tree)
  
  parms = externalRegularKM(obs=obs.train,delta=delta.train,xx.train, mtype=mtype, dtype=dtype)
  
  a1 <- parms[(num+1):(2*num)]
  b1 <- parms[(3*num+1):(4*num)]
  c1 <- parms[(5*num+1):(6*num)]
  y.imp.all <- a1 + b1 - c1
  
  # This for loop fits one tree in the forest at each iteration
  for(j in 1:n.tree){
    
    # bootstrap
    bs <- sample(1:num, size = num, replace = TRUE)
    
    # Creating the imputed dataset
    xx <- xx.train
    
    # Fitting a single tree sqrt(p) vs p/3
    # node size needs to be tuned
    tree.one = randomForest(xx[bs, ], y.imp.all[bs], ntree=1, mtry=ceiling(sqrt(p))) # 269 nodes for j = 1
    

    # Predicting the value of the observations in the test set in order to know in what
    # terminal node they fall. 
    predict.test <- predict(tree.one, newdata = data.test)
    predict.train <- predict(tree.one, newdata = xx)
    
    # Computing the Kaplan-Meier estimator for the simulated martinagel imputed tree
    for(i in 1:length(unique(predict.train))){
      # Finding observations in the bootstrap sample that fall in that node
      obs.term.node <- which(predict.train == unique(predict.train)[i])
      
      # Finding the observations falling in that terminal node
      delt <- delta.train[obs.term.node]
      obst <- obs.train[obs.term.node]
      
      # Finding the observations in the test set that fall in that node
      obs.orig <- which(predict.test == unique(predict.train)[i])
      # Calculating the probability P(T > t|W_test) using the Kaplan Meier estimator
      #for that bootstrap tree.
      
      # Finding P(T >tau1|W) for that terminal node using the Kaplan Meier product limit estimator on the dataset
      hazT=mapply(function(xx,dd){dd/sum(obst>=xx)},dd=delt,xx=obst)
      surv.est.sim[obs.orig, j] = prod(1-hazT[obst <= tau1])
    } #end i loop
  } # End j loop
  
  # Predicting P(T > tau1|W) by averaging over the KM estimator predictions
  pred.surv.sim <- apply(surv.est.sim, 1, mean)
  
  
  # Calculating the absolute survival difference.
  error.new.imp <-  pred.surv.sim
  return(error.new.imp)
}


BuckleyJamesL2 <- function(data.train, n.tree, tau1, data.test, mtype, dtype){
  
  obs.train = data.train[, 1]
  delta.train = data.train[, 2]
  xx.train = data.train[, -c(1,2)]
  p = dim(xx.train)[2]
  n <- length(obs.train)
  nu <- n
  num = n
  
  # surv.est[i,j] = P(T >tau1|W_i, tree j)
  surv.est.sim <- matrix(1,nrow = nrow(data.test), ncol = n.tree)
  
  parms = BJexternal(obs=obs.train,delta=delta.train,xx.train, mtype=mtype, dtype=dtype)
  
  y.imp.all <- parms
  
  # This for loop fits one tree in the forest at each iteration
  for(j in 1:n.tree){
    
    # bootstrap
    bs <- sample(1:num, size = num, replace = TRUE)
    
    # Creating the imputed dataset
    xx <- xx.train
    
    # Fitting a single tree
    tree.one = randomForest(xx[bs, ], y.imp.all[bs], ntree=1, mtry=ceiling(sqrt(p)))

    # Predicting the value of the observations in the test set in order to know in what
    # terminal node they fall. 
    predict.test <- predict(tree.one, newdata = data.test)
    predict.train <- predict(tree.one, newdata = xx.train)
    
    # Computing the Kaplan-Meier estimator for the simulated martinagel imputed tree
    for(i in 1:length(unique(predict.train))){
      # Finding observations in the bootstrap sample that fall in that node
      obs.term.node <- which(predict.train == unique(predict.train)[i])
      
      # Finding the observations falling in that terminal node
      delt <- delta.train[obs.term.node]
      obst <- obs.train[obs.term.node]
      
      # Finding the observations in the test set that fall in that node
      obs.orig <- which(predict.test == unique(predict.train)[i])
      # Calculating the probability P(T > t|W_test) using the Kaplan Meier estimator
      #for that bootstrap tree.
      
      # Finding P(T >tau1|W) for that terminal node using the Kaplan Meier product limit estimator on the dataset
      hazT=mapply(function(xx,dd){dd/sum(obst>=xx)},dd=delt,xx=obst)
      surv.est.sim[obs.orig, j] = prod(1-hazT[obst <= tau1])
    } #end i loop
  } # End j loop
  
  # Predicting P(T > tau1|W) by averaging over the KM estimator predictions
  pred.surv.sim <- apply(surv.est.sim, 1, mean)
  
  
  # Calculating the absulate survival difference.
  error.new.imp <-  pred.surv.sim
  return(error.new.imp)
}


#######################################################################################################
# main()
#######################################################################################################
# we will run 4 models: Weibull model; default random forest; doubly robust L2; and buckley-james L2.
run = function(data.used){  
  tau1 = 3 * 365
  nu <- dim(data.used)[1]
  
  obs <- data.used[, 14]
  delta = data.used[, 15]
  
  
  fitset = 1:nu  
  # Creating a test set and a training set
  # Choose level of v-fold cross validation 
  v <- 10
  # Creating the cross validated group. 
  get.sam <- rep(c(1:v),sum(delta[fitset]==TRUE)/v)
  if(length(get.sam)<sum(delta[fitset]==TRUE)){
    get.sam <- c(get.sam,c(1:(sum(delta[fitset]==TRUE)-length(get.sam))))
  }
  delta.1 <- sample(get.sam,sum(delta[fitset]==TRUE),replace=F)
  delta.0 <- sample(rep(c(1:v),nu-length(delta.1)),sum(delta[fitset]==FALSE),replace=F)
  k <- 1
  l <- 1
  grp.delt <- NULL
  for(m in 1:nu){
    if(delta[m]==0){
      grp.delt[m] <- delta.0[k]
      k <- k+1
    }
    if(delta[m]==1){
      grp.delt[m] <- delta.1[l]
      l <- l+1
    }
  }
  
  # Initializing Brier Scores for four models
  default.Wei.score = NULL
  default.rf.score = NULL
  dr.l2.score = NULL
  bj.l2.score = NULL
  
  
  for(k in 1:v){
    test.set <- grp.delt == k
    train.set <- grp.delt != k
    data.test <- data.used[test.set, ]
    data.train <- data.used[train.set, ]
    
    obs.train <- data.train[, 14]
    delta.train <- data.train[, 15]
    obs.test <- data.test[, 14]
    delta.test <- data.test[, 15]
    
    # covariates
    xx.train <- data.train[, -c(14,15)]
    xx.test <- data.test[, -c(14,15)]
    xx.all <- data.used[, -c(14, 15)]
    
    
    data.test = xx.test
    
    # Setting obs and delta to be first in data.train
    data.train = data.train[, c(14, 15, 1:13)]
    
    n <- length(obs)
    num <- length(obs)
    
    # Number of trees in the forest.
    n.tree <- 1000
    
    mtype = "rand.for"
    dtype = "b10" # truncation method
    obs.test.t = pmin(obs.test, tau1) # tau1 = 3 * 365 to predict S0(t|W) where t = 3 * 365, page 381
    delta.test.t = delta.test * (obs.test <= tau1) + (obs.test > tau1)
    a0 <- delta.test.t/GfuncKM(obs.test.t,delta.test.t,dtype)[[1]]
    
    default.w = default.Wei(data.train, data.test)
    default.Wei.score = c(default.Wei.score, mean(a0 * (default.w - (obs.test > tau1))^2))
    print(default.Wei.score)
    
    default.rf <-  default.RF(data.train, n.tree, tau1, data.test)
    default.rf.score = c(default.rf.score, mean(a0 * (default.rf - (obs.test > tau1))^2))
    print(default.rf.score)
    
    dr.l2 =DoublyRobustL2(data.train, n.tree, tau1, data.test, mtype, dtype)
    dr.l2.score = c(dr.l2.score, mean(a0 * (dr.l2 - (obs.test > tau1))^2))
    print(dr.l2.score)
    
    bj.l2 = BuckleyJamesL2(data.train, n.tree, tau1, data.test, mtype, dtype)
    bj.l2.score = c(bj.l2.score, mean(a0 * (bj.l2 - (obs.test > tau1))^2))
    print(bj.l2.score)
  }
  
  result = data.frame(default.Wei.score = mean(default.Wei.score), default.rf.score = mean(default.rf.score), dr.l2.score = mean(dr.l2.score), bj.l2.score = mean(bj.l2.score))
  return(result)
}


main <- function(){
  data(cost)
  
  # negligible variations to time in order to break ties
  cost$time <- cost$time + runif(length(cost$time), 0, 0.5)
  
  result <- run(cost)
  return(result)
}

result = main()

print(result)



