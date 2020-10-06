.1073 - 0.05*A + 0.1 + 0.1 *U + 0.2*A*Z + 0.2*A*U + 0*Z*U + 0.2*A*Z*U

sim<- function(j, k, l){
  set.seed(1)
  
  #============================
  # get.data: function to simulate the full data 
  #		(covariates, enrollment, txt, counterfactual outcomes, observed outcome)
  #	input: number of units, 
  #		biased.sample (if true, generate according to biased scheme in Lesko et al.; 
  #			if false, take a simple random sample)
  #		equal.allocation (if true, ensure n/2 units receive the treatment & n/2 receive the control;
  #			if false, completely randomize)
  # output: full data 
  #============================
  get.data<- function(n, biased.sample=F, equal.allocation=F){
    
    # generate baseline covariates
    Z <- rbinom(n, 1, .15)
    U <- rbinom(n, 1, .20)
    
    # if a biased sample, then selection probability depends on covariates
    if(biased.sample){
      strata1 <- sample(which(Z==0 & U==0), 320)
      strata2 <- sample(which(Z==1 & U==0), 480)
      strata3 <- sample(which(Z==0 & U==1), 480)
      strata4 <- sample(which(Z==1 & U==1), 720)
      
      S<- rep(0, n)
      S[c(strata1, strata2, strata3, strata4)]<- 1
    }else{
      # otherwise select all units
      S<- rep(1,n)
    }
    
    A<- rep(NA, n)
    # if randomize treatment with equal allocation
    #		ie ensure n/2 A=1 and n/2 with A=0
    if(equal.allocation){
      Atemp<- rbinom(sum(S)/2, 1, .5)
      Atemp2<- ifelse(Atemp==1, 0, 1)
      Atemp3<- sample( c(Atemp, Atemp2))
    }else{
      # randomize completely as done by Lesko et al.
      Atemp3<- rbinom(sum(S), 1, .5)
    }
    A[S==1]<- Atemp3
    
    # generate the unmeasured factor U_Y that determines the counterfactual & observed outcome
n<-750000
UY <- runif(n, 0, 1)
    
    # generate the counterfactual outcomes
    Y1<- get.Y(Z, U, A=1, UY)
    Y0<- get.Y(Z, U, A=0, UY)
    
    # the observed outcome is equal to the counterfactual outcome under the observed txt
    Y<- ifelse(A, Y1, Y0)
    
    # set observed exposure and outcome to 0 if not selected/enrolled 
    A[S==0] <- 0 
    Y[S==0] <- 0 
    data.frame(Z, U, S, A, Y1, Y0, Y)
  }
  
  # get.Y: function to generate the outcome
  get.Y <- function(Z, U, A, UY){
    as.numeric( UY < (0.1073 - 0.05*A + 0.2*Z + 0.2*U + 0.15*A*Z + 0.15*A*U + 0*Z*U + 0.15*A*Z*U))
    rbinom(n, 1, .15)
  }
  
  out<- get.data(n= 750000)
  PATE<- mean(out$Y1 - out$Y0)
  
  #--------------------------------
  # Under biased sampling - threat to external validity
  #----------------------------------
  # draw the population of 50000 units
  full.data <- get.data(n=50000, biased.sample=T)
  
  # calculate the SATE for the enrolled units (S==1)
  # Table 1 final row
  #SATE.biased[r]<- mean(full.data[full.data$S==1, 'Y1'] - full.data[full.data$S==1, 'Y0'])
  SATE.biased<- mean(full.data[full.data$S==1, 'Y1'] - full.data[full.data$S==1, 'Y0'])
  
  # subset on the observed data
  obs.data<- subset(full.data, select=c(Z,U,S,A,Y))
  
  #-----------------------------------------
  # unadjusted estimator (ave difference in mean outcomes)
  # unadj[r]<- mean(obs.data[obs.data$X==1 & obs.data$S==1, 'Y']) - 
    unadj<- mean(obs.data[obs.data$A==1 & obs.data$S==1, 'Y']) - 
    mean(obs.data[obs.data$A==0 & obs.data$S==1, 'Y']) 	
  
  #-----------------------------------------
  # Gcomputation with fully saturated parametric regression for the outcome 
  outcome.regression <- glm(Y ~ A*Z, family='binomial', data=obs.data[obs.data$S==1, ])
  A.11<- A.10<- obs.data
  A.11$A<- 1; A.11$S<- 1
  A.10$A<- 0; A.10$S<- 1
  predict.outcome.exp<- predict(outcome.regression, newdata=A.11, type='response')
  predict.outcome.unexp<- predict(outcome.regression, newdata=A.10, type='response')
  #Gcomp [r]<- mean(predict.outcome.exp - predict.outcome.unexp)
  Gcomp<- mean(predict.outcome.exp - predict.outcome.unexp)
  
  outcome.regression <- glm(Y ~ A*Z*U, family='binomial', data=obs.data[obs.data$S==1, ])
  A.11<- A.10<- obs.data
  A.11$A<- 1; A.11$S<- 1
  A.10$A<- 0; A.10$S<- 1
  predict.outcome.exp<- predict(outcome.regression, newdata=A.11, type='response')
  predict.outcome.unexp<- predict(outcome.regression, newdata=A.10, type='response')
  Gcomp.correct<- mean(predict.outcome.exp - predict.outcome.unexp)
  