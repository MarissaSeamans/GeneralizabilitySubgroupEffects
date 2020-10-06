 #============================
# R Code to implement simulations in the Invited Commentary:
# `All generalizations are dangerous, even this one.' - Dumas
# Written by Laura Balzer for Epidemiology 2017 
#
# Programmer: Laura Balzer
# lbbalzer@hsph.harvard.edu
#
# Full reference to Lesko et al: 
# C.R. Lesko, A.L. Buchanan, D. Westreich, J.K. Edwards, M.G. Hudgens, and S.R. Cole. 
# Generalizing study results: a potential outcomes perspective. Epidemiology, 2017.
#============================

args <- (commandArgs(TRUE))
xz <- c(-0.25, -0.15, -0.05, 0, 0.05, 0.15, 0.25)
xu <- c(-0.25, -0.15, -0.05, 0, 0.05, 0.15, 0.25)
xzu <- c(-0.25, -0.15, -0.05, 0, 0.05, 0.15, 0.25)

eg=expand.grid(j = xz,
               k = xu,
               l = xzu,
               StringsAsFactors = FALSE)
ieg=1
ieg=as.numeric(Sys.getenv("SGE_TASK_ID"))
if(is.na(ieg)) ieg =1

aa=eg$j[ieg]
bb=eg$k[ieg]
cc=eg$l[ieg]

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
    
    X<- rep(NA, n)
    # if randomize treatment with equal allocation
    #		ie ensure n/2 X=1 and n/2 with X=0
    if(equal.allocation){
      Xtemp<- rbinom(sum(S)/2, 1, .5)
      Xtemp2<- ifelse(Xtemp==1, 0, 1)
      Xtemp3<- sample( c(Xtemp, Xtemp2))
    }else{
      # randomize completely as done by Lesko et al.
      Xtemp3<- rbinom(sum(S), 1, .5)
    }
    X[S==1]<- Xtemp3
    
    # generate the unmeasured factor U_Y that determines the counterfactual & observed outcome
    UY <- runif(n, 0, 1)
    
    # generate the counterfactual outcomes
    Y1<- get.Y(Z, U, X=1, UY)
    Y0<- get.Y(Z, U, X=0, UY)
    
    # the observed outcome is equal to the counterfactual outcome under the observed txt
    Y<- ifelse(X, Y1, Y0)
    
    # set observed exposure and outcome to 0 if not selected/enrolled 
    X[S==0] <- 0 
    Y[S==0] <- 0 
    data.frame(Z, U, S, X, Y1, Y0, Y)
  }
  
  # get.Y: function to generate the outcome
  get.Y <- function(Z, U, X, UY){
    as.numeric( UY < (0.1073 - 0.05*X + 0.2*Z + 0.2*U + j*X*Z + k*X*U + 0*Z*U + l*X*Z*U))
  }
  
  #============================
  # Calculate the true value of the population effect with Monte Carlo simulations
  out<- get.data(n= 750000)
  PATE<- mean(out$Y1 - out$Y0)
  
  
  #============================
  # For varying sample sizes, repeatly generate the full data, calculate the SATE,
  #	and for the biased sample, implement the unadjusted, Gcomp & IPW estimator 
  nReps<- 5000
  SATE.biased <- Gcomp <- Gcomp.correct <- IPW <- IPW.correct <- unadj <- rep(NA, nReps)
  for(r in 1:nReps){
    
    #--------------------------------
    # Under biased sampling - threat to external validity
    #----------------------------------
    # draw the population of 50000 units
    full.data <- get.data(n=50000, biased.sample=T)
    
    # calculate the SATE for the enrolled units (S==1)
    # Table 1 final row
    SATE.biased[r]<- mean(full.data[full.data$S==1, 'Y1'] - full.data[full.data$S==1, 'Y0'])
    
    # subset on the observed data
    obs.data<- subset(full.data, select=c(Z,U,S,X,Y))
    
    #-----------------------------------------
    # unadjusted estimator (ave difference in mean outcomes)
    unadj[r]<- mean(obs.data[obs.data$X==1 & obs.data$S==1, 'Y']) - 
      mean(obs.data[obs.data$X==0 & obs.data$S==1, 'Y']) 	
    
    #-----------------------------------------
    # Gcomputation with fully saturated parametric regression for the outcome 
    outcome.regression <- glm(Y ~ X*Z, family='binomial', data=obs.data[obs.data$S==1, ])
    X.11<- X.10<- obs.data
    X.11$X<- 1; X.11$S<- 1
    X.10$X<- 0; X.10$S<- 1
    predict.outcome.exp<- predict(outcome.regression, newdata=X.11, type='response')
    predict.outcome.unexp<- predict(outcome.regression, newdata=X.10, type='response')
    Gcomp [r]<- mean(predict.outcome.exp - predict.outcome.unexp)
    
    outcome.regression <- glm(Y ~ X*Z*U, family='binomial', data=obs.data[obs.data$S==1, ])
    X.11<- X.10<- obs.data
    X.11$X<- 1; X.11$S<- 1
    X.10$X<- 0; X.10$S<- 1
    predict.outcome.exp<- predict(outcome.regression, newdata=X.11, type='response')
    predict.outcome.unexp<- predict(outcome.regression, newdata=X.10, type='response')
    Gcomp.correct [r]<- mean(predict.outcome.exp - predict.outcome.unexp)
    
    #-----------------------------------------
    # IPW with fully saturated parametric regression for the enrollment (selection) mechanism	
    # and with fully saturated parametric regression for the exposure mechanism
    # (alternative we know that in a randomized trial P(X=1|S=1, W) = 0.5)
    selection.regression<- glm(S~ Z*U, family='binomial', data=obs.data)
    predict.prob.select <- predict(selection.regression, type='response')	
    # exp.regression<- glm(X~ Z*U, family='binomial', data=obs.data[obs.data$S==1,])
    # predict.prob.exp <- predict(exp.regression, newdata=X.11, type='response')
    
    # calculate the denominator for the weights as the product of the probability of being selected
    #	and conditional probability of the obs. exposure
    den.11<- predict.prob.select*.5#*predict.prob.exp
    den.10<- predict.prob.select*.5#*(1-predict.prob.exp)
    IPW.correct[r] <- mean( 
      ( as.numeric(obs.data$X==1 & obs.data$S==1)/den.11 - 
          as.numeric(obs.data$X==0 & obs.data$S==1)/den.10 )*obs.data$Y )	
    
    selection.regression<- glm(S~ Z, family='binomial', data=obs.data)
    predict.prob.select <- predict(selection.regression, type='response')	
    # exp.regression<- glm(X~ Z*U, family='binomial', data=obs.data[obs.data$S==1,])
    # predict.prob.exp <- predict(exp.regression, newdata=X.11, type='response')
    
    # calculate the denominator for the weights as the product of the probability of being selected
    #	and conditional probability of the obs. exposure
    den.11<- predict.prob.select*.5#*predict.prob.exp
    den.10<- predict.prob.select*.5#*(1-predict.prob.exp)
    IPW[r] <- mean( 
      ( as.numeric(obs.data$X==1 & obs.data$S==1)/den.11 - 
          as.numeric(obs.data$X==0 & obs.data$S==1)/den.10 )*obs.data$Y )	
    
  }
  jtext <- gsub("\\.","pt",aa)
  ktext <- gsub("\\.","pt",bb)
  ltext <- gsub("\\.","pt",cc)
  # save the results
  save(PATE, SATE.biased, unadj, Gcomp, Gcomp.correct, IPW, IPW.correct, file=paste0(paste("Three_way", jtext, ktext, ltext, sep="_"),".Rdata"))
}
#j <- as.numeric(args[1])
jtext <- gsub("\\.","pt",aa)
ktext <- gsub("\\.","pt",bb)
ltext <- gsub("\\.","pt",cc)
results<- sim(j=aa, k=bb, l=cc)
#save(results, file=paste0(paste("Two_way", ktext, jtext, sep="_"), ".Rdata"))




