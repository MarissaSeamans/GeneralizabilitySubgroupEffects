#============================
# R Code to implement simulations in manuscript "Generalizability of Subgroup Effects"
# By Seamans MJ, Hong H, Ackerman B, Schmid I, Stuart EA for Epidemiology 2020
# 
# Simulation code adapted from Balzer LB. `All generalizations are dangerous, even this one.' - Alexandre Dumas. Epidemiology, 2017.
# https://github.com/LauraBalzer/On-Generalizability
# 
# Simulation setup adapted from Lesko CR, Buchanan AL, Westreich D, Edwards JK, Hudgens MG, Cole SR.
# Generalizing study results: a potential outcomes perspective. Epidemiology, 2017.
#
# Contact: Marissa Seamans
# mseamans@ph.ucla.edu
# Department of Epidemiology, UCLA
# #============================

# Create grid of values for parameters of two- and three-way interactions

args <- (commandArgs(TRUE))
az <- c(-0.15, -0.05, 0, 0.05, 0.15)
au <- c(-0.15, -0.05, 0, 0.05, 0.15)
azu <- c(-0.15, -0.05, 0, 0.05, 0.15)

eg=expand.grid(j = az,
               k = au,
               l = azu,
               StringsAsFactors = FALSE)
ieg=1
ieg=as.numeric(Sys.getenv("SGE_TASK_ID"))
if(is.na(ieg)) ieg =1

aa=eg$j[ieg]
bb=eg$k[ieg]
cc=eg$l[ieg]

sim<- function(j, k, l){
  set.seed(1)
  
#Code below adapted from Balzer 2017 https://github.com/LauraBalzer/On-Generalizability
  #============================
  # get.data: function to simulate the full data 
  #		(covariates, enrollment, txt, counterfactual outcomes, observed outcome)
  #	input: number of units, 
  #		biased.sample (if true, generate according to biased scheme in Lesko et al.; 
  #			if false, take a simple random sample)
  # output: full data 
  #============================
  get.data<- function(n, biased.sample=T){
    
    # generate baseline covariates Z and U
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
    }
    else{
    # otherwise select all units
      S<- rep(1,n)
    }
    
    # randomize treatment A completely as done by Lesko et al.
    A<- rep(NA, n)
    A[S==1]<- rbinom(sum(S), 1, .5)
    
    # generate the unmeasured factor U_Y that determines the counterfactual & observed outcome
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
    as.numeric( UY < (0.1073 - 0.05*A + 0.2*Z + 0.2*U + j*A*Z + k*A*U + 0*Z*U + l*A*Z*U))
  }
  
  #============================
  # Calculate the true value of the population effect with Monte Carlo simulations
  out<- get.data(n= 750000)
  PATE<- mean(out$Y1 - out$Y0)
  
  
  #============================
  # For nReps simulations, generate the full data, calculate the SATE,
  #	and for the biased sample, implement the unadjusted and Gcomp estimator 
  nReps<- 5000
  SATE.biased <- Gcomp <- Gcomp.correct <- unadj <- rep(NA, nReps)
  for(r in 1:nReps){
    
    #--------------------------------
    # Under biased sampling - threat to external validity
    #----------------------------------
    # draw the population of 50000 units
    full.data <- get.data(n=50000, biased.sample=T)
    
    # calculate the SATE for the enrolled units (S==1)
    SATE.biased[r]<- mean(full.data[full.data$S==1, 'Y1'] - full.data[full.data$S==1, 'Y0'])
    
    # subset on the observed data
    obs.data<- subset(full.data, select=c(Z,U,S,A,Y))
    
    #-----------------------------------------
    # unadjusted estimator (ave difference in mean outcomes)
    unadj[r]<- mean(obs.data[obs.data$A==1 & obs.data$S==1, 'Y']) - 
      mean(obs.data[obs.data$A==0 & obs.data$S==1, 'Y']) 	
    
    #-----------------------------------------
    # Gcomputation with fully saturated parametric regression for the outcome - correctly specified 
    outcome.regression <- glm(Y ~ A*Z*U, family='binomial', data=obs.data[obs.data$S==1, ])
    A.11<- A.10<- obs.data
    A.11$A<- 1; A.11$S<- 1
    A.10$A<- 0; A.10$S<- 1
    predict.outcome.exp<- predict(outcome.regression, newdata=A.11, type='response')
    predict.outcome.unexp<- predict(outcome.regression, newdata=A.10, type='response')
    Gcomp.correct[r]<- mean(predict.outcome.exp - predict.outcome.unexp)
    
    # Gcomputation with parametric regression for the outcome - exclude U, AU, and AZU to reflect scenario where U is unknown
    outcome.regression <- glm(Y ~ A*Z, family='binomial', data=obs.data[obs.data$S==1, ])
    A.11<- A.10<- obs.data
    A.11$A<- 1; A.11$S<- 1
    A.10$A<- 0; A.10$S<- 1
    predict.outcome.exp<- predict(outcome.regression, newdata=A.11, type='response')
    predict.outcome.unexp<- predict(outcome.regression, newdata=A.10, type='response')
    Gcomp[r]<- mean(predict.outcome.exp - predict.outcome.unexp)

  }
  
  jtext <- gsub("\\.","pt",aa)
  ktext <- gsub("\\.","pt",bb)
  ltext <- gsub("\\.","pt",cc)
  # save the results
  save(PATE, SATE.biased, unadj, Gcomp.correct, Gcomp, file=paste0(paste("Three_way", jtext, ktext, ltext, sep="_"),".Rdata"))
}

jtext <- gsub("\\.","pt",aa)
ktext <- gsub("\\.","pt",bb)
ltext <- gsub("\\.","pt",cc)
results<- sim(j=aa, k=bb, l=cc)