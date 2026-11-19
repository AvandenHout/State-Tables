# Simulate four-state model with and
# without missing data by design

# Ardo, UCL 2025

# Prelim:
library(msm)
set.seed(12345)
source("SimSummary.r")
digits <- 3

# Missing by design?
DESIGN   <- TRUE
censored <- -2

# Number of simulations
S <- 200

# True values:
beta   <- c(-3, -4, -2)
p.true <- c(beta)

# Simulate 4-state data:
grid  <- seq(0, 10, 1)
N     <- 1500

# Q matrix:
q12 <- exp(beta[1])
q14 <- exp(beta[2])
q23 <- exp(beta[3])
q24 <- q14
Q <- matrix( c(-(q12+q14), q12,0, q14,
               0, -(q23+q24), q23, q14,
               0,0,0,0,
               0,0,0,0), 4,4,byrow=TRUE)
cat("\nGiven Q:\n")
print(Q)

# Data frame:
L      <- length(grid)
sim.df <- data.frame(subject = rep(1:N, rep(L,N)), 
                     time = rep(grid,N))

# Run simulation:
simresults    <- matrix(NA,S,3)
simresults.SE <- matrix(NA,S,3)
simconv       <- rep(NA,S)
cat("\n---------------")
for(s in 1:S){

 # Monitor:
 cat("\ns = ",s,"\n")
  
 # Simulate data:
 dta <- simmulti.msm(sim.df, Q, start = rep(1,N), death = c(3,4))
 #print(statetable.msm(state, subject, data=dta))
 
 # Missing by Design:
 if(DESIGN){
   # Save original data
   dta.all <- dta
   # Control group:
   for(i in 1:(N/2)){
     dta.i <- dta[dta$subject==i,]
     n.i   <- nrow(dta.i)
     dta.i <- dta.i[c(1,n.i),]
     if(dta.i$state[2]%in%c(1,2)){dta.i$state[2] <- censored}
     if(i==1){dta.c <- dta.i}else{dta.c <- rbind(dta.c,dta.i)}
   }
   #cat("Statetable for controls:\n")
   #print(statetable.msm(state, subject, data=dta.c))
   # Screening group:
   for(i in (N/2+1):N){
     dta.i <- dta[dta$subject==i,]
     if(2%in%dta.i$state){
       n.i   <- min(which(2==dta.i$state))
       dta.i <- dta.i[1:n.i,]
     }
     if(i==(N/2+1)){dta.s <- dta.i}else{dta.s <- rbind(dta.s,dta.i)}
   }
   #cat("Statetable for screened:\n")
   #print(statetable.msm(state, subject, data=dta.s))
   # Combine the groups:
   dta <- rbind(dta.c,dta.s)
   cat("Statetable for combined sample:\n")
   print(statetable.msm(state, subject, data=dta))
   
 }
 
 # Fit model and check:
 if(DESIGN){
  model <- msm( state ~ time, subject=subject, data = dta,
                qconstraint = c(1,2,3,2),
                deathexact = c(3,4), qmatrix = Q, 
                censor=censored, censor.states = c(1,2))
 }else{
  model <- msm( state ~ time, subject=subject, data = dta,
                qconstraint = c(1,2,3,2),
                deathexact = c(3,4), qmatrix = Q)
 }
 cat("Convergence code = ",model$opt$convergence,"\n")
 print(round(model$opt$par,digits))
 cat("True values:\n")
 print(beta)
 cat("---------------")
 
 # Save:
 simresults[s,]    <- model$opt$par
 simresults.SE[s,] <- sqrt(diag(model$covmat[1:3,1:3]))
 simconv[s]        <- model$opt$convergence
}
                       
# Simulation results:
cat("\nS = ",S)
cat("\nN = ",N,"\n")
SimSummary(p.true=p.true, results=simresults,
           results.SE = simresults.SE, )

####################
#S =  100
#N =  500 
#Number of simulation S = 100 .
#Summary simulations results:
#true   mean  bias rel.bias  RMSE
#[1,]   -3 -3.000 0.000    0.000 0.080
#[2,]   -4 -4.000 0.000    0.000 0.107
#[3,]   -2 -1.993 0.007    0.004 0.106
# WITH DESIGN:
#true   mean   bias rel.bias  RMSE
#[1,]   -3 -3.002 -0.002    0.001 0.112
#[2,]   -4 -4.002 -0.002    0.001 0.108
#[3,]   -2 -1.991  0.009    0.005 0.206
#S =  200
#N =  1000 
#  true   mean  bias rel.bias  RMSE
#[1,]   -3 -2.999 0.001    0.000 0.078
#[2,]   -4 -3.996 0.004    0.001 0.076
#[3,]   -2 -1.997 0.003    0.001 0.154
#S =  200
#N =  1500 
#  true   mean  bias rel.bias  RMSE
#[1,]   -3 -2.999  0.001    0.000 0.061
#[2,]   -4 -3.996  0.004    0.001 0.063
#[3,]   -2 -2.009 -0.009    0.004 0.130