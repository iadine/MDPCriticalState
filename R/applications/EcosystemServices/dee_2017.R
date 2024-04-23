### Dee et al. SI 2: Numerical analyses for "To what extent can ecosystem services motivate protecting biodiversity?"
## Written by Laura Dee; last updated April 20, 2017 
# Please cite as Dee et al. To what extent can ecosystem services motivate protecting biodiversity? Ecology Letters
# load these libraries
library(ggplot2)
library(viridis)
#####################################################################################################
### Code for the analytical solution (equation [3] main text) to check with numerical analyses; #####
### see SI Materials and Methods 1 and equation [3] in the main text for details ####################
#####################################################################################################
## set parameter values
v= 90
c = 80
Beta = v/c
k = 10
S_0 = 250
discount = .95
## the analytical solution [3]
sBar <-ceiling((k*((v/c) - discount))/(1-discount))
print(sBar)
##############################################################################################################
## Numerical Analyses
########################################################################################
##############################################################################################################
## Part I. Performing value function iteration to solve the problem in which all k critical species are needed to
# provide the service & checking code the numerical solution against the analytical solution from equation [3]
##############################################################################################################
########################################################################################
## Value Iteration over an infinite time horizon ######################################
########################################################################################
# remove previously stored variables
rm(list=ls())
####### STEP 1 - define objectives (no need to code)
####### STEP 2 - define states
# state space limit 
S_0 <- 250 # total/maximum number of species
# Create vector of all possible states
states <- 0:S_0
######################################################################
## STEP 3 - Defining set of control actions ###########################
#######################################################################
# protect = 1, not protect = 0
D <- c(0,1)
#####################################################################
# STEP 4: DEFINE DYNAMICS ##########################################
#####################################################################
dynamic <- function(states,action) {
  nextpop <- states - (1 - action)
  nextpop1 <- pmax(nextpop, 0)
  return(nextpop1)
}
##################################################################################
# STEP 5: DEFINE UTILITY ##########################################################
##################################################################################
## Set model parameters within utility function: value (v), costs (c), number of critical species (k), and the discount factor.
# These parameters can be changed
v = 90 #value of ES
c = 80 #costs incurred by protection
k = 10 #number of critical species
discount = 0.95 ## Discount factor
## Utility function for the "all or nothing" cases where all k critical species are needed to obtain the service
get_utility <- function(s_t, a_t, v, c, k) {
  if(a_t==1) {
    if(s_t<k) {
      return(-c)
    } else {
      return(v-c)
    }
  } else if(a_t==0) {
    if(s_t < k) {
      return(0)
    } else {
      return(((s_t - k )/s_t)*v)
    }
  }
}
##################################################################################
# STEP 6: SOLVE BELLMAN EQUATION WITH VALUE FUNCTION ITERATION! ##################
##################################################################################
#Utility matrix
utility <- array(0, dim = c(length(states), length(D)))
#transition matrix (of states)
trans <- array(0,dim=c(length(states), length(D)))
#transition probability matrix 
p_kc <- numeric(length=length(states)) 
### Fill in the transition and utility matrix using a loop across all states and all actions 
##Loop over all states
for (s_t in 0:S_0) { 
  
  #Compute probability of keeping all critical species
  p_kc[s_t+1] <- max((s_t-k)/s_t,0)
  
  # Loop on all actions - here decision is just 0 or 1 (not protect, or protect)
  for (i in 0:1) {
    # Calculate the transition state at the next step, given the current state s_t and the decision to protect or not i 
    s_next <- dynamic(s_t,i)
    
    # Find next state index
    index <- which(states == s_next)
    # store the transition matrix of the possible states (sizes of s) 
    trans[s_t+1,i+1] <- index
    
    # Compute utility at each s for each action 
    # s_t is actual number of species, i is the actual action. We must store them as +1 due to indexing
    utility[s_t+1,i+1] <- get_utility(s_t, i, v, c, k)
    print(utility)
  } # end of action loop!
} # end of state loop!
## Action-value vector at t and t+1
Vt <- numeric(length(states))
Vtplus <- numeric(length(states))
# Optimal Policy vector
Optpol <- numeric(length(states))
# We define a factor and convergence cirterion to check for convergence time. Criterion used from Marescot et al (2013) 
did_converge <- FALSE 
discount <- 0.95
epsilon <- 0.0001
# infinite horizon - value iteration
while(did_converge==FALSE) { 
  # We define a matrix Q that stores the updated action values for #all states (rows) # actions (columns)!
  Q <- array(0, dim=c(length(states), length(D)))
  
  #Compute dynamic value of taking each action given current guess for value function
  #1. Compute expected value of not protecting. Note utility[,1] accounts for the effect of the probability of losing
  # a service providing species on the current period value, but we need to account for the effect on the 
  # future value here
  Q[,1] <- utility[,1] + discount*p_kc*Vtplus[trans[,1]]
  
  #2. compute value of protecting
  Q[,2] <- utility[,2] + discount*Vtplus[trans[,2]] 
  
  # Find the optimal action value at time t is the maximum of Q! ##
  Vt <- apply(Q, 1, max)
  if(max(abs(Vtplus-Vt)) <= epsilon*(1-discount)/(2*discount)) did_converge <- TRUE 
  
  # After filling vector Vt of the action values at all states, we update the vector: Vt+1 to Vt!
  # and we go to the next step standing for previous time t-1, since we iterate backwards.
  # this is our updated guess of the value function 
  Vtplus <- Vt
} # end of while loop for value iteration over infinite horizon
##################################################################################
# Find optimal action for each state ############################################
##################################################################################
for (s_t in 0:S_0) {
  # We look for each state which column of Q corresponds to the maximum of the last 
  #updated value of Vt (the one at time t+1). If the index vector is longer than 
  # 1(if there is more than one optimal value we chose the minimum) 
  Optpol[s_t+1] <- D[(max(which(Q[s_t+1,] == Vt[s_t+1])))]
}
########################################################################################
## Get & plot results ###################################################################
#########################################################################################
## Plot the optimal policy, where 1 means to protect at s_t species and 0 means not to protect at s_t species
# *Note that the first column of Optpol is actually corresponding to s = 0 *
plot(Optpol, xlab = "s_t", main = "Optimal Policy of species protection")
lines(Optpol, col = "blue")
## plot the expected net value resulting from this optimal policy 
plot(Vt, xlab = "s_t", ylab = "Vt, expected net value over time", main = "Expected 
     net ecosystem service value over infinite time horizon")
#########################################################################################################
########################################################################################################
## PART II: More general payoff functions where now the bellman equation is a function of s and r ######
### Using value iteration over infinite time horizon ###################################################
#########################################################################################################
#To clear previous variables: 
rm(list=ls())
####### STEP 1 - define objectives (no need to code)
####### STEP 2 - run functions that define control actions, the dynamics, and the utility function
####### STEP 3 - run wrapper function "ComputeNumericalSolutions" that solves the Bellman equation with value iteration
####### STEP 4 - specify parameters as inputs to the function "ComputeNumericalSolutions" and solve 
####### STEP 5 - plots for results; with code to create all figures in main text and SI 1 
###################################################################
# Defining set of control actions #################################
###################################################################
# protect = 1, not protect = 0
D <- c(0,1)
####################################################################
# DEFINE DYNAMICS ##################################################
####################################################################
dynamic <- function(states_s,states_r,action) {
  s_t = states_s
  r = states_r 
  
  if(s_t==0) {
    #have next s=0, next r=0 with probability 1
    return(matrix(c(0, 0, 1), nrow=1))
  }
  
  if(action==0) {
    nextpop_s <- pmax(s_t - (1 - action), 0)
    probLoss <- r/s_t
    probKeep <- (s_t - r)/s_t
    
    #2 possible outcomes: r declines or it doesn't
    #if r declines:
    nextpop_r <- pmax(r - (1 - action), 0)
    outcomeRLoss = c(nextpop_s, nextpop_r, probLoss)
    
    #if r does not decline
    if(s_t>r) {
      outcomeRKeep = c(nextpop_s, r, probKeep)
      
      #return a matrix with one row for each possible outcome
      #and 3 columns per row: next s, next r, probability of that outcome
      return(rbind(outcomeRLoss, outcomeRKeep))
    } else {
      return(matrix(outcomeRLoss, nrow=1)) ## this just says that when s=r that NOT losing a species in r is not a feasible state
      # one of the r must be lost with a prob = 1 when s=r and the action is not protecting so we return only that outcome. 
    }
  } else {
    outcomeProb = 1
    return(matrix(c(s_t, r, outcomeProb), nrow=1))
  }
}
#########################################################################################################
##DEFINE UTILITY ########################################################################################
##This part includes the different relationships between current-period ES value and number of species ###
##########################################################################################################
### Instantaneous payoff/current-period utility function ### 
# Modify this utility function to model different relationships between current-period ES value and number of species. 
# The relationship is based on a power function v_r = a*r_t^b where r are the number of species present that can provide the service.
curV = function(curR, a, b) {
  return(a*curR^b) 
  
  ### To verify the code/do check against analytical solution, using a step function as in the analytical case
  # this is an example where v = 90 (when all k critical species are present, i.e., r_0 = k), as in the Part 1 code: #
  # v <- rep(0, length(curR)) # use vector rather than "if" test since need a vector of r 
  # v[curR>=r_0] <- 90 #find elements in curR that are >= r_0 and change corresponding entries of v from 0 to 90 
  # return(v) 
}
###Specify utility function (current period payoffs) ######
# Utility function (states is a vector c(s_t,r)), 
#vFun refers to the current period utility function to use (call the function curV above, after specifying which relationship/params to use). 
get_utility <- function(states, a_t, c, vFun, ...) {
  ##states
  s_t = states[1]
  r = states[2]
  
  if(a_t==1) {
    return(vFun(r, ...)-c)
  } else if(a_t==0) { 
    states_next <- dynamic(s_t, r, a_t)
    probPerOutcome = states_next[,3]
    rPerOutcome = states_next[,2]
    #compute expected utility
    expUtility = sum(probPerOutcome * vFun(rPerOutcome, ...))
    return(expUtility)
  }
}
##################################################################################
# FUNCTION TO SOLVE BELLMAN EQUATION WITH VALUE ITERATION! ########################
##################################################################################
ComputeNumericalSolutions = function(a, b, discount, c, S_0, r_0){
  # vector of all possible states
  states_s <- 0:S_0
  states_r <- 0:r_0
  
  # Utility matrix
  utility <- array(NA, dim = c(length(states_s), length(states_r), length(D)))
  # check dimensions to make sure this works
  # dim(utility)
  
  ## Fill in the utility matrix using a loop across all states and all actions 
  ####Loop on all states 
  for (s_t in 0:S_0) { 
    #Loop over feasible values of r. r can't be bigger than s_t so we loop to the min of s_t and r.
    for(r in 0:min(s_t,r_0)) {
      # Loop on all actions - here it is just 0 or 1 so i did i in 0:1 rather than i in 1:length(D)
      for (i in 0:1) { 
        # Compute utility for each s for each action and r
        # s_t is actual number of species, i is the actual action. We must store them as +1 due to indexing
        #curV is the payoff function we are using.
        utility[s_t+1, r+1, i+1] <- get_utility(c(s_t, r), i, c, curV, a, b) 
      } # end of action loop!
    }
  } # end of state loop!
  
  ## Action-value vector at t and t+1
  Vt <- array(NA, dim = c(length(states_s), length(states_r)))
  Vtplus <- array(-50000, dim = c(length(states_s), length(states_r)))
  
  ## Optimal Policy vector
  Optpol <- array(NA, dim = c(length(states_s), length(states_r)))
  
  ## We define a factor and convergence cirterion to check for convergence time following Marescot et al (2013)
  did_converge <- FALSE 
  epsilon <- 0.0001
  
  ## infinite horizon - value iteration
  while(did_converge==FALSE){ 
    #did_converge==TRUE
    
    ## We define a matrix Q that stores the updated action values for #all states (rows) # actions (columns)!
    Q <- array(0, dim=c(length(states_s), length(states_r), length(D)))
    
    #Compute dynamic value of taking each action given current guess for value function
    #1. Compute expected value of not protecting. Note utility[,1] accounts for the effect of the probability of losing
    # a service providing species on the current period value, but we need to account for the effect on the 
    # future value here
    
    ### calc expected future value from not protecting #########
    # the loop inputs a specific s an a specific r. We are calculating only the states and value that result from not protecting (i.e., a=0)
    for(s in 0:S_0) {
      for(r in 0:min(s,r_0)) {
        # print(paste(s,r,sep=",")) # check to make sure it goes through each combination of s and r
        outcomes = dynamic(s, r, action=0)
        #print(outcomes) #check to see what transitions look like
        vnp <- utility[s+1,r+1,1] # we start off with the current value which has already been computed and stored.
        # next we need to compute and store the expected future value
        # starting with the expected future value of not protecting 
        for(j in 1:nrow(outcomes)) {
          outcome_s = outcomes[j,1]
          outcome_r = outcomes[j,2]
          outcome_p = outcomes[j,3]
          vnp <- vnp + discount*(outcome_p*Vtplus[outcome_s+1, outcome_r+1]) #sequentially 
          # adding the value of not protecting's expected future value for all possible future states
        }
        Q[s+1,r+1,1] <- vnp #store the value of not protecting while indexing with s+1, r+1 and 1 = NP
      }
    }
    
    #2. compute the value of protecting (this is easier because there is only one possible outcome :) ) 
    Q[,,2] <- utility[,,2] + discount*Vtplus[,] 
    
    # Find the optimal action value at time t is the maximum of Q!
    # Apply for each given s,r pair (the c(1,2)) and find the max across actions. We provide c(1,2) because we have a 3-dimensional
    # matrix. apply it for each combination of s and r. 
    Vt <- apply(Q, c(1,2), max)
    
    #na.rm=T says to ignore NAs in computing max
    if(max(abs(Vtplus-Vt), na.rm=T) <= epsilon*(1-discount)/(2*discount)) did_converge <- TRUE 
    
    # After filling vector Vt of the action values at all states, we update the vector … Vt+1 to Vt!
    # and we go to the next step standing for previous time t-1, since we iterate backward!
    Vtplus <- Vt
    
  } # end of while loop for value iteration over infinite horizon
  
  # Find optimal action for each state
  for (s_t in 0:S_0) {
    for(r in 0:min(s_t,r_0)) {
      # We look for each state which column of Q corresponds to the maximum of the last 
      # updated value of Vt (the one at time t+1).
      
      Optpol[s_t+1, r+1] <- D[(max(which(Q[s_t+1,r+1, ] == Vt[s_t+1, r+1])))]
      #note that the first column of optpol is actually corresponding to r = 0 and first row is s = 0
    }
  }
  ## To Return the value function (uncomment this next line:)
  # Vt
  
  ### Return the optimal policy ###
  # *Note that the first column of Optpol is actually corresponding to r = 0 and first row is s = 0*
  # columns refer to r_t from 0 to k 
  # rows refer to s_t from 0 to s_0
  # 1 = Protect and 0 = Not protect
  
  # Optpol ## Uncomment this line to return Optpol
  
  ## pull out the sbar for each r_t
  sbars = apply(Optpol, 2, function(col) {return(suppressWarnings(max(which(col==1))-1))})
  sbars = sbars[-1] # remove meaningless sbar corresponding to r=0
  # if there are cases where no protection is optimal, the lines above will result in infinity
  # since there will be no 1s in that column. rewrite those as zeros to be consistent
  # with Sbar=0 indicating no protection is optimal.
  sbars[is.infinite(sbars)] = 0 
  results.df = data.frame(c = c,
                          discount = discount,
                          a = a,
                          b = b,
                          S_0 = S_0,
                          r_0 = r_0,
                          r_t = 1:(length(sbars)), 
                          sbar = sbars)
  return(list(value.fun = Vt, optimal.policy = results.df))
}
##############################################################################
# Use our numerical solutions function #######################################
###############################################################################
##########################
# Specify Parameters #####
##########################
discount = 0.95 # Discount factor
c = 10 # costs incurred by protection 
# state space limit
S_0 <- 250 # total/maximum number of species
r_0 <- 10 # this is k, the maximum number of species that can provide the service in a pool of size s
### Set parameters for current period payoffs from an ecosystem service ### 
# Power function v(r_t) = a*r_t^b where r_t are the number of species present that can provide the service.
# the b parameter determines the shape of this biodiversity-ecosystem service relationship 
##linear: b = 1 in a power function
##concave or saturating: b< 1
##convex or accelerating: b>1
# Change these parameters to modify the functional relationship between biodiversity and current period payoffs from ecosystem services 
a = 20
# b = .51 #upper estimate from Reich et al (2012) 
# b = .17 #lower estimate from Reich et al (2012) 
b = .26 #O'Connor et al (2016) & Liang et al (2016)
# b = .22 #Cardinale et al (2011)
# b = 1 #linear relationship for comparison 
#b = 2.6 #Mora et al (2014) 
# b = 1.8 #Mora et al (2014) 
################################################################################################
#### Compute solution for different parameter combinations; generate results from figures ######
#########################
### Fig. 5 ##############
########################
results = ComputeNumericalSolutions(a = 20,
                                    b = 0.17,
                                    discount = 0.95,
                                    c = 10,
                                    S_0 = 250,
                                    r_0 = 10)
results2 = ComputeNumericalSolutions(a = 20,
                                     b = 0.22,
                                     discount = 0.95,
                                     c = 10,
                                     S_0 = 250,
                                     r_0 = 10)
results3 = ComputeNumericalSolutions(a = 20,
                                     b = 0.26,
                                     discount = 0.95,
                                     c = 10,
                                     S_0 = 250,
                                     r_0 = 10)
results4 = ComputeNumericalSolutions(a = 20,
                                     b = 0.51,
                                     discount = 0.95,
                                     c = 10,
                                     S_0 = 250,
                                     r_0 = 10)
results5 = ComputeNumericalSolutions(a = 20,
                                     b = 1,
                                     discount = 0.95,
                                     c = 10,
                                     S_0 = 250,
                                     r_0 = 10)
results6 = ComputeNumericalSolutions(a = 20,
                                     b = 1.8,
                                     discount = 0.95,
                                     c = 10,
                                     S_0 = 250,
                                     r_0 = 10)
results7 = ComputeNumericalSolutions(a = 20,
                                     b = 2.6,
                                     discount = 0.95,
                                     c = 10,
                                     S_0 = 250,
                                     r_0 = 10)
results.to.plot = rbind(results$optimal.policy,
                        results2$optimal.policy,
                        results3$optimal.policy,
                        results4$optimal.policy,
                        results5$optimal.policy,
                        results6$optimal.policy,
                        results7$optimal.policy)
results.to.plot$b = as.factor(results.to.plot$b)
ggplot(results.to.plot, aes(x=r_t, y=sbar, group=b, colour=b, fill=b)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  theme_bw() +
  theme(legend.text = element_text(size = 16), legend.title=element_text(size=18)) + 
  ylab("Optimal level of species protection") +
  xlab("Number of critical species (r_t)") +
  scale_x_continuous(breaks=seq(from=0, to = 10, by=1)) +
  theme(axis.title = element_text(size=16)) 
####################################################################
#### Fig. 6 #vary costs to see how results change with b =.26 #####
###################################################################
results = ComputeNumericalSolutions(a = 20,
                                    b = 0.26,
                                    discount = 0.95,
                                    c = 20,
                                    S_0 = 250,
                                    r_0 = 10)
results2 = ComputeNumericalSolutions(a = 20,
                                     b = 0.26,
                                     discount = 0.95,
                                     c = 10,
                                     S_0 = 250,
                                     r_0 = 10)
results3 = ComputeNumericalSolutions(a = 20,
                                     b = 0.26,
                                     discount = 0.95,
                                     c = 5,
                                     S_0 = 250,
                                     r_0 = 10)
results4 = ComputeNumericalSolutions(a = 20,
                                     b = 0.26,
                                     discount = 0.95,
                                     c = 25,
                                     S_0 = 250,
                                     r_0 = 10)
results.to.plot = rbind(results$optimal.policy,
                        results2$optimal.policy, 
                        results3$optimal.policy, 
                        results4$optimal.policy)
results.to.plot$c = as.factor(results.to.plot$c)
ggplot(results.to.plot, aes(x=r_t, y=sbar, group=c, colour=c, fill=c)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  theme_bw() +
  theme(legend.text = element_text(size = 16), legend.title=element_text(size=18)) + 
  ylab("Optimal level of species protection") +
  xlab("Number of critical species (r_t)") +
  scale_x_continuous(breaks=seq(from=0, to = 10, by=1)) +
  theme(axis.title = element_text(size=16)) 
########################################################################
## Fig S5 panel (A) comparing results with different discount factor ###
# to produce FigS5 panel (B), change discount to 0.8 for each below ####
########################################################################
results = ComputeNumericalSolutions(a = 20,
                                    b = 0.17,
                                    discount = 0.99,
                                    c = 10,
                                    S_0 = 250,
                                    r_0 = 10)
results2 = ComputeNumericalSolutions(a = 20,
                                     b = 0.26,
                                     discount = 0.99,
                                     c = 10,
                                     S_0 = 250,
                                     r_0 = 10)
results3 = ComputeNumericalSolutions(a = 20,
                                     b = 1,
                                     discount = 0.99,
                                     c = 10,
                                     S_0 = 250,
                                     r_0 = 10)
results4 = ComputeNumericalSolutions(a = 20,
                                     b = 2.6,
                                     discount = 0.99,
                                     c = 10,
                                     S_0 = 250,
                                     r_0 = 10)
results.to.plot = rbind(results$optimal.policy,
                        results2$optimal.policy,
                        results3$optimal.policy,
                        results4$optimal.policy)
results.to.plot$b = as.factor(results.to.plot$b)
ggplot(results.to.plot, aes(x=r_t, y=sbar, group=b, colour=b, fill=b)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  theme_bw() +
  theme(legend.text = element_text(size = 16), legend.title=element_text(size=18)) + 
  ylab("Optimal level of species protection") +
  xlab("Number of critical species (r_t)") +
  scale_x_continuous(breaks=seq(from=0, to = 10, by=1)) +
  theme(axis.title = element_text(size=16)) 
##############################################################
### Fig S6 panel B: costs > a for different b parameters ####
##############################################################
results = ComputeNumericalSolutions(a = 20,
                                    b = 0.17,
                                    discount = 0.95,
                                    c = 25,
                                    S_0 = 250,
                                    r_0 = 10)
results2 = ComputeNumericalSolutions(a = 20,
                                     b = 0.22,
                                     discount = 0.95,
                                     c = 25,
                                     S_0 = 250,
                                     r_0 = 10)
results3 = ComputeNumericalSolutions(a = 20,
                                     b = 0.26,
                                     discount = 0.95,
                                     c = 25,
                                     S_0 = 250,
                                     r_0 = 10)
results4 = ComputeNumericalSolutions(a = 20,
                                     b = 0.51,
                                     discount = 0.95,
                                     c = 25,
                                     S_0 = 250,
                                     r_0 = 10)
results5 = ComputeNumericalSolutions(a = 20,
                                     b = 1,
                                     discount = 0.95,
                                     c = 25,
                                     S_0 = 250,
                                     r_0 = 10)
results6 = ComputeNumericalSolutions(a = 20,
                                     b = 1.8,
                                     discount = 0.95,
                                     c = 25,
                                     S_0 = 250,
                                     r_0 = 10)
results7 = ComputeNumericalSolutions(a = 20,
                                     b = 2.6,
                                     discount = 0.95,
                                     c = 25,
                                     S_0 = 250,
                                     r_0 = 10)
results.to.plot = rbind(results$optimal.policy,
                        results2$optimal.policy,
                        results3$optimal.policy,
                        results4$optimal.policy,
                        results5$optimal.policy,
                        results6$optimal.policy,
                        results7$optimal.policy)
results.to.plot$b = as.factor(results.to.plot$b)
ggplot(results.to.plot, aes(x=r_t, y=sbar, group=b, colour=b, fill=b)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  theme_bw() +
  theme(legend.text = element_text(size = 16), legend.title=element_text(size=18)) + 
  ylab("Optimal level of species protection") +
  xlab("Number of critical species (r_t)") +
  scale_x_continuous(breaks=seq(from=0, to = 10, by=1)) +
  theme(axis.title = element_text(size=16)) 
##################################################################################################
## Fig S7 panel (A) - vary the a parameter to see how results change with b =.26 #############
# to create Fig S7 panel (B) change costs (c) to 20 #############################################
##################################################################################################
results = ComputeNumericalSolutions(a = 5,
                                    b = 0.26,
                                    discount = 0.95,
                                    c = 10,
                                    S_0 = 250,
                                    r_0 = 10)
results2 = ComputeNumericalSolutions(a = 10,
                                     b = 0.26,
                                     discount = 0.95,
                                     c = 10,
                                     S_0 = 250,
                                     r_0 = 10)
results3 = ComputeNumericalSolutions(a = 20,
                                     b = 0.26,
                                     discount = 0.95,
                                     c = 10,
                                     S_0 = 250,
                                     r_0 = 10)
results4 = ComputeNumericalSolutions(a = 30,
                                     b = 0.26,
                                     discount = 0.95,
                                     c = 10,
                                     S_0 = 250,
                                     r_0 = 10)
results.to.plot = rbind(results$optimal.policy,
                        results2$optimal.policy, 
                        results3$optimal.policy, 
                        results4$optimal.policy)
results.to.plot$a = as.factor(results.to.plot$a)
ggplot(results.to.plot, aes(x=r_t, y=sbar, group=a, colour=a, fill=a)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  theme_bw() +
  theme(legend.text = element_text(size = 16), legend.title=element_text(size=18)) + 
  ylab("Optimal level of species protection") +
  xlab("Number of critical species (r_t)") +
  scale_x_continuous(breaks=seq(from=0, to = 10, by=1)) +
  theme(axis.title = element_text(size=16)) 
#######################################################################################################
## Fig S8: look results for same c/a ratio in Fig. 6 but different a & c values with b =.26 ##########
#######################################################################################################
results = ComputeNumericalSolutions(a = 10,
                                    b = 0.26,
                                    discount = 0.95,
                                    c = 2.5,
                                    S_0 = 250,
                                    r_0 = 10)
results2 = ComputeNumericalSolutions(a = 10,
                                     b = 0.26,
                                     discount = 0.95,
                                     c = 5,
                                     S_0 = 250,
                                     r_0 = 10)
results3 = ComputeNumericalSolutions(a = 10,
                                     b = 0.26,
                                     discount = 0.95,
                                     c = 5.5,
                                     S_0 = 250,
                                     r_0 = 10)
results4 = ComputeNumericalSolutions(a = 10,
                                     b = 0.26,
                                     discount = 0.95,
                                     c = 12.5,
                                     S_0 = 250,
                                     r_0 = 10)
results.to.plot = rbind(results$optimal.policy,
                        results2$optimal.policy, 
                        results3$optimal.policy, 
                        results4$optimal.policy)
results.to.plot$c = as.factor(results.to.plot$c)
ggplot(results.to.plot, aes(x=r_t, y=sbar, group=c, colour=c, fill=c)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  theme_bw() +
  theme(legend.text = element_text(size = 16), legend.title=element_text(size=18)) + 
  ylab("Optimal level of species protection") +
  xlab("Number of critical species (r_t)") +
  scale_x_continuous(breaks=seq(from=0, to = 10, by=1)) +
  theme(axis.title = element_text(size=16)) 
###########################################################################
## References #############################################################
###########################################################################
# Dee et al. To what extent can ecosystem services motivate protecting biodiversity?
# Liang, J., Crowther, T.W., Picard, N., Wiser, S., Zhou, M., Alberti, G., et al. (2016). Positive
# biodiversity-productivity relationship predominant in global forests. Science (80-. )., 354.
# Marescot, L., Chapron, G., Chadès, I., Fackler, P.L., Duchamp, C., Marboutin, E., et al. (2013). Complex
# decisions made simple: a primer on stochastic dynamic programming. Methods Ecol. Evol., 4, 872–884.
# Mora, C., Danovaro, R. & Loreau, M. (2014). Alternative hypotheses to explain why biodiversity-ecosystem
# functioning relationships are concave-up in some natural ecosystems but concave-down in manipulative
# experiments. Sci. Rep., 4, 1–9.
# O’Connor, M.I., Gonzalez, A., Byrnes, J.E.K., Cardinale, B.J., Duffy, J.E., Gamfeldt, L., et al. (2016). A
# general biodiversity-function relationship is mediated by trophic level. Oikos, 1–18.
# Reich, P.B., Tilman, D., Isbell, F., Mueller, K., Hobbie, S.E., Flynn, D.F.B., et al. (2012). Impacts of
# biodiversity loss escalate through time as redundancy fades. Science, 336, 589–92.