#mdp_critical_state_minus <- function(P,R) {
  library(MDPtoolbox)
  source("mdp_example_rand.r")
  source("mdp_eval_policy_iterative_q.r")
  source("mdp_value_iteration.r")
  source("mdp_bellman_operator_q.r")  
  PR=mdp_example_rand(10,4) # S and A
  P = PR$P
  R = PR$R
  discount = 0.96
  
  sol=mdp_value_iteration(P,R,discount,0.001)
  Pol = sol$policy
  
  
  # need to write these 2 functions in R
  y = mdp_eval_policy_iterative_q(P, R, discount, Pol); 
  # check that this is doing what it is supposed to; Get the value and Q values
  Q <- y$Q
  Vpolicy <- y$Vpolicy
  x <- dim(Q)
  nbS <- x[1]
  nbA <- x[2]
  LP <- array(0,c(nbS,nbA))
  MinusLP <- array(0,c(nbS))
  PlusLP <- array(0,c(nbS))
  for (i in 1:nbS){
    # for all states
    optA <- Pol[i]
    SetA <- setdiff(1:nbA, optA)
    LP[i,optA] <- NA
    for (a in SetA){
      LP[i,a] <- (Vpolicy[i] - Q[i,a])/Vpolicy[i]
    }
    MinusLP[i] <- min(LP[i,], na.rm= TRUE)
    PlusLP[i] <- max(LP[i,], na.rm= TRUE)
  }
  

  
  #plot(MinusLP)
  #plot(PlusLP)
  # Calculate the range for y-axis
  yrange <- range(c(MinusLP, PlusLP)) * 100
  
  plot(1:nbS, MinusLP*100, type="p", ylim=yrange, xlab='States',ylab='Loss (vs optimal) (%)', main="Plot of MinusLP and PlusLP")
  points(1:nbS,PlusLP*100, col="red")
  # Draw lines between points
  segments(x0=1:nbS, y0=MinusLP*100, x1=1:nbS, y1=PlusLP*100, col="blue")
  legend("topright", legend = c("LP- (best case)", "LP+ (worst case)"), col = c("black", "red"), pch = 1)

  
 # }