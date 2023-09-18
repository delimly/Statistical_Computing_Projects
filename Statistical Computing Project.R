# MATH6173 STATISTICAL COMPUTING
# SUPPLEMENTARY (REFERRAL) COURSEWORK
# STUDENT I.D. : 30077966

# The code is split into sections which can be accessed via the switchboard below

# SECTION 1 #######

set.seed(2022)

# Section 1 , Question 1 : Ridge Regression ########

# Part A [5 Marks]: ######

# Functions Purpose:
# To calculate Beta_Hat by providing the close form solution of the  
# optimization problem.

myRidgeReg <- function(X = NULL, Y = NULL, lambda = 0, scale = FALSE){
  
  meanColX <- NA 
  sdColX <- NA
  
  # Conditional Statement for response to Scale Object
  
  if(scale == TRUE){
    
    #Standardizing Both X and Y
    
    Y <- scale(Y)
    X <- scale(X)
    
  }else if(scale == FALSE){
    
    X <- X
    Y <- Y
    
  }else{
    
    stop("Incorect value entered. Please, enter TRUE or FALSE")
  }
  
  # Carrying Out the Calculations of the closed form solution to the 
  # Optimization Problem

  
  coefs <- solve(t(X) %*% X + lambda * diag(ncol(X))) %*% t(X) %*% Y
  
  # Converting to a vector
  
  coefs <- as.vector(coefs)
  
  return(coefs)
}

# Part B [2 Marks]:#####

#Training data ; with n = 40, p = 6

X <- matrix(rnorm(40*6),ncol=6) # X matrix

beta <- c(-5,-3,-1, 1, 3, 5) # true coefficients

error <- rnorm(40, 0, 13) # errors

Y <- X %*% beta + error # Y vector

# Carrying Out the Function

Test1 <- myRidgeReg(X = X, Y = Y)

# Fitting 

XY.df <- data.frame(X, Y)

fit <- lm(Y ~ X - 1, data = XY.df)

summary(fit)

Test1

# Upon checking.  they are the same results (To three decimal places) this is expected as 
# the lm function carry out 
# a similar function in regards to what we have calculated. 

# Part C [3 Marks]:######

lambda.vec <- seq(0,200, length.out = 41)

# First Plot

p0.vec <- rep(0,6)

beta0.vec <- myRidgeReg(X=X, Y=Y, lambda = lambda.vec[1])

plot(p0.vec, beta0.vec, 
            main = "The relationship between \n Beta Estimates and Lambda",
            xlab = "lambda", 
            ylab = "Beta Estimates",
            pch = 20, 
            col = c(1,2,3,4,5,6),
            xlim = c(0,200))

# Subsequent Plots
 for(i in 2:41){
   plambda.vec <- rep(lambda.vec[i], 6)
   
   betalambda.vec <- myRidgeReg(X=X, Y=Y, lambda = lambda.vec[i])
   
   points(plambda.vec, betalambda.vec, 
          pch = 20,
          col = c(1,2,3,4,5,6),
          xlim = c(lambda.vec[0],lambda.vec[i]),
          ylim = c(min(betalambda.vec),max(betalambda.vec))
          )
 }

legend(175,5, pch = rep(20,6), col = c(1,2,3,4,5,6), 
       legend = c("Beta 1", "Beta 2", "Beta 3", "Beta 4", "Beta 5", "Beta 6"), cex = 0.66)


# Part D [3 Marks]:#####


lambda.vec # sequence Vector
SPE.vec <- NA #SPE VECTOR
AVR.SPE.vec <- NA # AVR SPE VECTOR

for (i in 1:41) {
  
  beta.vec <- myRidgeReg(X=X, Y=Y, lambda = lambda.vec[i])
  
  for(j in 1:100){
    for(k in 1:10){
      #Test Data

      # Using m=10, p= 6 
      
      X.star <- matrix(rnorm(10*6),ncol=6) # X* matrix
      
      error.star <- rnorm(10, 0, 13) # errors
      
      Y.star <- X.star %*% beta + error.star # Y* vector
      
      # Calculating the  Squared Mean Error Vector of Length 100
      
      SPE.vec[j] <- mean((Y.star[k] - X.star[k, ]*beta.vec)^2) 
      
      # Calculating the Average Mean Squared Error Vector of Length 41
      
      AVR.SPE.vec[i] <- mean(SPE.vec)
      }
    }
  }
 

plot(lambda.vec, AVR.SPE.vec, type = "l", lwd = 2.25, main = "The relationship between Lambda 
     and the averaged SPE", ylab = "Avaraged SPE", xlab = "Lambda")

# Finding the Optimum Values

lambda.vec[which.min(AVR.SPE.vec)]
min(AVR.SPE.vec)

abline(h = min(AVR.SPE.vec), col ="Blue", lwd = 2.25)
abline(v = lambda.vec[which.min(AVR.SPE.vec)], col = "red", lwd = 2.25)



# Section 1 , Question 2 : A-R, Importance and M-H sampling ########

# Part A [5 Marks]: ########

# Functions Purpose
# To implement the A-R Sampling Algorithm and be able to return a sequence of 
# n-number of samples 

r.myHFLN.ar <- function(a = NULL, n = NULL){
  
  X <- NULL
  
  for(i in 1:n){
  
    repeat{
      x <- rlnorm(1)
      y <- runif(1)
      
      if(y <= (1 + a*sin(2 * pi * log(x)))/(1 + a)){
      
        X[i] <- x
        break}
    }
  }
  return(X)
  }
  
# Producing 2000 Samples
samples1 <- r.myHFLN.ar(0.7 , 2000)
  
# The Monte Carlo Estimation of P(X > 4)
mean(samples1>4) 


# Part B [5 Marks]: ##########

# Functions Purpose
# To be able to implement the Importance Sampling Algorithm to be ablr to 
# calculate an estimate of the tail probabilty

r.myHFLN.is <- function(a = NULL, n = NULL){
  
  # Producing samples 
  # Calculating the importance weights vector and excluding negative samples 
  
  X <- rlnorm(n)
  X[X <= 0] <- 0
  
  # Calculating the Weights
  
  v.weights <- (1 + a * sin( 2 * pi * log(X)))*(X > 0)
  
  # Calculating the Estimator
  
  p.t <- mean(v.weights*(X>4))
  return(p.t)
}

# Tail probability estimator

r.myHFLN.is(0.7, 2000)

# Part C [5 Marks]: ########

# Functions Purpose 
# To be able to implement the Metropolis Hastings Algorithm to generate a 
# sequence of n number of samples 

r.myHFLN.mh <- function(x0 = NULL, a = NULL, n = NULL, m = NULL){
  
  X <- rep(NA, n + m + 1)
  
  # Initial Value
  
  X[1] <- x0
  
  for (j in 2 : (n + m + 1)) {
    Y <- rlnorm(1)
    
    alpha <- min((1 + a * sin(2 * pi * log(Y)))/(1 + a * sin(2 * pi * log(X[j-1])))
                  ,1)
    
    if(runif(1) <= alpha){
      
      X[j] <- Y
    }else{
      X[j] <- X[j-1]
    }
  }
  return(X[(m+2): (n+m+1)])
}

# Producing 2000 Samples
samples2 <- r.myHFLN.mh(0.8, 0.7, 2000, 500)

# The Monte Carlo Estimation of P(X > 4)
mean(samples2>4)

# Section 1 , Question 3 : Inverse transformation and Gibbs sampler ########

# Part A [3 Marks]: ######

# Answer Is Within The Report 

# Part B [3 Marks]:######

# Answer Is Within The Report

# Part C [3 Marks]: #########

# Functions Purpose
# To produce one random sample based upon the condition provided by the argument
# the Function uses the inverse transformation algorithm to produce a sample


r.myconditional <- function(con = NULL){
  
  xc <- con
  yc <- con
  
  # Generating X, given Y = y
  
  u1 <- runif(1)
  x <- - (yc+2) + sqrt((yc+2)^2 + u1*(2*yc + 3))
  
  # Generating Y, given X = x
  
  u2 <- runif(1)
  y <- - (xc+2) + sqrt((xc+2)^2 + u2*(2*xc + 3))
  
  #one random sample of X and Y as a vector
  
  X <- c(x,y)
  
  return(X)
}
  
  
# Part D [3 Marks]:

# Functions Purpose
# To be able to carry out the Gibbs Sampler Algorithm to produce a sequrnce of
# n random samples. 

gibbs.bipoly <- function(x0 = NULL, y0 = NULL
                          , n = NULL, m = NULL){
  
  x.seq<-y.seq<-rep(NA,m+n+1)
  
  x.seq[1] <- x0
  y.seq[1] <- y0 
  
  for (i in 2:(m+n+1)) {
    
    x.seq[i] <- r.myconditional(con = y.seq[i-1])[1]
    
    y.seq[i] <- r.myconditional(con = x.seq[i])[2]
    
  }
  
  result <- list(X = x.seq[(m+2):(m+n+1)], Y = y.seq[(m+2):(m+n+1)])
  
  return(result)
} 

# Part E [2 Marks]: ##########



samples3 <- gibbs.bipoly(0.5, 0.5, 1000, 500)

par(mfrow = c(1,2))

points = seq(0,1, length.out = 2000)

plot(density(samples3$X), 
     main = "simulated and actual pdf of X",
     xlab = "X", lwd = 3, col = "blue", xlim = c(-0.5, max(points)))

lines(x = points, y = 2 * (2 * points + 1)/(2 * points + 3), lwd = 2, lty = 3,
      col = "red")

plot(density(samples3$Y) , 
     main = "simulated and actual pdf of Y",
     xlab = "Y", lwd = 3, col = "blue" , xlim = c(-0.5, max(points)))

lines(x = points, y = 2 * (2 * points + 1)/(2 * points + 3), lwd = 2, lty = 3,
      col = "red")

dev.off() # Shut me off if you can

# Section 1 , Question 4 : Sampling N(0,1) ########

# Part A [4 Marks]: ########

# Answer is inside the report

# Part B [3 Marks]:

r.mynormal.polar <- function(n = NULL){
  
  # LCG for Uniform Distribution of Theta U(0,2 pi)
  
  # Initialization 
  
  M <- 2^32
  Y0 <- 2022
  a <- 07432281838
  b <- 12345
  Y <- NULL
  Theta  <- NULL
  X1 <- NULL
  X2 <- NULL
  
  # Loops and Conditional Statements to get values of Theta, X1, X2
  
  for (i in 1: n){ 
    if (i==1){
    Y[i] <- (a*Y0 +b) %% M 
    
    # Case of the first element
    
    Theta[1] <- Y[1]/(2*pi *M)
    
    # Sampling radius^2 
    
    R_squared <- rchisq(1, df =2)
    
    # Calculating radius 
    
    R <- sqrt(R_squared)
    
    # X1[1] and x2[1] pairs for 1st element
    
    X1[1] <- R* cos(Theta[1])
    X2[1] <- R* sin(Theta[1])
    
    }else {
      
      # Subsequent Samples 
      
      # X1[1] and x2[1] pairs for ith element
      
      Y[i] <- (a*Y[i-1] +b) %% M
      
      Theta[i] <- Y[i]/(2*pi *M)
      
      # Sampling radius^2 
      
      R_squared <- rchisq(1, df =2)
      
      # Calculating radius 
      
      R <- sqrt(R_squared)
      
      X1[i] <- R* cos(Theta[i])
      X2[i] <- R* sin(Theta[i])
    }}
    
  # Combining X1 and X2 Vectors to 
  
  X <- rbind(X1, X2)
    X <- as.vector(X)
    return(X)
  }

# Part C [2 Marks]: 

# Samples

sample4 <- r.mynormal.polar(1000)

# Plots , adding curve and Legend

plot(density(sample4),
     main = "Empirical and actual pdfs of N(0,1)",
     lwd=3,
     col = "blue")

curve(dnorm(x), add = T, lwd = 3, lty = 3, col = "red")

legend(x=3,y=0.6 ,legend = c("Empirical pdf","Actual pdf"), lty=c(1,3),col=c("blue","red"),cex=1)


# SECTION 2 #########


# Section 2, Question 1 :  Integration by the trapezoidal rule ##########

# Part A [6 Marks]########

# Functions Purpose 
# To calculate the integral of a cubic spline function via the trapezium rule

cubSP1 <- function(x = NULL){
  
  # Creating our function vector, f[i] will be replaced by f(x[i]) as outlined below
  
  f <- rep(0, length(x))
  
  # Conditional Statements and Loop to be able to carry out cubic spline function
   
 for (i in 1:length(x)){
   
   # Here we calculate by vector element to replace the values in our function vector
   
   if(x[i] <= -1.75){
     
     f[i] <- -0.25 + 0.45 * x[i] + 0.15 * (x[i])^2 + 0.01 * (x[i])^3
     
   }else if(-1.75 < x[i] && x[i] <= 2.25){
     
     f[i] <- -0.25 + 0.45 * x[i] + 0.15 * (x[i])^2 + 0.01 * (x[i])^3 - 0.02 * (x[i] + 1.75)^3
     
   }else{
     
     f[i] <- -0.25 + 0.45 * x[i] + 0.15 * (x[i])^2 + 0.01 * (x[i])^3 - 0.02 * (x[i] + 1.75)^3 - 0.025 * (x[i] - 2.25)^3
   }
 } 
  
return(f)
  
}

# Part B [7 Marks]#######

#Functions Purpose
# To calculate the Area of the cubic spline function via trapezium rule

calcCubSplArea <- function(xVec = NULL){
  
  
  
  Area <- 0
  Int <- 0
  
  f <- cubSP1(xVec)
  
  f[which(f < 0)] <- f[which(f < 0)]*(-1)
  
  for (j in 2:length(xVec)) {
    Int <- (xVec[j] - xVec[j-1])
    Area <- Area + (( f[j] + f[j-1])/2) * Int 
    
  }
  
  return(Area)
  
}

# Testing to see algorithm works
testX <- seq(-5.5, 5.5, by = 0.1)
testX

length(testX)

cubSP1(testX) # function vector 

calcCubSplArea(testX) # Area 

# Part C [4 Marks]######

# Functions Purpose
# To carry out the function in Part A but for a generic case


cubSP1Gen <- function(x = NULL, beta = NULL, knots = NULL){
  
  # Conditionals to ensure right arguments are provided
  
  if(length(beta) != 6){
    stop("Wrong number of beta values entered please try again")
  }
  
  if(length(knots) != 2){
    stop("Wrong number of knots values entered please try again")
  }
  
  # Creating our function vector, f[i] will be replaced by f(x[i]) as outlined below
  
  f <- rep(0, length(x))
  
  # Conditional Statements and Loop to be able to carry out cubic spline function

  
  for (i in 1:length(x)){
    
    # Here we calculate by vector element to replace the values in our function vector
    
    if(x[i] <= knots[1]){
      
      f[i] <- beta[1] + beta[2] * x[i] + beta[3] * (x[i])^2 + beta[4] * (x[i])^3
      
    }else if(knots[1] < x[i] && x[i] <= knots[2]){
      
      f[i] <- beta[1] + beta[2] * x[i] + beta[3] * (x[i])^2 + beta[4] * (x[i])^3 - beta[5] * (x[i] - knots[1])^3
    }else{
      
      f[i] <- beta[1] + beta[2] * x[i] + beta[3] * (x[i])^2 + beta[4] * (x[i])^3 - beta[5] * (x[i] - knots[1])^3 - beta[6] * (x[i] - knots[2])^3
    }
  } 
  
  return(f)

}



# Part D [5 Marks]#######

# Functions Purpose
# To carry out the function in part b but for a generic case.

calcCubSplGenArea <- function(xVec = NULL, beta = NULL, knots = NULL){
  
  # Calculate the F_x
  
  Area <- 0
  Int <- 0
  
  f <- cubSP1Gen(xVec, beta = beta, knots = knots)
  
  f[which(f < 0)] <- f[which(f < 0)]*(-1)
  
  for (j in 2:length(xVec)) {
    Int <- (xVec[j] - xVec[j-1])
    Area <- Area + (( f[j] + f[j-1])/2) * Int 
    
  }
  
  return(Area)
  
}
  

# Part E [3 Marks] #########

#Maintaining values of Xvec, beta and knots

testX <- seq(-5.5, 5.5, by = 0.1)
testX
length(testX)

testBeta <- c(-0.5, 1.2, 0.35, 0.01, -0.04, -0.03)
testBeta

testKnots <- c(-2,2)
testKnots

# Calculating the Area 

cubSP1Gen(testX, testBeta, testKnots) # function vector

calcCubSplGenArea(testX, testBeta, testKnots) # This gives the value of the integral


# Section 2, Question 2: EM algorithm for a Mixture of Three Gamma Distributions ##########

# One Question : [10 Marks] ######

# This Question Was answered in the report. 

# Section 2, Question 3 : Bootstrapping for Estimation of Uncertainty ##########

# One Question : [15 Marks] #####

bootStats <- function(x = NULL, w = NULL, bootCount = NULL){
  
  # Conditional Statements to ensure the algorithm runs smoothly
  
  if(length(x) != length(w)){
    
    stop("The vectors are not of same length. Please, enter vectors of same length")
    
  }
 
  if(any(w < 0 & 1 < w)){
    
    stop("The weight vector has values that cannot be accepted. Please enter a 
         vector whose values are between 0 and 1 (inclusive)")
  }
  
  # The matrix that will contain the bootstrap estimates 
  
  bootCoefs <- matrix(nrow = bootCount, ncol = 2)
  
  # Looping through each bootstrap estimate
  
  for (repIndex in 1:bootCount) {
    
    # Creating bootstrap samples 
    
    repSample <- sample(1:length(x), size = length(x), replace = T)
    
    xB <- x[repSample]
    wB <- w[repSample]
    
    # Calculating the weighted means
    
    meanWtdxB <- mean(wB * xB)
    
    # Calculating the weighted variance 
    
    varWtdxB <- sum(((wB*xB)-meanWtdxB)^2)/(bootCount-1)
    
    # The bootstrap estimate of beta
    
    betaBoot <- meanWtdxB/varWtdxB
    
    # The bootstrap estimate of alpha
    
    alphaBoot <- meanWtdxB * betaBoot
    
    bootCoefs[repIndex,] <- c(alphaBoot, betaBoot)
    
  }
  
  # Calculating the standard error
  
  bootSE <- apply(bootCoefs, 2, function(bootCoef){
    sqrt(sum((bootCoef - mean(bootCoef))^2)/(bootCount - 1))
    })
  
  bootSE <- as.vector(bootSE)
  
  return(bootSE)
}

#This is the end of the script 