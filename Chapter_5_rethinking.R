# Chapter 5: The many variables & the spurious waffles

# Megan Verma 
# 10/3/2022

# A= age at marriage 
# M= marriage rate 
# D= divorce rate (outcome)

# age v. divorce linear regression (univariate) -----

# load data and copy
library(rethinking)
data(WaffleDivorce)
d   <- WaffleDivorce

# standardize variables
d$A <- scale(d$MedianAgeMarriage)
d$D <- scale(d$Divorce)

sd(d$MedianAgeMarriage)
# 1.24 

# compute approximate quadratic posterior dist
m5.1 <- quap(alist(
  D ~ dnorm(mu , sigma) ,
  mu <- a + bA * A ,
  a ~ dnorm(0 , 0.2) ,
  bA ~ dnorm(0 , 0.5) ,
  sigma ~ dexp(1)
) ,
data = d)

# simulate from priors 
# plot lines over the range of 2 standard deviations for both outcome & predictor
set.seed(10)
prior <- extract.prior( m5.1 )
mu <- link( m5.1 , post=prior , data=list( A=c(-2,2) ) )
plot( NULL , xlim=c(-2,2) , ylim=c(-2,2) )
for ( i in 1:50 ) lines( c(-2,2) , mu[i,] , col=col.alpha("black",0.4) )

# posterior predictions 
# compute percentile interval of mean 5.5
A_seq <- seq(from           = -3 ,
             to             = 3.2 ,
             length.out     = 30)
mu <- link(m5.1 , data      = list(A = A_seq))
mu.mean <- apply(mu , 2, mean)
mu.PI <- apply(mu , 2 , PI)

# plot it all
plot(D ~ A ,
     data  = d ,
     col   = rangi2)
lines(A_seq ,
      mu.mean ,
      lwd  = 2)
shade(mu.PI ,
      A_seq)

# marriage v. divorce linear regression (univariate)----
# standardize marriage rate 
d$M  <- scale(d$Marriage)

# model
m5.2 <- quap(alist(
  D ~ dnorm(mu , sigma) ,
  mu <- a + bM * M ,
  a ~ dnorm(0 , 0.2) ,
  bM ~ dnorm(0 , 0.5) ,
  sigma ~ dexp(1)
) ,
data = d)

# simulate from priors 
# plot lines over the range of 2 standard deviations for both outcome & predictor
set.seed(10)
prior <- extract.prior( m5.2 )
mu <- link( m5.2 , post=prior , data=list( M=c(-2,2) ) )
plot( NULL , xlim=c(-2,2) , ylim=c(-2,2) )
for ( i in 1:50 ) lines( c(-2,2) , mu[i,] , col=col.alpha("black",0.4) )

# posterior predictions 
# compute percentile interval of mean 5.5
M_seq <- seq(from           = -3 ,
             to             = 3.2 ,
             length.out     = 30)
mu <- link(m5.2 , data      = list(M = M_seq))
mu.mean <- apply(mu , 2, mean)
mu.PI <- apply(mu , 2 , PI)

# plot it all
plot(D ~ M ,
     data  = d ,
     col   = rangi2)
lines(M_seq ,
      mu.mean ,
      lwd  = 2)
shade(mu.PI ,
      M_seq)
