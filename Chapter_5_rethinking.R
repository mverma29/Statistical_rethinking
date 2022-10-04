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




















# make a DAG -----
# install.packages("daggity")
library(dagitty)

# (all variables related)
dag5.1 <- dagitty("dag {
A -> D
A -> M
M -> D
}")

coordinates(dag5.1) <- list(x = c(A = 0, D = 1, M = 2) , y = c(A = 0, D =
                                                                 1, M = 0))
drawdag(dag5.1)

# (no relationship between M & D)
DMA_dag2 <- dagitty('dag{ D <- A -> M }')
impliedConditionalIndependencies(DMA_dag2)

# check the conditional dependences of first dag DMA_dag2 <- dagitty('dag{ D <- A -> M }')
DMA_dag1 <- dagitty('dag{ D <- A -> M -> D }')
impliedConditionalIndependencies(DMA_dag1)
# no conditional independencies, so no output

# multivariate regression (A,M,D)----

# approximate posterior distribution 
m5.3 <- quap(alist(
  D ~ dnorm(mu , sigma) ,
  mu <- a + bM * M + bA * A ,
  a ~ dnorm(0 , 0.2) ,
  bM ~ dnorm(0 , 0.5) ,
  bA ~ dnorm(0 , 0.5) ,
  sigma ~ dexp(1)
) ,
data = d)
precis(m5.3)

# plot the posterior dist for all 3 models (to compare)
# only slopes for bM & bA
plot( coeftab(m5.1,m5.2,m5.3), par=c("bA","bM") )


# posterior residual plots----
# model for marriage rate using age at marriage 
m5.4 <- quap(alist(
  M ~ dnorm(mu , sigma) ,
  mu <- a + bAM * A ,
  a ~ dnorm(0 , 0.2) ,
  bAM ~ dnorm(0 , 0.5) ,
  sigma ~ dexp(1)
) ,
data = d)

# compute residuals by subtracting the observed marriage 
# rate (in each state) from the predicted rate
mu <- link(m5.4)
mu_mean <- apply( mu , 2 , mean )
mu_resid <- d$M - mu_mean

# posterior prediction plots----
# call link without specifying new data so it uses original data
mu <- link(m5.3)

# summarize samples across cases
mu_mean <- apply(mu , 2 , mean)
mu_PI <- apply(mu , 2 , PI)

# simulate observations
# again no new data, so uses original data
D_sim <- sim(m5.3 , n = 1e4)
D_PI <- apply(D_sim , 2 , PI)

# plot predictions against observed 
# add a line to show perfect prediction
# add line segments for the CI of each prediction
plot(
  mu_mean ~ d$D ,
  col  = rangi2 ,
  ylim = range(mu_PI) ,
  xlab = "Observed divorce" ,
  ylab = "Predicted divorce"
)
abline(a = 0 , b = 1 , lty = 2)
for (i in 1:nrow(d))
  lines(rep(d$D[i], 2) , mu_PI[, i] , col = rangi2)

# label a few points (the ones I click on)
identify( x=d$D , y=mu_mean , labels=d$Loc )

# counterfactual plots-----

# use the DAG with all arrows: 

# (all variables related)
drawdag(dag5.1)

# set of functions to tell how each variable is generated 
# use Gaussian dist for each variable, for simplicity
# to estimate the influence of A on M, we regress A on M 
  # just add this to quap model 

data(WaffleDivorce)
d <- list()
d$A <- standardize(WaffleDivorce$MedianAgeMarriage)
d$D <- standardize(WaffleDivorce$Divorce)
d$M <- standardize(WaffleDivorce$Marriage) 

m5.3_A <- quap(
  alist(
    ## A -> D <- M
    D ~ dnorm(mu , sigma) ,
    mu <- a + bM * M + bA * A ,
    a ~ dnorm(0 , 0.2) ,
    bM ~ dnorm(0 , 0.5) ,
    bA ~ dnorm(0 , 0.5) ,
    sigma ~ dexp(1),
    ## A -> M
    M ~ dnorm(mu_M , sigma_M),
    mu_M <- aM + bAM * A,
    aM ~ dnorm(0 , 0.2),
    bAM ~ dnorm(0 , 0.5),
    sigma_M ~ dexp(1)
  ) ,
  data = d
)

precis(m5.3_A)

# simulate 

# define a range of values for A (to manipulate so M is reduced)
A_seq <- seq(from       = -2 ,
             to         = 2 ,
             length.out = 30) #30 obs between |2| SD from mean

# prep data
sim_dat <- data.frame(A = A_seq)

# simulate M and then D, using A_seq
s <- sim(m5.3_A , 
         data = sim_dat , 
         vars = c("M", "D"))

# display counterfactual predictions (A + M on D)
plot(
  sim_dat$A ,
  colMeans(s$D) ,
  ylim = c(-2, 2) ,
  type = "l" ,
  xlab = "manipulated A" ,
  ylab = "counterfactual D"
)
shade(apply(s$D, 2, PI) , sim_dat$A)
mtext("Total counterfactual effect of A on D")

# display counterfactual predictions (A on M)
plot(
  sim_dat$A ,
  colMeans(s$M) ,
  ylim = c(-2, 2) ,
  type = "l" ,
  xlab = "manipulated A" ,
  ylab = "counterfactual M"
)
shade(apply(s$M, 2, PI) , sim_dat$A)
mtext("Total counterfactual effect of A on M")

# new counterfactual-- manipulate M 
# new DAG bc if we manipulate M, A doesn't influence it
DMA_dag3 <- dagitty('dag{ M -> D <- A }')
coordinates(DMA_dag3) <- list(x = c(A = 0, D = 1, M = 2) , y = c(A = 0, D =
                                                                 1, M = 0))
drawdag(DMA_dag3)

# simulate counterfactual result for an average state (A=0) & change M
sim_dat <- data.frame(M = seq(
  from       = -2,
  to         = 2,
  length.out = 30
) , A        = 0)

s <- sim(m5.3_A , data = sim_dat , vars = "D")

plot(
  sim_dat$M ,
  colMeans(s) ,
  ylim = c(-2, 2) ,
  type = "l" ,
  xlab = "manipulated M" ,
  ylab = "counterfactual D"
)
shade(apply(s, 2, PI) , sim_dat$M)
mtext("Total counterfactual effect of M on D")


# masked relationship---- 

library(rethinking)
data(milk)
d <- milk
str(d)

# standardize vars in question
# outcome is kcal per gram milk  
d$K <- scale(d$kcal.per.g) # kcal per g of milk 
d$N <- scale(d$neocortex.perc) # neocortex percent of brain mass
d$M <- scale(log(d$mass)) # avg female body mass

# first try to run the model with vague priors 
m5.5_draft <- quap(alist(
  K ~ dnorm(mu , sigma) ,
  mu <- a + bN * N ,
  a ~ dnorm(0 , 1) ,
  bN ~ dnorm(0 , 1) ,
  sigma ~ dexp(1)
) ,
data = d)
# error message-- missing values in N variable 
# complete case analysis 
dcc <- d[complete.cases(d$K, d$N, d$M) ,]

# try model again using cc analysis data frame
m5.5_draft <- quap(alist(
  K ~ dnorm(mu , sigma) ,
  mu <- a + bN * N ,
  a ~ dnorm(0 , 1) ,
  bN ~ dnorm(0 , 1) ,
  sigma ~ dexp(1)
) ,
data = dcc)

# simulate & plot 50 prior regression lines 
prior <- extract.prior(m5.5_draft)
xseq <- c(-2, 2)
mu <- link(m5.5_draft , post = prior , data = list(N = xseq))
plot(NULL , xlim = xseq , ylim = xseq)
for (i in 1:50)
  lines(xseq , mu[i, ] , col = col.alpha("black", 0.3))
# crazy priors! 

# try again, making a prior closer to 0 
# (so outcome will be closer to 0) and tighten
# slope of bN so it doesn't produce impossibly strong relationships
m5.5 <- quap(alist(
  K ~ dnorm(mu , sigma) ,
  mu <- a + bN * N ,
  a ~ dnorm(0 , 0.2) ,
  bN ~ dnorm(0 , 0.5) ,
  sigma ~ dexp(1)
) ,
data = dcc)

# plot these revised priors
prior <- extract.prior(m5.5)
xseq <- c(-2, 2)
mu <- link(m5.5 , post = prior , data = list(N = xseq))
plot(NULL , xlim = xseq , ylim = xseq)
for (i in 1:50)
  lines(xseq , mu[i, ] , col = col.alpha("black", 0.3))
# better (but still vague) priors
# stays within high prob region of observable data 

# look at model output table
precis(m5.5)

# plot predicted mean and 89% compatibility interval for mean
xseq    <-
  seq(
    from       = min(dcc$N) - 0.15 ,
    to         = max(dcc$N) + 0.15 ,
    length.out = 30
  )
mu      <- link(m5.5 , data = list(N = xseq))
mu_mean <- apply(mu, 2, mean)
mu_PI   <- apply(mu, 2, PI)

plot(K ~ N , data = dcc)
lines(xseq , mu_mean , lwd = 2)
shade( mu_PI , xseq )
# weakly positive post mean
# highly imprecise 

# now, let's consider bivariate relationship btw. kcal & body mass
m5.6 <- quap(alist(
  K ~ dnorm(mu , sigma) ,
  mu <- a + bM * M ,
  a ~ dnorm(0 , 0.2) ,
  bM ~ dnorm(0 , 0.5) ,
  sigma ~ dexp(1)
) ,
data = dcc)

precis(m5.6)

# plot predicted mean and 89% compatibility interval for mean
xseq    <-
  seq(
    from       = min(dcc$M) - 0.15 ,
    to         = max(dcc$M) + 0.15 ,
    length.out = 30
  )
mu      <- link(m5.6 , data = list(M = xseq))
mu_mean <- apply(mu, 2, mean)
mu_PI   <- apply(mu, 2, PI)

plot(K ~ M , data = dcc)
lines(xseq , mu_mean , lwd = 2)
shade( mu_PI , xseq )


# add BOTH predictor variables at same time to regression 
# approx posterior dist
m5.7 <- quap(alist(
  K ~ dnorm(mu , sigma) ,
  mu <- a + bN * N + bM * M ,
  a ~ dnorm(0 , 0.2) ,
  bN ~ dnorm(0 , 0.5) ,
  bM ~ dnorm(0 , 0.5) ,
  sigma ~ dexp(1)
) ,
data = dcc)

precis(m5.7)

# compare this posterior to that of models m5.5 and m5.6
plot(coeftab(m5.5 , m5.6 , m5.7) , pars = c("bM", "bN"))

# are N & M correlated? 
pairs (~K + M + N, dcc) #yes 

# make counterfactual plots (shows how model sees problem)

# holding N=0
xseq <-
  seq(
    from       = min(dcc$M) - 0.15 ,
    to         = max(dcc$M) + 0.15 ,
    length.out = 30
  )
mu <- link(m5.7 , data = data.frame(M = xseq , N = 0))
mu_mean <- apply(mu, 2, mean)
mu_PI <- apply(mu, 2, PI)

plot(NULL , xlim = range(dcc$M) , ylim = range(dcc$K))
lines(xseq , mu_mean , lwd = 2)
shade(mu_PI , xseq)

# holding M=0
xseq <-
  seq(
    from       = min(dcc$N) - 0.15 ,
    to         = max(dcc$N) + 0.15 ,
    length.out = 30
  )
mu <- link(m5.7 , data = data.frame(N = xseq , M = 0))
mu_mean <- apply(mu, 2, mean)
mu_PI <- apply(mu, 2, PI)

plot(NULL , xlim = range(dcc$N) , ylim = range(dcc$K))
lines(xseq , mu_mean , lwd = 2)
shade(mu_PI , xseq)



