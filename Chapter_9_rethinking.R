# Chp 9 Exercises

# Megan Verma
# 10/11/2022

library(rethinking)

# Markov's island hopping ex:----
num_weeks <- 1e5
positions <- rep(0, num_weeks)
current <- 10
for (i in 1:num_weeks) {
  # record current position
  positions[i] <- current
  # flip coin to generate proposal
  proposal <- current + sample(c(-1, 1) , size = 1)
  # now make sure he loops around the archipelago
  if (proposal < 1)
    proposal <- 10
  if (proposal > 10)
    proposal <- 1
  # move?
  prob_move <- proposal / current
  current <- ifelse(runif(1) < prob_move , proposal , current)
}

plot(1:100 , positions[1:100])
plot(table(positions))

# sample randomly from a high dimension dist (01 dimensions is enough)
# plot radial dist of the points
D <- 10 
T <- 1e3
Y <- rmvnorm(T, rep(0, D), diag(D))
rad_dist <- function(Y)
  sqrt(sum(Y ^ 2))
Rd <- sapply(1:T , function(i)
  rad_dist(Y[i, ]))
dens(Rd)

# HMC ulam example:----

# return to terrain ruggedness example 
library(rethinking)
data(rugged)
d <- rugged
# reduce to cases (nations) that have outcome variable of interest 
d$log_gdp <- log(d$rgdppc_2000)
dd <- d[complete.cases(d$rgdppc_2000) ,]
# standardize GDP, center ruggedness
dd$log_gdp_std <- dd$log_gdp / mean(dd$log_gdp)
dd$rugged_std <- dd$rugged / max(dd$rugged)
# index variable for africa (or not)
dd$cid <- ifelse(dd$cont_africa == 1 , 1 , 2)

# old model-- predict log GDP with terrain ruggedness, continent & interaction
m8.3 <- quap(alist(
  log_gdp_std ~ dnorm(mu , sigma) ,
  mu <- a[cid] + b[cid] * (rugged_std - 0.215) ,
  a[cid] ~ dnorm(1 , 0.1) ,
  b[cid] ~ dnorm(0 , 0.3) ,
  sigma ~ dexp(1)
) ,
data = dd)
precis(m8.3 , depth = 2)

# now do same model using Hamiltonian monte carlo 

# first, pre-transform variables
# we've already done this 

# then, make a new list with only what's in model 
dat_slim <- list(
  log_gdp_std = dd$log_gdp_std,
  rugged_std = dd$rugged_std,
  cid = as.integer(dd$cid)
)
str(dat_slim)

# sampling from the posterior using ulam & stan
m9.1 <- ulam(
  alist(
    log_gdp_std ~ dnorm(mu , sigma) ,
    mu <- a[cid] + b[cid] * (rugged_std - 0.215) ,
    a[cid] ~ dnorm(1 , 0.1) ,
    b[cid] ~ dnorm(0 , 0.3) ,
    sigma ~ dexp(1)
  ) ,
  data = dat_slim ,
  chains = 1
)

precis(m9.1 , depth = 2)

# try again, this time running 4 chains across the 4 cores of the computer 
m9.1 <- ulam(
  alist(
    log_gdp_std ~ dnorm(mu , sigma) ,
    mu <- a[cid] + b[cid] * (rugged_std - 0.215) ,
    a[cid] ~ dnorm(1 , 0.1) ,
    b[cid] ~ dnorm(0 , 0.3) ,
    sigma ~ dexp(1)
  ) ,
  data = dat_slim ,
  chains = 4 ,
  cores = 4
)

# how long did each chain take?
show(m9.1)

precis(m9.1 , 2)

# pairs plot shows hte smoothed histogram of each parameter
pairs(m9.1)

# trace plot of the model
traceplot(m9.1)
# explore the chains individually
traceplot(m9.1, chains=4)

# trace-rank plot (trank plot)
trankplot(m9.1 , n_cols = 2)




# taming a wild chain example: -----
y <- c(-1, 1)
set.seed(11)

m9.2 <- ulam(
  alist(
    y ~ dnorm(mu , sigma) ,
    mu <- alpha ,
    alpha ~ dnorm(0 , 1000) ,
    sigma ~ dexp(0.0001)
  ) ,
  data = list(y = y) ,
  chains = 3
)

precis(m9.2) # this is crazy 

traceplot(m9.2) # sick Markov chains
trankplot(m9.2) 

# try again with weakly informative priors
set.seed(11)
m9.3 <-  ulam(
  alist(
    y ~ dnorm(mu , sigma) ,
    mu <- alpha ,
    alpha ~ dnorm(1 , 10) ,
    sigma ~ dexp(1)
  ) ,
  data = list(y = y) ,
  chains = 3
)
precis(m9.3)
traceplot(m9.3) # better Markov chains
trankplot(m9.3) 



# non-identifiable ex:----

# first sim 100 obs from a Gaussian dist with mean zero and SD 1 
set.seed(41)
y <- rnorm(100 , mean = 0 , sd = 1)

# fit the linear model (2 parameters)
# only the sum of the parameters can be identified 

# run the Markov chain 
set.seed(384)
m9.4 <- ulam(
  alist(
    y ~ dnorm(mu , sigma) ,
    mu <- a1 + a2 ,
    a1 ~ dnorm(0 , 1000),
    a2 ~ dnorm(0 , 1000),
    sigma ~ dexp(1)
  ) ,
  data = list(y = y) ,
  chains = 3
)
precis(m9.4)
traceplot(m9.4)
trankplot(m9.4) 
# clearly a disaster (not mixing, not stationary)

# try again with weakly regularizing priors 
m9.5 <- ulam(
  alist(
    y ~ dnorm(mu , sigma) ,
    mu <- a1 + a2 ,
    a1 ~ dnorm(0 , 10),
    a2 ~ dnorm(0 , 10),
    sigma ~ dexp(1)
  ) ,
  data = list(y = y) ,
  chains = 3
)
precis(m9.5)

traceplot(m9.5) 
trankplot(m9.5) 
# a lot better 




