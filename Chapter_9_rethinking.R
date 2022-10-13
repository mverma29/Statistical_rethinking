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








# week 4 HW: 
# Q1------

library(rethinking)
d <- sim_happiness( seed=1977 , N_years=1000 )
precis(d)

# rescale age so that the range 18-65 is one unit 
d2 <- d[d$age > 17 ,] # only adults
d2$A <- (d2$age - 18) / (65 - 18) # ranges from 0-1, 0=18 1=65

d2$mid <- d2$married + 1 # mid= married index variable (1=single, 2=married)

m6.9 <- quap(alist(
  happiness ~ dnorm(mu , sigma),
  mu <- a[mid] + bA*A,
  a[mid] ~ dnorm(0 , 1), # 95% of mass in the -2 to +2 interval of happiness
  bA ~ dnorm(0 , 2), #95% of plausible slopes are less than maximally strong
  sigma ~ dexp(1)
) ,
data = d2)

precis(m6.9, depth = 2)

# try a model that omits marriage status 
m6.10 <- quap(alist(
  happiness ~ dnorm(mu , sigma),
  mu <- a + bA * A,
  a ~ dnorm(0 , 1),
  bA ~ dnorm(0 , 2),
  sigma ~ dexp(1)
) ,
data = d2)
precis(m6.10) # no association between age and happiness

# PSIS and WAIC of m6.9 and m6.10
compare(m6.9 , m6.10)
compare(m6.9 , m6.10, func = PSIS)
# according to these criteria, m6.9 is expected to make better predictions
# would interpret as negative (non-causal) effect of age on happiness


# Q2-----

library(rethinking)
data(foxes)
d <- foxes

# standardize variables
d$F <- scale(d$avgfood)
d$G <- scale(d$groupsize)
d$A <- scale(d$area)
d$W <- scale(d$weight)

summary(d)

# different models to predict weight 

# m1: F on W
m1<- quap(
  alist(
    W ~ dnorm(mu , sigma) ,
    mu <- a + bF*F,
    a ~ dnorm(0 , 0.2) ,
    bF ~ dnorm(0 ,0.5) ,
    sigma ~ dexp(1)
  ),
  data = d
)
precis(m1) # very small effect(-0.2)

# m2: G on W 
m2<- quap(
  alist(
    W ~ dnorm(mu , sigma) ,
    mu <- a + bG*G,
    a ~ dnorm(0 , 0.2) ,
    bG ~ dnorm(0 ,0.5) ,
    sigma ~ dexp(1)
  ),
  data = d
)
precis(m2) # 


# m3: A on W 
m3<- quap(
  alist(
    W ~ dnorm(mu , sigma) ,
    mu <- a + bA*A,
    a ~ dnorm(0 , 0.2) ,
    bA ~ dnorm(0 ,0.5) ,
    sigma ~ dexp(1)
  ),
  data = d
)
precis(m3) #


# m4: F and G on W 
m4<- quap(
  alist(
    W ~ dnorm(mu , sigma) ,
    mu <- a + bF*F + bG*G,
    a ~ dnorm(0 , 0.2) ,
    bF ~ dnorm(0 ,0.5) ,
    bG ~ dnorm(0 ,0.5) ,
    sigma ~ dexp(1)
  ),
  data = d
)
precis(m4) 


# m5: F, G and A on W 
m5<- quap(
  alist(
    W ~ dnorm(mu , sigma) ,
    mu <- a + bF*F + bA*A + bG*G,
    a ~ dnorm(0 , 0.2) ,
    bF ~ dnorm(0 ,0.5) ,
    bA ~ dnorm(0 ,0.5) ,
    bG ~ dnorm(0 ,0.5) ,
    sigma ~ dexp(1)
  ),
  data = d
)
precis(m5)


# m6: A and G on W 
m6<- quap(
  alist(
    W ~ dnorm(mu , sigma) ,
    mu <- a + bA*A + bG*G,
    a ~ dnorm(0 , 0.2) ,
    bA ~ dnorm(0 ,0.5) ,
    bG ~ dnorm(0 ,0.5) ,
    sigma ~ dexp(1)
  ),
  data = d
)
precis(m6)

compare(m1, m2, m3, m4, m5, m6)
compare(m1, m2, m3, m4, m5, m6, func=PSIS)
# m5 is best scoring model: F and G on W 
# F has a pos relationship with W 
# A has a pos relationship with W 
# G has a negative relationship with W (greater magnitude than F's relationship)

# Q3-----

data("cherry_blossoms")
d <- cherry_blossoms
summary(d)

# predict doy with temp (consider 2 functions)

# standardize variables
d$D <- standardize(d$doy)
d$T <- standardize(d$temp)

# drop cases with NA 
dd <- d[complete.cases(d$D, d$T), ]

summary(dd)

# different models to predict doy 

# m1: intercept only 
m1<- quap(
  alist(
    D ~ dnorm(mu , sigma) ,
    mu <- a,
    a ~ dnorm(0 , 10) ,
    sigma ~ dexp(1)
  ),
  data = list(D=dd$D, T=dd$T)
)
precis(m1)

# m2: T on D
m2<- quap(
  alist(
    D ~ dnorm(mu , sigma) ,
    mu <- a + bT*T,
    a ~ dnorm(0 , 10) ,
    bT ~ dnorm(0 , 10) ,
    sigma ~ dexp(1)
  ),
  data = list(D=dd$D, T=dd$T)
)
precis(m2)

# m3: T^2 on D 
m3<- quap(
  alist(
    D ~ dnorm(mu , sigma) ,
    mu <- a + b1*T + + b2*T,
    a ~ dnorm(0 , 10) ,
    c(b1,b2) ~ dnorm(0 , 10) ,
    sigma ~ dexp(1)
  ),
  data = list(D=dd$D, T=dd$T)
)
precis(m3)

compare(m1, m2, m3)
compare(m1, m2, m3, func=PSIS)
# m3 is better (on basis of WAIC, PSIS)

Tval <- (9- mean(d$temp, na.rm=TRUE))/sd(d$temp, na.rm=TRUE)
D_sim <- sim(m2, data = list(T=Tval))
# put back on natural scale 
doy_sim <- D_sim*sd(d$doy, na.rm = TRUE) + mean(d$doy, na.rm=TRUE)
dens(doy_sim, lwd=4, col=2, xlab="day in year 1st bloom")
abline(v=89, lty=1)
dens(d$doy, add=TRUE, lwd=3)
