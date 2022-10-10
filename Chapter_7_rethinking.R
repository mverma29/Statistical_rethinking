# Chp 7 Exercises

# Megan Verma
# 10/7/2022

library(rethinking)

# overfitting ex: hominim species-----

sppnames <- c(
  "afarensis",
  "africanus",
  "habilis",
  "boisei",
  "rudolfensis",
  "ergaster",
  "sapiens"
)
brainvolcc <- c(438 , 452 , 612, 521, 752, 871, 1350)
masskg <- c(37.0 , 35.5 , 34.5 , 41.5 , 55.5 , 61.0 , 53.5)
d <- data.frame(species = sppnames ,
                brain = brainvolcc ,
                mass = masskg)

# standardize body mass
# rescale brain volume so largest value observed is 1
d$mass_std <- (d$mass - mean(d$mass)) / sd(d$mass)
d$brain_std <- d$brain / max(d$brain)

# linear model: 
m7.1 <- quap(alist(
  brain_std ~ dnorm(mu, exp(log_sigma)),
  mu <- a + b*mass_std, 
  a ~ dnorm(0.5, 1), # centered on mean brain volume, rescaled
  b ~ dnorm(0,10), # very flat, centered on 0 
  log_sigma ~ dnorm(0.1)
), data=d )

# variation of the residuals-- R2

#compute posterior predictive dist for each obs
set.seed(12)
s <- sim( m7.1 )

# subtract each obs from its prediction to get residual 
r <- apply(s,2,mean) - d$brain_std
resid_var <- var2(r)

# need variance of these residuals and outcome variable (avg squared
# SD of mean)
outcome_var <- var2( d$brain_std )

# R2 at mean prediction: 
1 - resid_var/outcome_var #0.477

# write function to make this process repeatable
R2_is_bad <- function(quap_fit) {
  s <- sim(quap_fit , refresh = 0)
  r <- apply(s, 2, mean) - d$brain_std
  1 - var2(r) / var2(d$brain_std)
}


# quadratic model: 

# define B as a vector-- tell quap how long that is using a start list
m7.2 <- quap(
  alist(
    brain_std ~ dnorm(mu , exp(log_sigma)),
    mu <- a + b[1] * mass_std + b[2] * mass_std ^ 2,
    a ~ dnorm(0.5 , 1),
    b ~ dnorm(0 , 10),
    log_sigma ~ dnorm(0 , 1)
  ),
  data = d ,
  start = list(b = rep(0, 2))
)

# 3rd degree polynomial model

m7.3 <- quap(
  alist(
    brain_std ~ dnorm(mu , exp(log_sigma)),
    mu <- a + b[1] * mass_std + b[2] * mass_std ^ 2 +
      b[3] * mass_std ^ 3,
    a ~ dnorm(0.5 , 1),
    b ~ dnorm(0 , 10),
    log_sigma ~ dnorm(0 , 1)
  ),
  data = d ,
  start = list(b = rep(0, 3))
)

# 4th degree polynomial model

m7.4 <- quap(
  alist(
    brain_std ~ dnorm(mu , exp(log_sigma)),
    mu <- a + b[1] * mass_std + b[2] * mass_std ^ 2 +
      b[3] * mass_std ^ 3 + b[4] * mass_std ^ 4,
    a ~ dnorm(0.5 , 1),
    b ~ dnorm(0 , 10),
    log_sigma ~ dnorm(0 , 1)
  ),
  data = d ,
  start = list(b = rep(0, 4))
)

# 5th degree polynomial model

m7.5 <- quap(
  alist(
    brain_std ~ dnorm(mu , exp(log_sigma)),
    mu <- a + b[1] * mass_std + b[2] * mass_std ^ 2 +
      b[3] * mass_std ^ 3 + b[4] * mass_std ^ 4 +
      b[5] * mass_std ^ 5,
    a ~ dnorm(0.5 , 1),
    b ~ dnorm(0 , 10),
    log_sigma ~ dnorm(0 , 1)
  ),
  data = d ,
  start = list(b = rep(0, 5))
)

# 6th degree polynomial model-- replace SD with constant 0.001
m7.6 <- quap(
  alist(
    brain_std ~ dnorm(mu , 0.001),
    mu <- a + b[1] * mass_std + b[2] * mass_std ^ 2 +
      b[3] * mass_std ^ 3 + b[4] * mass_std ^ 4 +
      b[5] * mass_std ^ 5 + b[6] * mass_std ^ 6,
    a ~ dnorm(0.5 , 1),
    b ~ dnorm(0 , 10)
  ),
  data = d ,
  start = list(b = rep(0, 6))
)

# plot all models: 
# extract samples from posterior, compute posterior predictive dist,
# summarize, then plot:

# m7.1
post <- extract.samples(m7.1) 
mass_seq <-
  seq(
    from = min(d$mass_std) ,
    to = max(d$mass_std) ,
    length.out = 100
  )
l <- link(m7.1 , data = list(mass_std = mass_seq))
mu <- apply(l , 2 , mean)
ci <- apply(l , 2 , PI)
plot(brain_std ~ mass_std , data = d)
lines(mass_seq , mu)
shade(ci , mass_seq)

# m7.2
post <- extract.samples(m7.2) 
mass_seq <-
  seq(
    from = min(d$mass_std) ,
    to = max(d$mass_std) ,
    length.out = 100
  )
l <- link(m7.2 , data = list(mass_std = mass_seq))
mu <- apply(l , 2 , mean)
ci <- apply(l , 2 , PI)
plot(brain_std ~ mass_std , data = d)
lines(mass_seq , mu)
shade(ci , mass_seq)

# m7.3
post <- extract.samples(m7.3) 
mass_seq <-
  seq(
    from = min(d$mass_std) ,
    to = max(d$mass_std) ,
    length.out = 100
  )
l <- link(m7.3 , data = list(mass_std = mass_seq))
mu <- apply(l , 2 , mean)
ci <- apply(l , 2 , PI)
plot(brain_std ~ mass_std , data = d)
lines(mass_seq , mu)
shade(ci , mass_seq)

# m7.4
post <- extract.samples(m7.4) 
mass_seq <-
  seq(
    from = min(d$mass_std) ,
    to = max(d$mass_std) ,
    length.out = 100
  )
l <- link(m7.4 , data = list(mass_std = mass_seq))
mu <- apply(l , 2 , mean)
ci <- apply(l , 2 , PI)
plot(brain_std ~ mass_std , data = d)
lines(mass_seq , mu)
shade(ci , mass_seq)

# m7.5
post <- extract.samples(m7.5) 
mass_seq <-
  seq(
    from = min(d$mass_std) ,
    to = max(d$mass_std) ,
    length.out = 100
  )
l <- link(m7.5 , data = list(mass_std = mass_seq))
mu <- apply(l , 2 , mean)
ci <- apply(l , 2 , PI)
plot(brain_std ~ mass_std , data = d)
lines(mass_seq , mu)
shade(ci , mass_seq)

# m7.6
post <- extract.samples(m7.6) 
mass_seq <-
  seq(
    from = min(d$mass_std) ,
    to = max(d$mass_std) ,
    length.out = 100
  )
l <- link(m7.6 , data = list(mass_std = mass_seq))
mu <- apply(l , 2 , mean)
ci <- apply(l , 2 , PI)
plot(brain_std ~ mass_std , data = d)
lines(mass_seq , mu)
shade(ci , mass_seq)
# R2 is 1-- perfect fit 

# log-probability scores for models-----

set.seed(1)
lppd(m7.1 , n = 1e4)

# for all models in this chapter 
set.seed(1)
sapply(list(m7.1, m7.2, m7.3, m7.4, m7.5, m7.6) , function(m)
  sum(lppd(m)))
# more complex models have larger scores! 


# model comparison examples:-----

# from chp 6 post-treatment effect: plant fungus example----

set.seed(71)
# number of plants
N <- 100
# simulate initial heights
h0 <- rnorm(N, 10, 2)
# assign treatments and simulate fungus and growth
treatment <- rep(0:1 , each = N / 2)
fungus <- rbinom(N , size = 1 , prob = 0.5 - treatment * 0.4)
h1 <- h0 + rnorm(N, 5 - 3 * fungus)

# compose a clean data frame
d <-
  data.frame(
    h0 = h0 ,
    h1 = h1 ,
    treatment = treatment ,
    fungus = fungus
  )

precis(d)

# prior dist of proportion parameter: log-normal(0,0.25)
# because we expect the proportion of plant's initial height to be always positive
# (as it is a proportion) and generally larger than 1 
sim_p <- rlnorm(1e4 , 0 , 0.25)
precis(data.frame(sim_p))
# expects ~ 40% shrinkage to 50% growth, mean is 4% growth 

# fit this linear height of plant model (p is, proportion of plant initial height)
m6.6 <- quap(alist(h1 ~ dnorm(mu , sigma),
                   mu <- h0 * p,
                   p ~ dlnorm(0 , 0.25),
                   sigma ~ dexp(1)), data = d)
precis(m6.6) # on average, 40% growth


# fit the model with the treatment & fungus parameters too
# approx the posterior 
m6.7 <- quap(
  alist(
    h1 ~ dnorm(mu , sigma),
    mu <- h0 * p,
    p <- a + bt * treatment + bf * fungus,
    a ~ dlnorm(0 , 0.2) ,
    bt ~ dnorm(0 , 0.5),
    bf ~ dnorm(0 , 0.5),
    sigma ~ dexp(1)
  ),
  data = d
)
precis(m6.7) # no effect of treatment on height

# try again, omitting fungus
m6.8 <- quap(
  alist(
    h1 ~ dnorm(mu , sigma),
    mu <- h0 * p,
    p <- a + bt * treatment,
    a ~ dlnorm(0 , 0.2),
    bt ~ dnorm(0 , 0.5),
    sigma ~ dexp(1)
  ),
  data = d
)
precis(m6.8) # impact of treatment is now positive

set.seed(11)
WAIC(m6.7)

set.seed(77)
compare(m6.6 , m6.7 , m6.8)

# 99% interval (equiv to z-score of 2.6) of the difference between model 6.7 & 6.8 
40.9 + c(-1,1)*10.49*2.6 # large diff in expected out of sample accuracy 

plot(compare(m6.6 , m6.7 , m6.8))

# WAIC comparison says m6.7 is the best predictive model, but we know it's causally wrong

# inspect dSE of 6.6 and 6.8 
set.seed(93)
compare(m6.6 , m6.7 , m6.8)@dSE #4.93
# dSE is larger than the difference in WAIC
# treatment effect has little relative impact on outcome but is causally related


# outliers: from chp 5 divorce/waffle house example----

data(WaffleDivorce)
d <- WaffleDivorce
d$A <- standardize(d$MedianAgeMarriage)
d$D <- standardize(d$Divorce)
d$M <- standardize(d$Marriage)

# effect of age on divorce rate
m5.1 <- quap(alist(
  D ~ dnorm(mu , sigma) ,
  mu <- a + bA * A ,
  a ~ dnorm(0 , 0.2) ,
  bA ~ dnorm(0 , 0.5) ,
  sigma ~ dexp(1)
) ,
data = d)

# effect of marriage rate on divorce rate
m5.2 <- quap(alist(
  D ~ dnorm(mu , sigma) ,
  mu <- a + bM * M ,
  a ~ dnorm(0 , 0.2) ,
  bM ~ dnorm(0 , 0.5) ,
  sigma ~ dexp(1)
) ,
data = d)

# effect of age + marriage rate on divorce rate
m5.3 <- quap(alist(
  D ~ dnorm(mu , sigma) ,
  mu <- a + bM * M + bA * A ,
  a ~ dnorm(0 , 0.2) ,
  bM ~ dnorm(0 , 0.5) ,
  bA ~ dnorm(0 , 0.5) ,
  sigma ~ dexp(1)
) ,
data = d)

precis(m5.1)
precis(m5.2)
precis(m5.3) 
# marriage rate has little association w/ divorce rate once age is included

# compare models using PSIS 
set.seed(24071847)
compare(m5.1 , m5.2 , m5.3 , func = PSIS)

# look at individual states in m5.3
set.seed(24071847)
PSIS_m5.3 <- PSIS(m5.3, pointwise = TRUE)

# plot the individual penalty values from WAIC
# shows relationship btw Pareto k and the information theoretic prediction penalty
set.seed(24071847)
WAIC_m5.3 <- WAIC(m5.3, pointwise = TRUE)
plot(
  PSIS_m5.3$k ,
  WAIC_m5.3$penalty ,
  xlab = "PSIS Pareto k" ,
  ylab = "WAIC penalty" ,
  col = rangi2 ,
  lwd = 2
)

WAIC(m5.3)

# re-estimate the divorce model using the student-T dist with v=2 
# to reduce influence of outliers
m5.3t <- quap(
  alist(
    D ~ dstudent(2 , mu , sigma) ,
    mu <- a + bM * M + bA * A ,
    a ~ dnorm(0 , 0.2) ,
    bM ~ dnorm(0 , 0.5) ,
    bA ~ dnorm(0 , 0.5) ,
    sigma ~ dexp(1)
  ) ,
  data = d
)

# look at PSIS in m5.3t
set.seed(24071847)
PSIS(m5.3t)

precis(m5.3)
precis(m5.3t)
# coefficient bA has gotten farther from 0 in student-T dist because
# outlier (Idaho) caused model to under-estimate association 
