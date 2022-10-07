# Chp 7 Exercises

# Megan Verma
# 10/7/2022

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




