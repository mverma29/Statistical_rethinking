# Chp 13 Exercises

# Megan Verma
# 10/19/2022

# multilevel tadpoles example------ 
library(rethinking)
data(reedfrogs)
d <- reedfrogs
str(d)

# make the tank cluster variable
d$tank <- 1:nrow(d)
dat <- list(S    = d$surv,
            N    = d$density,
            tank = d$tank)
# approximate posterior (regular model)
m13.1 <- ulam(
  alist(S ~ dbinom(N , p) ,
        logit(p) <- a[tank] ,
        a[tank] ~ dnorm(0 , 1.5)),
  data    = dat ,
  chains  = 4 ,
  log_lik = TRUE
)
precis(m13.1,depth=2) # 48 diff intercepts, one for each tank 

# multilevel model 
m13.2 <- ulam(
  alist(
    S ~ dbinom(N , p) ,
    logit(p) <- a[tank] ,
    a[tank] ~ dnorm(a_bar , sigma) ,
    a_bar ~ dnorm(0 , 1.5) ,
    sigma ~ dexp(1)
  ),
  data    = dat ,
  chains  = 4 ,
  log_lik = TRUE
)

compare(m13.1 , m13.2)


# extract Stan samples
post <- extract.samples(m13.2)
# compute median intercept for each tank
# also transform to probability with logistic
d$propsurv.est <- logistic(apply(post$a , 2 , mean))
# display raw proportions surviving in each tank
plot(
  d$propsurv ,
  ylim = c(0, 1) ,
  pch  = 16 ,
  xaxt = "n" ,
  xlab = "tank" ,
  ylab = "proportion survival" ,
  col  = rangi2
)
axis(1 , at = c(1, 16, 32, 48) , labels = c(1, 16, 32, 48))
# overlay posterior means
points(d$propsurv.est)
# mark posterior mean probability across tanks
abline(h = mean(inv_logit(post$a_bar)) , lty = 2)
# draw vertical dividers between tank densities
abline(v = 16.5 , lwd = 0.5)
abline(v = 32.5 , lwd = 0.5)
text(8 , 0 , "small tanks")
text(16 + 8 , 0 , "medium tanks")
text(32 + 8 , 0 , "large tanks")


# what does inferred population distribution of survival look like? 
# visualize by sampling from posterior dist 
# show first 100 populations in the posterior
plot(
  NULL ,
  xlim = c(-3, 4) ,
  ylim = c(0, 0.35) ,
  xlab = "log-odds survive" ,
  ylab = "Density"
)
for (i in 1:100)
  curve(dnorm(x, post$a_bar[i], post$sigma[i]) ,
        add = TRUE ,
        col = col.alpha("black", 0.2))
# sample 8000 imaginary tanks from the posterior distribution
sim_tanks <- rnorm(8000 , post$a_bar , post$sigma)
# transform to probability and visualize
dens(inv_logit(sim_tanks) , lwd = 2 , adj = 0.1)


# simulate same scenario but with ponds so prediction is goal 
a_bar  <- 1.5
sigma  <- 1.5
nponds <- 60
Ni     <- as.integer(rep(c(5, 10, 25, 35) , each = 15))

set.seed(5005)
a_pond <- rnorm(nponds , mean = a_bar , sd = sigma)

# put simulation in a df to organize it 
dsim <- data.frame(pond   = 1:nponds ,
                   Ni     = Ni ,
                   true_a = a_pond)

# generate a simulated survival count for each pond 
# probability is defined by the logistic function 
dsim$Si <-
  rbinom(nponds , prob = logistic(dsim$true_a) , size = dsim$Ni)

# compute the no-pooling estimates 
# empirical proportions of survivors from each pond
dsim$p_nopool <- dsim$Si / dsim$Ni

# compute partial pooling estimates 
# basic varying intercept model 
dat <- list(Si   = dsim$Si ,
            Ni   = dsim$Ni ,
            pond = dsim$pond)

m13.3 <- ulam(
  alist(
    Si ~ dbinom(Ni , p),
    logit(p) <- a_pond[pond],
    a_pond[pond] ~ dnorm(a_bar , sigma),
    a_bar ~ dnorm(0 , 1.5),
    sigma ~ dexp(1)
  ),
  data   = dat ,
  chains = 4
)

precis(m13.3 , depth = 2)

# compute the predicted survival proportions and add to simulation df
post            <- extract.samples(m13.3)
dsim$p_partpool <- apply(inv_logit(post$a_pond), 2, mean)

# compute true per-pond survival probabilities used to generate the data 
dsim$p_true     <- inv_logit(dsim$true_a)

# compute absolute error between estimates and true varying effects
nopool_error    <- abs(dsim$p_nopool - dsim$p_true)
partpool_error  <- abs(dsim$p_partpool - dsim$p_true)

# plot the results 
plot(
  1:60 ,
  nopool_error ,
  xlab = "pond" ,
  ylab = "absolute error" ,
  col  = rangi2 ,
  pch  = 16
)
points(1:60 , partpool_error)

# calculate the average error rates 
nopool_avg   <- aggregate(nopool_error, list(dsim$Ni), mean)
partpool_avg <- aggregate(partpool_error, list(dsim$Ni), mean)

# multilevel chimpanzees example: more than one cluster-----

library(rethinking)
data(chimpanzees)
d           <- chimpanzees
d$treatment <- 1 + d$prosoc_left + 2 * d$condition
dat_list    <- list(
  pulled_left = d$pulled_left,
  actor       = d$actor,
  block_id    = d$block,
  treatment   = as.integer(d$treatment)
)

set.seed(13)

m13.4 <- ulam(
  alist(
    pulled_left ~ dbinom(1 , p) ,
    logit(p) <- a[actor] + g[block_id] + b[treatment] ,
    b[treatment] ~ dnorm(0 , 0.5),
    ## adaptive priors
    a[actor] ~ dnorm(a_bar , sigma_a),
    g[block_id] ~ dnorm(0 , sigma_g),
    ## hyper-priors
    a_bar ~ dnorm(0 , 1.5), # only one global mean parameter! 
    sigma_a ~ dexp(1),
    sigma_g ~ dexp(1)
  ) ,
  data = dat_list ,
  chains = 4 ,
  cores = 4 ,
  log_lik = TRUE
)

precis(m13.4 , depth = 2)
plot(precis(m13.4, depth = 2)) # also plot

# compare to model with only varying intercepts on actor (ignoring block)

set.seed(14)
m13.5 <- ulam(
  alist(
    pulled_left ~ dbinom(1 , p) ,
    logit(p) <- a[actor] + b[treatment] ,
    b[treatment] ~ dnorm(0 , 0.5),
    a[actor] ~ dnorm(a_bar , sigma_a),
    a_bar ~ dnorm(0 , 1.5),
    sigma_a ~ dexp(1)
  ) ,
  data    = dat_list ,
  chains  = 4 ,
  cores   = 4 ,
  log_lik = TRUE
)

# compare to the model with both clusters 
compare(m13.4 , m13.5)

# try partial pooling on the treatment effects (controversial)
set.seed(15)
m13.6 <- ulam(
  alist(
    pulled_left ~ dbinom(1 , p) ,
    logit(p) <- a[actor] + g[block_id] + b[treatment] ,
    b[treatment] ~ dnorm(0 , sigma_b),
    a[actor] ~ dnorm(a_bar , sigma_a),
    g[block_id] ~ dnorm(0 , sigma_g),
    a_bar ~ dnorm(0 , 1.5),
    sigma_a ~ dexp(1),
    sigma_g ~ dexp(1),
    sigma_b ~ dexp(1)
  ) ,
  data    = dat_list ,
  chains  = 4 ,
  cores   = 4 ,
  log_lik = TRUE
)
coeftab(m13.4 , m13.6)
compare(m13.4 , m13.6)


# divergent transitions ex-----

m13.7 <- ulam(alist(v ~ normal(0, 3),
                    x ~ normal(0, exp(v))),
              data = list(N = 1) ,
              chains = 4)

precis(m13.7) # horrible n_eff & Rhat

traceplot(m13.7) # also horrible 

# re-parameterization (non-centered)


m13.7nc <- ulam(
  alist(v ~ normal(0, 3),
        z ~ normal(0, 1),
        gq > real[1]:x <<- z * exp(v)),
  data = list(N = 1) ,
  chains = 4
)
precis(m13.7nc)

# divergent transitions: multilevel chimpanzees-----

# change the warm-up phase acceptance rate to 0.99 (default is 0.95)
set.seed(13)
m13.4b <- ulam( m13.4 , chains=4 , cores=4 , control=list(adapt_delta=0.99) )

# check divergence 
divergent(m13.4b) # 2 divergent samples 

precis(m13.4b) # n_eff still way below 2000 used

# non-centered reparameterization of model
set.seed(13)
m13.4nc <- ulam(
  alist(
    pulled_left ~ dbinom(1 , p) ,
    logit(p) <- a_bar + z[actor] * sigma_a + # actor intercepts
      x[block_id] * sigma_g +
      b[treatment] ,
    b[treatment] ~ dnorm(0 , 0.5),
    z[actor] ~ dnorm(0 , 1),
    x[block_id] ~ dnorm(0 , 1),
    a_bar ~ dnorm(0 , 1.5),
    sigma_a ~ dexp(1),
    # block intercepts
    sigma_g ~ dexp(1),
    gq > vector[actor]:a <<- a_bar + z * sigma_a,
    gq > vector[block_id]:g <<- x * sigma_g
  ) ,
  data = dat_list ,
  chains = 4 ,
  cores = 4
)

# examine matching n_eff parameters for the models and plot
precis_c <- precis(m13.4 , depth = 2)
precis_nc <- precis(m13.4nc , depth = 2)

pars <- c(
  paste("a[", 1:7, "]", sep = "") ,
  paste("g[", 1:6, "]", sep = "") ,
  paste("b[", 1:4, "]", sep = "") ,
  "a_bar" ,
  "sigma_a" ,
  "sigma_g"
)

neff_table <-
  cbind(precis_c[pars, "n_eff"] , precis_nc[pars, "n_eff"])

plot(
  neff_table ,
  xlim = range(neff_table) ,
  ylim = range(neff_table) ,
  xlab = "n_eff (centered)" ,
  ylab = "n_eff (non-centered)" ,
  lwd = 2
)
abline(a = 0 , b = 1 , lty = 2)

# multilevel posterior predictions ex: chimpanzees, same cluster-----

# looking at chimp 2,  using link function 
chimp  <- 2
d_pred <- list(
  actor = rep(chimp, 4),
  treatment = 1:4,
  block_id = rep(1, 4)
)

p      <- link(m13.4 , data = d_pred)
p_mu   <- apply(p , 2 , mean)
p_ci   <- apply(p , 2 , PI)


# using model definition
post <- extract.samples(m13.4)
str(post) # a is a matrix of with samples on the rows and actors on columns

# plot density for actor 5 
dens(post$a[, 5])

# build our own link function 
p_link <- function(treatment ,
                   actor = 1 ,
                   block_id = 1) {
  logodds <- with(post ,
                  a[, actor] + g[, block_id] + b[, treatment])
  return(inv_logit(logodds))
}

# compute predictions
p_raw <-
  sapply(1:4 , function(i)
    p_link(i , actor = 2 , block_id = 1))
p_mu <- apply(p_raw , 2 , mean)
p_ci <- apply(p_raw , 2 , PI)

# multilevel posterior predictions ex: chimpanzees, new cluster-----

# need to simulate an imaginary chimp at the average intercept 
# new link function, with a twist 
p_link_abar <- function(treatment) {
  logodds <-
    with(post , a_bar + b[, treatment]) # ignores block bc generalizing to new block
  return(inv_logit(logodds))
}

post <- extract.samples(m13.4)
p_raw <- sapply(1:4 , function(i)
  p_link_abar(i))
p_mu <- apply(p_raw , 2 , mean)
p_ci <- apply(p_raw , 2 , PI)
plot(
  NULL ,
  xlab = "treatment" ,
  ylab = "proportion pulled left" ,
  ylim = c(0, 1) ,
  xaxt = "n" ,
  xlim = c(1, 4)
)
axis(1 ,
     at = 1:4 ,
     labels = c("R/N", "L/N", "R/P", "L/P"))
lines(1:4 , p_mu)
shade(p_ci , 1:4)

# sample some random chimpanzees 
a_sim <- with(post , rnorm(length(post$a_bar) , a_bar , sigma_a))

# link function that references those sampled chimps
p_link_asim <- function(treatment) {
  logodds <- with(post , a_sim + b[, treatment])
  return(inv_logit(logodds))
}
p_raw_asim <- sapply(1:4 , function(i)
  p_link_asim(i))

post <- extract.samples(m13.4)
p_raw <- sapply(1:4 , function(i)
  p_link_asim(i))
p_mu <- apply(p_raw , 2 , mean)
p_ci <- apply(p_raw , 2 , PI)
plot(
  NULL ,
  xlab = "treatment" ,
  ylab = "proportion pulled left" ,
  ylim = c(0, 1) ,
  xaxt = "n" ,
  xlim = c(1, 4)
)
axis(1 ,
     at = 1:4 ,
     labels = c("R/N", "L/N", "R/P", "L/P"))
lines(1:4 , p_mu)
shade(p_ci , 1:4)

# show both the variation between actors and summarize 
# plot new "average actors" for each of the treatments

plot(
  NULL ,
  xlab = "treatment" ,
  ylab = "proportion pulled left" ,
  ylim = c(0, 1) ,
  xaxt = "n" ,
  xlim = c(1, 4)
)
axis(1 ,
     at = 1:4 ,
     labels = c("R/N", "L/N", "R/P", "L/P"))
for (i in 1:100)
  lines(1:4 , p_raw_asim[i, ] , col = grau(0.25) , lwd = 2)
