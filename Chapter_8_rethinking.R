# Chp 8 Exercises

# Megan Verma
# 10/11/2022

library(rethinking)

# single interation: Africa ruggedness GDP ex----- 

data(rugged)
d <- rugged

# make log version of outcome
d$log_gdp <- log(d$rgdppc_2000)

# extract countries with GDP data
dd <- d[complete.cases(d$rgdppc_2000) ,]

# rescale variables
dd$log_gdp_std <- dd$log_gdp / mean(dd$log_gdp)
dd$rugged_std <- dd$rugged / max(dd$rugged) 

summary(dd$log_gdp_std) # range: 0.72-1.28, mean=1
summary(dd$rugged_std) # range: 0.0004-0.316, mean=0.215

# model with very broad priors
m8.1 <- quap(alist(
  log_gdp_std ~ dnorm(mu , sigma) ,
  mu <- a + b * (rugged_std - 0.215) ,
  a ~ dnorm(1 , 1) ,
  b ~ dnorm(0 , 1) ,
  sigma ~ dexp(1)
) ,
data = dd)

#prior predictions
set.seed(7)
prior <- extract.prior(m8.1)

# set up the plot dimensions
plot(
  NULL ,
  xlim = c(0, 1) ,
  ylim = c(0.5, 1.5) ,
  xlab = "ruggedness" ,
  ylab = "log GDP"
)
abline(h = min(dd$log_gdp_std) , lty = 2)
abline(h = max(dd$log_gdp_std) , lty = 2)
# draw 50 lines from the prior
rugged_seq <- seq(from = -0.1 ,
                  to = 1.1 ,
                  length.out = 30)
mu <-
  link(m8.1 ,
       post = prior ,
       data = data.frame(rugged_std = rugged_seq))
for (i in 1:50)
  lines(rugged_seq , mu[i, ] , col = col.alpha("black", 0.3))

# how many slopes are extreme under this prior (abs val > 0.6?)
sum(abs(prior$b) > 0.6) / length(prior$b) # more than half 


# try again with tighter priors
# a ~ N(0,0.1)
# b ~ N(0,0.3)

m8.1 <- quap(alist(
  log_gdp_std ~ dnorm(mu , sigma) ,
  mu <- a + b * (rugged_std - 0.215) ,
  a ~ dnorm(1 , 0.1) ,
  b ~ dnorm(0 , 0.3) ,
  sigma ~ dexp(1)
) ,
data = dd)

#extract prior & plot as before

#prior predictions
set.seed(7)
prior <- extract.prior(m8.1)

# set up the plot dimensions
plot(
  NULL ,
  xlim = c(0, 1) ,
  ylim = c(0.5, 1.5) ,
  xlab = "ruggedness" ,
  ylab = "log GDP"
)
abline(h = min(dd$log_gdp_std) , lty = 2)
abline(h = max(dd$log_gdp_std) , lty = 2)
# draw 50 lines from the prior
rugged_seq <- seq(from = -0.1 ,
                  to = 1.1 ,
                  length.out = 30)
mu <-
  link(m8.1 ,
       post = prior ,
       data = data.frame(rugged_std = rugged_seq))
for (i in 1:50)
  lines(rugged_seq , mu[i, ] , col = col.alpha("black", 0.3))
# much better priors, though still some implausibly strong slopes

# posterior dist 
precis(m8.1)
# no association between ruggedness of terrain and log GDP 

# make variable to index Africa (1) or not (2)
dd$cid <- ifelse(dd$cont_africa == 1 , 1 , 2)

# model with index variable for African continent or not
m8.2 <- quap(alist(
  log_gdp_std ~ dnorm(mu , sigma) ,
  mu <- a[cid] + b * (rugged_std - 0.215) ,
  a[cid] ~ dnorm(1 , 0.1) ,
  b ~ dnorm(0 , 0.3) ,
  sigma ~ dexp(1)
) ,
data = dd)

# compare models using WAIC 
compare(m8.1 , m8.2)
# difference is 64, so including continent seems to be picking up 
# important variation in the sample

precis(m8.2 , depth = 2)
# intercept for africa is lower then for other continents 

# posterior contrast between the 2 intercepts
post <- extract.samples(m8.2)
diff_a1_a2 <- post$a[, 1] - post$a[, 2]
PI(diff_a1_a2)
# difference is def below 0

# sample from the posterior and compute the predicted means & 
# intervals for African & non-African nations 
rugged.seq <- seq(from = -0.1 ,
                  to = 1.1 ,
                  length.out = 30)
# compute mu over samples, fixing cid=2
mu.NotAfrica <- link(m8.2 ,
                     data = data.frame(cid = 2 , rugged_std = rugged.seq))
# compute mu over samples, fixing cid=1
mu.Africa <- link(m8.2 ,
                  data = data.frame(cid = 1 , rugged_std = rugged.seq))
# summarize to means and intervals
mu.NotAfrica_mu <- apply(mu.NotAfrica , 2 , mean)
mu.NotAfrica_ci <- apply(mu.NotAfrica , 2 , PI , prob = 0.97)
mu.Africa_mu <- apply(mu.Africa , 2 , mean)
mu.Africa_ci <- apply(mu.Africa , 2 , PI , prob = 0.97)


# include an interaction between ruggedness and being in Africa to get a new slope
m8.3 <- quap(alist(
  log_gdp_std ~ dnorm(mu , sigma) ,
  mu <- a[cid] + b[cid] * (rugged_std - 0.215) ,
  a[cid] ~ dnorm(1 , 0.1) ,
  b[cid] ~ dnorm(0 , 0.3) ,
  sigma ~ dexp(1)
) ,
data = dd)

# inspect marginal posterior dist 
precis(m8.3, depth=2)
# slope is essentially reversed in Africa 

# use PSIS to compare this model to the previous 2
compare(m8.1, m8.2, m8.3, func = PSIS)

# plot PSIS Pareto k values for m8.3
plot( PSIS( m8.3 , pointwise=TRUE )$k )

# plot the interaction
# plot Africa - cid=1
d.A1 <- dd[dd$cid == 1 ,]
plot(
  d.A1$rugged_std ,
  d.A1$log_gdp_std ,
  pch = 16 ,
  col = rangi2 ,
  xlab = "ruggedness (standardized)" ,
  ylab = "log GDP (as proportion of mean)" ,
  xlim = c(0, 1)
)
mu <-
  link(m8.3 , data = data.frame(cid = 1 , rugged_std = rugged_seq))
mu_mean <- apply(mu , 2 , mean)
mu_ci <- apply(mu , 2 , PI , prob = 0.97)
lines(rugged_seq , mu_mean , lwd = 2)
shade(mu_ci , rugged_seq , col = col.alpha(rangi2, 0.3))
mtext("African nations")

# plot non-Africa - cid=2
d.A0 <- dd[dd$cid == 2 ,]
plot(
  d.A0$rugged_std ,
  d.A0$log_gdp_std ,
  pch = 1 ,
  col = "black" ,
  xlab = "ruggedness (standardized)" ,
  ylab = "log GDP (as proportion of mean)" ,
  xlim = c(0, 1)
)
mu <-
  link(m8.3 , data = data.frame(cid = 2 , rugged_std = rugged_seq))
mu_mean <- apply(mu , 2 , mean)
mu_ci <- apply(mu , 2 , PI , prob = 0.97)
lines(rugged_seq , mu_mean , lwd = 2)
shade(mu_ci , rugged_seq)
mtext("Non-African nations")

# plot reverse interpretation: association of being in Africa with
# log GDP depends on terrain ruggedness

# compute diff between a nation in Africa and outside Africa, holding 
# ruggedness constantrugged_seq <- seq(from=-0.2,to=1.2,length.out=30)
muA <- link(m8.3 , data = data.frame(cid = 1, rugged_std = rugged_seq))
muN <- link(m8.3 , data = data.frame(cid = 2, rugged_std = rugged_seq))
delta <- muA - muN

# plot delta (the counterfactual)
plot(
  delta,
  pch = 1 ,
  col = "black" ,
  xlab = "ruggedness" ,
  ylab = "expected diff of log GDP" ,
  xlim = c(0, 1),
  ylim = c(-0.3, 0.2)
)
delta_mean <- apply(delta , 2 , mean)
delta_ci <- apply(delta , 2 , PI , prob = 0.97)
lines(rugged_seq , delta_mean , lwd = 2)
shade(delta_ci , rugged_seq , col = col.alpha(rangi2, 0.3))
mtext("African nations")
abline(h = 0 , lty = 2)


# continuous variable interaction: winter flower ex-----

data(tulips)
d <- tulips
str(d)
# blooms is outcome 
# water & shade are predictor vars 
# water is ordered categorical: low, med, high (1-3)
# shade is light exposure: high, med, low (1-3)
# bed is a cluster of plants from same section of greenhouse 

# center W and S, scale B by its maximum 
d$blooms_std <- d$blooms / max(d$blooms)
d$water_cent <- d$water - mean(d$water)
d$shade_cent <- d$shade - mean(d$shade)


# model w/o interaction
m8.4 <- quap(
  alist(
    blooms_std ~ dnorm(mu , sigma) ,
    mu <- a + bw * water_cent + bs * shade_cent ,
    a ~ dnorm(0.5 , 0.25),
    #95% of obs within 0-1 of scaled blooms var
    bw ~ dnorm(0 , 0.25),
    #95% of prior slopes are from -0.5 to 0.5
    bs ~ dnorm(0 , 0.25),
    #95% of prior slopes are from -0.5 to 0.5
    sigma ~ dexp(1)
  ) ,
  data = d
  
)

# model with interaction
m8.5 <- quap(
  alist(
    blooms_std ~ dnorm(mu , sigma) ,
    mu <-
      a + bw * water_cent + bs * shade_cent + bws * water_cent * shade_cent ,
    a ~ dnorm(0.5 , 0.25) ,
    bw ~ dnorm(0 , 0.25) ,
    bs ~ dnorm(0 , 0.25) ,
    bws ~ dnorm(0 , 0.25) ,
    sigma ~ dexp(1)
  ) ,
  data = d
)


# draw posterior predictions for m8.4
par(mfrow = c(1, 3)) # 3 plots in 1 row
for (s in-1:1) {
  idx <- which(d$shade_cent == s)
  plot(
    d$water_cent[idx] ,
    d$blooms_std[idx] ,
    xlim = c(-1, 1) ,
    ylim = c(0, 1) ,
    xlab = "water" ,
    ylab = "blooms" ,
    pch = 16 ,
    col = rangi2
  )
  mu <-
    link(m8.4 , data = data.frame(shade_cent = s , water_cent = -1:1))
  for (i in 1:20)
    lines(-1:1 , mu[i, ] , col = col.alpha("black", 0.3))
}

# prior predictions
set.seed(7) 
prior <- extract.prior(m8.5)
