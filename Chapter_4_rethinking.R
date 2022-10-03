# Chapter 4: Geocentric models in R and Stan

# install packages-----

install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# we recommend running this is a fresh R session or restarting your current session
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
install.packages(c("coda","mvtnorm","devtools","loo","dagitty","shape"))
devtools::install_github("rmcelreath/rethinking")

# linear regression---- 

library(rethinking)
data("Howell1")
d <- Howell1
str(d)

# summary function in rethinking pkg is precis
precis(d)

# filter df to individuals 18+ 
d2 <- d[d$age >=18, ]

# plot the dist of heights 
dens(d2$height)

# plot priors
#mu
curve(dnorm(x, 178, 20), from=100, to=250)
# sigma
curve(dunif(x, 0, 50), from=-10, to=60)

# simulate heights by sampling from the prior: prior predictive simulation
sample_mu <- rnorm(1e4, 178, 20)
sample_sigma <- runif(1e4, 0, 50)
prior_h <- rnorm(1e4, sample_mu, sample_sigma)
dens(prior_h)

# prior predictive simulation with large SD
sample_mu <- rnorm(1e4, 178, 100)
prior_h <- rnorm(1e4, sample_mu, sample_sigma)
dens(prior_h)

# grid approx of the posterior dist 
mu.list <- seq( from=150, to=160 , length.out=100 )
sigma.list <- seq( from=7 , to=9 , length.out=100 )
post <- expand.grid( mu=mu.list , sigma=sigma.list )
post$LL <- sapply( 1:nrow(post) , function(i) sum(
  dnorm( d2$height , post$mu[i] , post$sigma[i] , log=TRUE ) ) )
post$prod <- post$LL + dnorm( post$mu , 178 , 20 , TRUE ) +
  dunif( post$sigma , 0 , 50 , TRUE )
post$prob <- exp( post$prod - max(post$prod) )
# contour plot 
contour_xyz( post$mu , post$sigma , post$prob )
# heat map 
image_xyz( post$mu , post$sigma , post$prob )
# sample from the posterior (combinations of the two parameters)
sample.rows <- sample( 1:nrow(post) , size=1e4 , replace=TRUE ,
                       prob=post$prob )
sample.mu <- post$mu[ sample.rows ]
sample.sigma <- post$sigma[ sample.rows ]
# plot samples 
plot( sample.mu , sample.sigma , cex=0.5 , pch=16 , col=col.alpha(rangi2,0.1) )
dens( sample.mu )
dens( sample.sigma )
# summarize the widths of these densities with posterior compatability intervals 
PI( sample.mu )
PI( sample.sigma )

# quadratic approximation of the posterior dist-----
# formula for models
flist <- alist(
  height ~ dnorm( mu , sigma ) ,
  mu ~ dnorm( 178 , 20 ) ,
  sigma ~ dunif( 0 , 50 )
)
# fit the model 
m4.1 <- quap( flist , data=d2 )
# examine posterior dist 
precis(m4.1)

# try again with more informative prior -- SD to 0.1 
m4.2 <- quap(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu ~ dnorm( 178 , 0.1 ) ,
    sigma ~ dunif( 0 , 50 )
  ) , data=d2 )
precis( m4.2 ) # sigma increases a lot to compensate for a very narrow prior
# see the matrix of variances and covariances for model m4.1
vcov(m4.1)
# break it down 
diag( vcov( m4.1 ) )
cov2cor( vcov( m4.1 ) )

# sample from this multi-dimensional posterior
post <- extract.samples( m4.1 , n=1e4 )
head(post)
precis(post)







# linear prediction---- 

# plot height and weight against one another 
plot(d2$height ~ d2$weight)

# plot priors for the simulation of height & weight model
set.seed(2971)
N <- 100
a <- rnorm( N , 178 , 20 )
b <- rnorm( N , 0 , 10 )

plot(
  NULL ,
  xlim = range(d2$weight) ,
  ylim = c(-100, 400) ,
  xlab = "weight" ,
  ylab = "height"
)
abline(h = 0 , lty = 2)
abline(h = 272 , lty = 1 , lwd = 0.5)
mtext("b ~ dnorm(0,10)")
xbar <- mean(d2$weight)
for (i in 1:N)
  curve(
    a[i] + b[i] * (x - xbar) ,
    from = min(d2$weight) ,
    to = max(d2$weight) ,
    add = TRUE ,
    col = col.alpha("black", 0.2)
  )

# restrict to a positive relationship
b <- rlnorm( 1e4 , 0 , 1 )
dens( b , xlim=c(0,5) , adj=0.1 )

set.seed(2971)
N <- 100                   # 100 lines
a <- rnorm( N , 178 , 20 )
b <- rlnorm( N , 0 , 1 )
plot(
  NULL ,
  xlim = range(d2$weight) ,
  ylim = c(-100, 400) ,
  xlab = "weight" ,
  ylab = "height"
)
abline(h = 0 , lty = 2)
abline(h = 272 , lty = 1 , lwd = 0.5)
mtext("b ~ dnorm(0,10)")
xbar <- mean(d2$weight)
for (i in 1:N)
  curve(
    a[i] + b[i] * (x - xbar) ,
    from = min(d2$weight) ,
    to = max(d2$weight) ,
    add = TRUE ,
    col = col.alpha("black", 0.2)
  )

# approx the posterior dist 

# define the average weight, x-bar
xbar <- mean(d2$weight)
# fit model
m4.3 <- quap(alist(
  height ~ dnorm(mu , sigma) ,
  mu <- a + b * (weight - xbar) ,
  a ~ dnorm(178 , 20) ,
  b ~ dlnorm(0 , 1) ,
  sigma ~ dunif(0 , 50)
),
data = d2)

# redo but with logB 
m4.3b <- quap(alist(
  height ~ dnorm(mu , sigma) ,
  mu <- a + exp(log_b) * (weight - xbar),
  a ~ dnorm(178 , 20) ,
  log_b ~ dnorm(0 , 1) ,
  sigma ~ dunif(0 , 50)
),
data = d2)

# describe the quadratic posterior 
precis (m4.3)
round(vcov(m4.3), 3) # very little covariation 

# plot the posterior means over the height and weight data 
plot(height ~ weight , data = d2 , col = rangi2)
post <- extract.samples(m4.3)
a_map <- mean(post$a)
b_map <- mean(post$b)
curve(a_map + b_map * (x - xbar) , add = TRUE)
# extract first 10 cases & re-estimate model 
N <- 10
dN <- d2[1:N ,]
mN <- quap(alist(
  height ~ dnorm(mu , sigma) ,
  mu <- a + b * (weight - mean(weight)) ,
  a ~ dnorm(178 , 20) ,
  b ~ dlnorm(0 , 1) ,
  sigma ~ dunif(0 , 50)
) ,
data = dN)

# plot 20 of the lines:
# extract 20 samples from the posterior
post <- extract.samples(mN , n = 20)
# display raw data and sample size
plot(
  dN$weight ,
  dN$height ,
  xlim = range(d2$weight) ,
  ylim = range(d2$height) ,
  col = rangi2 ,
  xlab = "weight" ,
  ylab = "height"
)
mtext(concat("N = ", N))
# plot the lines, with transparency
for (i in 1:20)
  curve(post$a[i] + post$b[i] * (x - mean(dN$weight)) ,
        col = col.alpha("black", 0.3) ,
        add = TRUE)

# let's look at a single weight value (50 kg)
# make a list of 10,000 values of mu for an individual who weighs 50 kg using samples from post
post <- extract.samples(m4.3)
mu_at_50 <- post$a + post$b * (50 - xbar)
# plot density for this vector
dens(mu_at_50 ,
     col = rangi2 ,
     lwd = 2 ,
     xlab = "mu|weight=50")
# find 89% compatibility interval of mu at 50 kg
PI(mu_at_50 , prob = 0.89)
# repeat for every value of weight in the dataset 
mu <- link(m4.3)
str(mu)
# define sequence of weights to compute predictions for these values will be on the hor axis
weight.seq <- seq(from = 25 , to = 70 , by = 1)
# use link to compute mu for each sample from posterior and for each weight in weight.seq
mu <- link(m4.3 , data = data.frame(weight = weight.seq))
str(mu)

# let's plot dist of mu values at each height
# use type="n" to hide raw data
plot(height ~ weight , d2 , type = "n")
# loop over samples and plot each mu value
for (i in 1:100)
  points(weight.seq , mu[i, ] , pch = 16 , col = col.alpha(rangi2, 0.1))

# summarize the distribution of mu
mu.mean <- apply(mu , 2 , mean)
mu.PI <- apply(mu , 2 , PI , prob = 0.89)

# plot the summaries on top of the data 
# plot raw data, fading out points to make line and interval more visible
plot(height ~ weight , data = d2 , col = col.alpha(rangi2, 0.5))
# plot the MAP line, aka the mean mu for each weight
lines(weight.seq , mu.mean)
# plot a shaded region for 89% PI
shade(mu.PI , weight.seq)

# express uncertainty in both mu and posterior 
sim.height <- sim(m4.3 , data = list(weight = weight.seq))
str(sim.height)
# summarize simulated heights
height.PI <- apply(sim.height , 2 , PI , prob = 0.89)

# plot everything: avg line, shaded region of 89% plausible mu, boundaries of simulated heights
# plot raw data
plot(height ~ weight , d2 , col = col.alpha(rangi2, 0.5))
# draw MAP line
lines(weight.seq , mu.mean)
# draw PI region for line
shade(mu.PI , weight.seq)
# draw PI region for simulated heights
shade(height.PI , weight.seq)

# curves from lines------

# polynomials -----


plot(d$height ~ d$weight)

# build the square of weight.s as a separate variable

d$weight_s <- (d$weight - mean(d$weight)) / sd(d$weight)
d$weight_s2 <- d$weight_s ^ 2

# quadratic model
m4.5   <- quap(
  alist(
    height ~ dnorm(mu , sigma) ,
    mu <- a + b1 * weight_s + b2 * weight_s2 ,
    a ~ dnorm(178 , 20) ,
    b1 ~ dlnorm(0 , 1) ,
    b2 ~ dnorm(0 , 1) ,
    sigma ~ dunif(0 , 50)
  ) ,
  data = d
)

# model summary 
precis(m4.5)

# must plot to understand 

weight.seq <- seq(from       = -2.2 ,
                  to         = 2 ,
                  length.out = 30) 
pred_dat <- list(weight_s = weight.seq , weight_s2 = weight.seq ^ 2)
mu <- link(m4.5 , data = pred_dat)
mu.mean <- apply(mu , 2 , mean)
mu.PI <- apply(mu , 2 , PI , prob = 0.89)
sim.height <- sim(m4.5 , data = pred_dat)
height.PI <- apply(sim.height , 2 , PI , prob = 0.89)

plot(height ~ weight_s , d , col = col.alpha(rangi2, 0.5))
lines(weight.seq , mu.mean)
shade(mu.PI , weight.seq)
shade(height.PI , weight.seq)

# higher order -- cubic regression on weight

d$weight_s3 <- d$weight_s ^ 3

m4.6   <- quap(
  alist(
    height ~ dnorm(mu , sigma) ,
    mu <- a + b1 * weight_s + b2 * weight_s2 + b3 * weight_s3 ,
    a ~ dnorm(178 , 20) ,
    b1 ~ dlnorm(0 , 1) ,
    b2 ~ dnorm(0 , 10) ,
    b3 ~ dnorm(0 , 10) ,
    sigma ~ dunif(0 , 50)
  ) ,
  data = d
)

# plot cubic linear function 
weight.seq <- seq(from       = -2.2 ,
                  to         = 2 ,
                  length.out = 30) 
pred_dat <- list(weight_s = weight.seq , weight_s2 = weight.seq ^ 2,  weight_s3 = weight.seq ^ 3)
mu <- link(m4.6 , data = pred_dat)
mu.mean <- apply(mu , 2 , mean)
mu.PI <- apply(mu , 2 , PI , prob = 0.89)
sim.height <- sim(m4.6 , data = pred_dat)
height.PI <- apply(sim.height , 2 , PI , prob = 0.89)

plot(height ~ weight_s , d , col = col.alpha(rangi2, 0.5))
lines(weight.seq , mu.mean)
shade(mu.PI , weight.seq)
shade(height.PI , weight.seq)

# B-splines-----

library(rethinking)
data(cherry_blossoms)
d <- cherry_blossoms
precis(d)

# plot temp against year 
plot(d$temp ~ d$year)

# define a list of knot positions in "year" variable
d2 <- d[complete.cases(d$temp) ,] # complete cases on temp
num_knots <- 15
knot_list <-
  quantile(d2$year , probs = seq(0, 1, length.out = num_knots))

# build a basis funciton for any list of knots & degree you choose
# 3rd degreee function used here 

library(splines)
B <- bs(d2$year,
        knots     = knot_list[-c(1, num_knots)] ,
        degree    = 3 ,
        intercept = TRUE)

# plot each column against year to display basis functions

plot(
  NULL ,
  xlim = range(d2$year) ,
  ylim = c(0, 1) ,
  xlab = "year" ,
  ylab = "basis value"
)
for (i in 1:ncol(B))
  lines(d2$year , B[, i])

# build model 
m4.7   <- quap(
  alist(
    T ~ dnorm(mu , sigma) ,
    mu <- a + B %*% w , # matrix multiplication
    a ~ dnorm(6, 10),
    w ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ),
  data = list(T = d2$temp , B = B) ,
  start = list(w = rep(0 , ncol(B)))
)

# plot the posterior predictions 

# weighted basis functions 
post <- extract.samples(m4.7)
w <- apply(post$w , 2 , mean)
plot(
  NULL ,
  xlim = range(d2$year) ,
  ylim = c(-2, 2) ,
  xlab = "year" ,
  ylab = "basis * weight"
)
for (i in 1:ncol(B))
  lines(d2$year , w[i] * B[, i])

# 97% posterior interval 
mu    <- link(m4.7)
mu_PI <- apply(mu, 2, PI, 0.97)
plot(d2$year , d2$temp , col = col.alpha(rangi2, 0.3) , pch = 16)
shade(mu_PI , d2$year , col = col.alpha("black", 0.5))
