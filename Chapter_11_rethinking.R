# Chp 11 Exercises

# Megan Verma
# 10/14/2022


# logistic regression: prosocial chimpanzees ex------

library(rethinking)
data(chimpanzees)
d <- chimpanzees

# outcome is pulled_left (0/1)
# prosoc_left (0/1) and condition (partner's presence, 0/1) are predictor vars

# build an index variable containing values 1-4 for 4 conditions possible:
# prosoc_left=0, condition=0: food items on right, no partner
# prosoc_left=1, condition=0: food items on left, no partner
# prosoc_left=0, condition=1: food items on right, partner
# prosoc_left=1, condition=1: food items on left, partner

d$treatment <- 1 + d$prosoc_left + 2 * d$condition

# verify using cross-tabs
xtabs( ~ treatment + prosoc_left + condition , d)

# prior predictive simulation
# linear model with just alpha (intercept) parameter

# flat prior-- SD=10
m11.1 <- quap(alist(pulled_left ~ dbinom(1 , p) ,
                    logit(p) <- a ,
                    a ~ dnorm(0 , 10)) , data = d)

# sample from prior
set.seed(1999)
prior <- extract.prior(m11.1 , n = 1e4)

# convert parameter to outcome scale! use inv_logit
p <- inv_logit(prior$a)
dens(p , adj = 0.1)

# try again, with SD of 1.5
m11.1b <- quap(alist(pulled_left ~ dbinom(1 , p) ,
                     logit(p) <- a ,
                     a ~ dnorm(0 , 1.5)) , data = d)

# sample from prior
set.seed(1999)
prior2 <- extract.prior(m11.1b , n = 1e4)

# convert parameter to outcome scale! use inv_logit
p2 <- inv_logit(prior2$a)
dens(p , adj = 0.1)
dens(p2 , adj = 0.1, add = TRUE)
# better than before but still not ideal. we'll still use it

# now prior for beta (slope)
# again try flat SD of 10
m11.2 <- quap(alist(
  pulled_left ~ dbinom(1 , p) ,
  logit(p) <- a + b[treatment] ,
  a ~ dnorm(0 , 1.5),
  b[treatment] ~ dnorm(0 , 10)
) ,
data = d)

# sample priors
set.seed(1999)
prior <- extract.prior(m11.2 , n = 1e4)
# prob for pulling left for each treatment
p <- sapply(1:4 , function(k)
  inv_logit(prior$a + prior$b[, k]))
# plot priors-- looking at absolute prior difference between first two treatments
dens(abs(p[, 1] - p[, 2]) , adj = 0.1)

# try again, prior with SD of 0.5 for beta
m11.3 <- quap(alist(
  pulled_left ~ dbinom(1 , p) ,
  logit(p) <- a + b[treatment] ,
  a ~ dnorm(0 , 1.5),
  b[treatment] ~ dnorm(0 , 0.5)
) ,
data = d)
set.seed(1999)
prior <- extract.prior(m11.3 , n = 1e4)
p <- sapply(1:4 , function(k)
  inv_logit(prior$a + prior$b[, k]))
# what is average prior difference?
mean(abs(p[, 1] - p[, 2]))
# 0.098 --> 10% avg diff

# run model using MCMC

# prior trimmed data list
dat_list <- list(
  pulled_left = d$pulled_left,
  actor = d$actor,
  treatment = as.integer(d$treatment)
)

m11.4 <- ulam(
  alist(
    pulled_left ~ dbinom(1 , p) ,
    logit(p) <- a[actor] + b[treatment] ,
    a[actor] ~ dnorm(0 , 1.5),
    b[treatment] ~ dnorm(0 , 0.5)
  ) ,
  data = dat_list ,
  chains = 4 ,
  log_lik = TRUE
)
precis(m11.4 , depth = 2)

# look at the intercepts on the outcome scale
# intercepts are unique to each chimp
# represent the tendency of each individual to pull the left lever
post <- extract.samples(m11.4)
p_left <- inv_logit(post$a)
plot(precis(as.data.frame(p_left)) , xlim = c(0, 1))

# consider the treatment effects
labs <-
  c("R/N", "L/N", "R/P", "L/P") # R/N means food on right, no partner
plot(precis(m11.4 , depth = 2 , pars = "b") , labels = labs)

# calc diff between no-partner/partner
# contrasts on the log odds scale
diffs <- list(db13 = post$b[, 1] - post$b[, 3],
              db24 = post$b[, 2] - post$b[, 4])
plot(precis(diffs))
# no compelling evidence

# posterior prediction check: summ the proportions of left pulls for each actor
# in each treatment and then plot against the posterior predictions

# calc the proportion in each combo of actor and treatment
pl <- by(d$pulled_left , list(d$actor , d$treatment) , mean)
# first actor results: pulling left for each treatment
pl[1, ]

# plot these values against the posterior predictions
plot(
  NULL ,
  xlim = c(1, 28) ,
  ylim = c(0, 1) ,
  xlab = "" ,
  ylab = "proportion left lever" ,
  xaxt = "n" ,
  yaxt = "n"
)
axis(2 , at = c(0, 0.5, 1) , labels = c(0, 0.5, 1))
abline(h = 0.5 , lty = 2)
for (j in 1:7)
  abline(v = (j - 1) * 4 + 4.5 , lwd = 0.5)
for (j in 1:7)
  text((j - 1) * 4 + 2.5 , 1.1 , concat("actor ", j) , xpd = TRUE)
for (j in (1:7)[-2]) {
  lines((j - 1) * 4 + c(1, 3) , pl[j, c(1, 3)] , lwd = 2 , col = rangi2)
  lines((j - 1) * 4 + c(2, 4) , pl[j, c(2, 4)] , lwd = 2 , col = rangi2)
}
points(1:28 ,
       t(pl) ,
       pch = 16 ,
       col = "white" ,
       cex = 1.7)
points(1:28 ,
       t(pl) ,
       pch = c(1, 1, 16, 16) ,
       col = rangi2 ,
       lwd = 2)
yoff <- 0.01
text(1 , pl[1, 1] - yoff , "R/N" , pos = 1 , cex = 0.8)
text(2 , pl[1, 2] + yoff , "L/N" , pos = 3 , cex = 0.8)
text(3 , pl[1, 3] - yoff , "R/P" , pos = 1 , cex = 0.8)
text(4 , pl[1, 4] + yoff , "L/P" , pos = 3 , cex = 0.8)
mtext("observed proportions\n")


# posterior predictions
dat <- list(actor = rep(1:7, each = 4) ,
            treatment = rep(1:4, times = 7))
p_post <- link(m11.4 , data = dat)
p_mu <- apply(p_post , 2 , mean)
p_ci <- apply(p_post , 2 , PI)


# try a model that splits the location of the prosocial option and the
# presence of partner into separate index variables
# driving hypothesis is an interaction between prosocial and partner

# index variables
d$side <- d$prosoc_left + 1 # right 1, left 2
d$cond <- d$condition + 1 # no partner 1, partner 2

# model
dat_list2 <- list(
  pulled_left = d$pulled_left,
  actor = d$actor,
  side = d$side,
  cond = d$cond
)
m11.5 <- ulam(
  alist(
    pulled_left ~ dbinom(1 , p) ,
    logit(p) <- a[actor] + bs[side] + bc[cond] ,
    a[actor] ~ dnorm(0 , 1.5),
    bs[side] ~ dnorm(0 , 0.5),
    bc[cond] ~ dnorm(0 , 0.5)
  ) ,
  data = dat_list2 ,
  chains = 4 ,
  log_lik = TRUE
)

compare(m11.5 , m11.4 , func = PSIS)


# relative effects (change in odds)
# switching from treatment 2 to 4
post <- extract.samples(m11.4)
mean(exp(post$b[, 4] - post$b[, 2]))

# aggregated binomial
data(chimpanzees)
d <- chimpanzees
d$treatment <- 1 + d$prosoc_left + 2 * d$condition
d$side <- d$prosoc_left + 1 # right 1, left 2
d$cond <- d$condition + 1 # no partner 1, partner 2
d_aggregated <- aggregate(
  d$pulled_left ,
  list(
    treatment = d$treatment ,
    actor = d$actor ,
    side = d$side ,
    cond = d$cond
  ) ,
  sum
)
colnames(d_aggregated)[5] <- "left_pulls"

# exact same inference as before, with aggregated binomial model
dat <- with(
  d_aggregated ,
  list(
    left_pulls = left_pulls,
    treatment = treatment,
    actor = actor,
    side = side,
    cond = cond
  )
)

m11.6 <- ulam(
  alist(
    left_pulls ~ dbinom(18 , p) ,
    logit(p) <- a[actor] + b[treatment] ,
    a[actor] ~ dnorm(0 , 1.5) ,
    b[treatment] ~ dnorm(0 , 0.5)
  ) ,
  data = dat ,
  chains = 4 ,
  log_lik = TRUE
)

precis(m11.6, depth = 2)

compare(m11.6 , m11.4 , func = PSIS)
# PSIS/WAIC will be smaller for aggregated data bc multiplicity term

# try with simple example
# deviance of aggregated 6-in-9-2 * dbinom(6, 9, 0.2, log = TRUE)
# deviance of dis-aggregated-2 * sum(dbern(c(1, 1, 1, 1, 1, 1, 0, 0, 0), 0.2, log =
TRUE))



# grad school admissions ex------

library(rethinking)
data(UCBadmit)
d <- UCBadmit
# data is aggregated by department and gender

# make data list for model
dat_list <- list(
  admit = d$admit,
  applications = d$applications,
  gid = ifelse(d$applicant.gender == "male" , 1 , 2)
)

# model
m11.7 <- ulam(
  alist(admit ~ dbinom(applications , p) ,
        logit(p) <- a[gid] ,
        a[gid] ~ dnorm(0 , 1.5)) ,
  data = dat_list ,
  chains = 4,
  cores = 4
)

precis(m11.7 , depth = 2)

# compute the contrast
# logit scale
post <- extract.samples(m11.7)
diff_a <- post$a[, 1] - post$a[, 2]
# outcome scale (relative, ORs)
diff_p <- inv_logit(post$a[, 1]) - inv_logit(post$a[, 2])
precis(list(diff_a = diff_a , diff_p = diff_p))

# check posterior dist & plot them
postcheck(m11.7)
# draw lines connecting points from same dept
for (i in 1:6) {
  x <- 1 + 2 * (i - 1)
  y1 <- d$admit[x] / d$applications[x]
  y2 <- d$admit[x + 1] / d$applications[x + 1]
  lines(c(x, x + 1) , c(y1, y2) , col = rangi2 , lwd = 2)
  text(x + 0.5 ,
       (y1 + y2) / 2 + 0.05 ,
       d$dept[x] ,
       cex = 0.8 ,
       col = rangi2)
}
# bad prediction-- females are only admitted at lower rates in 2 depts

# try again, adding department into model as an index variable (1-6)
# this will still predict for overall diff in gender
# construct the department index
dat_list$dept_id <- rep(1:6, each = 2)

# model
m11.8 <- ulam(
  alist(
    admit ~ dbinom(applications , p) ,
    logit(p) <- a[gid] + delta[dept_id] ,
    a[gid] ~ dnorm(0 , 1.5) ,
    delta[dept_id] ~ dnorm(0 , 1.5)
  ) ,
  data = dat_list ,
  chains = 4 ,
  cores = 4,
  iter = 4000
)

precis(m11.8 , depth = 2)

# calculate contrasts
post <- extract.samples(m11.8)

# absolute
diff_a <- post$a[, 1] - post$a[, 2]
# relative
diff_p <- inv_logit(post$a[, 1]) - inv_logit(post$a[, 2])
precis(list(diff_a = diff_a , diff_p = diff_p))
# if males have it worse, it's only by 2%

# tabulation to show that females and males tend to apply to diff departments
pg <- with(dat_list , sapply(1:6 , function(k)
  applications[dept_id == k] / sum(applications[dept_id == k])))
rownames(pg) <- c("male", "female")
colnames(pg) <- unique(d$dept)
round(pg , 2)

# Poisson model ex: Oceanic tool complexity-----

library(rethinking)
data(Kline)
d <- Kline
d

# total_tools is the outcome variable 

# make data list 
d$P <- scale( log(d$population) )
d$contact_id <- ifelse( d$contact=="high" , 2 , 1 )

# example of a prior with Normal (0,10) for alpha
# lambda has log-normal dist
curve( dlnorm( x , 0 , 10 ) , from=0 , to=100 , n=200 )

# example of a prior with Normal (3, 0.5) for alpha
curve( dlnorm( x , 3 , 0.5 ) , from=0 , to=100 , n=200 )
# prior mean is about 20 

# example of a prior predictive dist with Normal (0, 10) for beta
N <- 100
a <- rnorm( N , 3 , 0.5 )
b <- rnorm( N , 0 , 10 )
plot( NULL , xlim=c(-2,2) , ylim=c(0,100) )
for ( i in 1:N ) curve( exp( a[i] + b[i]*x ) , add=TRUE , col=grau() )

# try something much tighter 
# beta ~ N(0, 0.2)
set.seed(10)
N <- 100
a <- rnorm( N , 3 , 0.5 )
b <- rnorm( N , 0 , 0.2 )
plot( NULL , xlim=c(-2,2) , ylim=c(0,100) )
for ( i in 1:N ) curve( exp( a[i] + b[i]*x ) , add=TRUE , col=grau() )

# let's look at priors of total tools & un-standardized log population 
x_seq <- seq( from=log(100) , to=log(200000) , length.out=100 )
lambda <- sapply( x_seq , function(x) exp( a + b*x ) )
plot( NULL , xlim=range(x_seq) , ylim=c(0,500) , xlab="log population" ,
      ylab="total tools" )
for ( i in 1:N ) lines( x_seq , lambda[i,] , col=grau() , lwd=1.5 )

# let's look at priors of total tools & natural population 
plot( NULL , xlim=range(exp(x_seq)) , ylim=c(0,500) , xlab="population" ,
      ylab="total tools" )
for ( i in 1:N ) lines( exp(x_seq) , lambda[i,] , col=grau() , lwd=1.5 )

# approx posterior dists 
dat <- list(
  T = d$total_tools ,
  P = d$P ,
  cid = d$contact_id )

# intercept only
m11.9 <- ulam(
  alist(
    T ~ dpois( lambda ),
    log(lambda) <- a,
    a ~ dnorm(3,0.5)
  ), data=dat , chains=4 , cores=4, log_lik=TRUE )

# interaction model
m11.10 <- ulam(
  alist(T ~ dpois( lambda ),
        log(lambda) <- a[cid] + b[cid]*P,
        a[cid] ~ dnorm( 3 , 0.5 ),
        b[cid] ~ dnorm( 0 , 0.2 )
  ), data=dat , chains=4 , cores=4, log_lik=TRUE )

compare( m11.9 , m11.10 , func=PSIS )
# highly influential points warning 

# plot highly influential points 
k <- PSIS( m11.10 , pointwise=TRUE )$k
plot( dat$P , dat$T , xlab="log population (std)" , ylab="total tools" ,
      col=rangi2 , pch=ifelse( dat$cid==1 , 1 , 16 ) , lwd=2 ,
      ylim=c(0,75) , cex=1+normalize(k) )

# set up the horizontal axis values to compute predictions at
ns <- 100
P_seq <- seq( from=-1.4 , to=3 , length.out=ns )

# predictions for cid=1 (low contact)
lambda <- link( m11.10 , data=data.frame( P=P_seq , cid=1 ) )
lmu <- apply( lambda , 2 , mean )
lci <- apply( lambda , 2 , PI )
lines( P_seq , lmu , lty=2 , lwd=1.5 )
shade( lci , P_seq , xpd=TRUE )

# predictions for cid=2 (high contact)
lambda <- link( m11.10 , data=data.frame( P=P_seq , cid=2 ) )
lmu <- apply( lambda , 2 , mean )
lci <- apply( lambda , 2 , PI )
lines( P_seq , lmu , lty=1 , lwd=1.5 )
shade( lci , P_seq , xpd=TRUE )

# plot again, but this time on natural scale 
plot( d$population , d$total_tools , xlab="population" , ylab="total tools" ,
      col=rangi2 , pch=ifelse( dat$cid==1 , 1 , 16 ) , lwd=2 ,
      ylim=c(0,75) , cex=1+normalize(k) )
ns <- 100
P_seq <- seq( from=-5 , to=3 , length.out=ns )
# 1.53 is sd of log(population)
# 9 is mean of log(population)
pop_seq <- exp( P_seq*1.53 + 9 )
lambda <- link( m11.10 , data=data.frame( P=P_seq , cid=1 ) )
lmu <- apply( lambda , 2 , mean )
lci <- apply( lambda , 2 , PI )
lines( pop_seq , lmu , lty=2 , lwd=1.5 )
shade( lci , pop_seq , xpd=TRUE )
lambda <- link( m11.10 , data=data.frame( P=P_seq , cid=2 ) )
lmu <- apply( lambda , 2 , mean )
lci <- apply( lambda , 2 , PI )
lines( pop_seq , lmu , lty=1 , lwd=1.5 )
shade( lci , pop_seq , xpd=TRUE )

# Poisson monastery rate example------

# simulate a month with 1.5 manuscripts per day
num_days <- 30
y <- rpois( num_days , 1.5 )

# simulate new monastery rate 
num_weeks <- 4
y_new <- rpois( num_weeks , 0.5*7 )

# df to organize counts and see exposure for each case 
y_all <- c( y , y_new )
exposure <- c( rep(1,30) , rep(7,4) )
monastery <- c( rep(0,30) , rep(1,4) )
d <- data.frame( y=y_all , days=exposure , monastery=monastery )

# compute the offset
d$log_days <- log( d$days )

# fit the model
m11.12 <- quap(
  alist(
    y ~ dpois( lambda ),
    log(lambda) <- log_days + a + b*monastery,
    a ~ dnorm( 0 , 1 ),
    b ~ dnorm( 0 , 1 )
  ), data=d )

# compute the posterior dist 
# don't include offset bc already on daily scale
post <- extract.samples( m11.12 )
lambda_old <- exp( post$a )
lambda_new <- exp( post$a + post$b )

precis( data.frame( lambda_old , lambda_new ) )


# multinomial models ex: career choice (predictors mateched to outcomes)------

# simulate career choices among 500 individuals
N <- 500             # number of individuals
income <- c(1, 2, 5)   # expected income of each career
score <- 0.5 * income  # scores for each career, based on income
# next line converts scores to probabilities
p <- softmax(score[1], score[2], score[3])
# now simulate choice
# outcome career holds event type values, not counts
career <- rep(NA, N)  # empty vector of choices for each individual
# sample chosen career for each individual
set.seed(34302)
for (i in 1:N)
  career[i] <- sample(1:3 , size = 1 , prob = p)

# fit the model (in stan)
code_m11.13 <- "
data{
    int N; // number of individuals
    int K; // number of possible careers
    int career[N]; // outcome
    vector[K] career_income;
}
parameters{
    vector[K-1] a; // intercepts
    real<lower=0> b; // association of income with choice
} model{
    vector[K] p;
    vector[K] s;
    a ~ normal( 0 , 1 );
    b ~ normal( 0 , 0.5 );
    s[1] = a[1] + b*career_income[1];
    s[2] = a[2] + b*career_income[2];
    s[3] = 0; // pivot
    p = softmax( s );
    career ~ categorical( p );
} "

# set up data list and invoke stan
dat_list <- list(
  N = N ,
  K = 3 ,
  career = career ,
  career_income = income
)
m11.13 <-
  stan(
    model_code = code_m11.13 ,
    data = dat_list ,
    chains = 4,
    cores = 4
  )
precis(m11.13 , 2)

# conduct a counterfactual simulation
# extract the samples and make our own-- imagine doubling the income of career 2

post <- extract.samples(m11.13)
# set up logit scores
s1 <- with(post , a[, 1] + b * income[1])
s2_orig <- with(post , a[, 2] + b * income[2])
s2_new <- with(post , a[, 2] + b * income[2] * 2)
# compute probabilities for original and counterfactual
p_orig <- sapply(1:length(post$b) , function(i)
  softmax(c(s1[i], s2_orig[i], 0)))
p_new <- sapply(1:length(post$b) , function(i)
  softmax(c(s1[i], s2_new[i], 0)))
# summarize
p_diff <- p_new[2, ] - p_orig[2, ]
precis(p_diff)

# multinomial models ex: career choice (predictors matched to observations)------

# effect of family income on a person's choice of career
# predictor variabel has the same value for each linear model, for each row of data
# unique parameter multiplying it in each linear model 

N <- 500
# simulate family incomes for each individual
family_income <- runif(N)
# assign a unique coefficient for each type of event
b <- c(-2, 0, 2)
career <- rep(NA, N)  # empty vector of choices for each individual
for (i in 1:N) {
  score <- 0.5 * (1:3) + b * family_income[i]
  p <- softmax(score[1], score[2], score[3])
  career[i] <- sample(1:3 , size = 1 , prob = p)
}


# fit model
code_m11.14 <- "
data{
    int N; // number of observations
    int K; // number of outcome values
    int career[N]; // outcome
    real family_income[N];
}
parameters{
    vector[K-1] a; // intercepts
    vector[K-1] b; // coefficients on family income
} model{
    vector[K] p;
    vector[K] s;
    a ~ normal(0,1.5);
    b ~ normal(0,1);
    for ( i in 1:N ) {for ( j in 1:(K-1) ) s[j] = a[j] + b[j]*family_income[i];
s[K] = 0; // the pivot
p = softmax( s );
career[i] ~ categorical( p );
} }
"

dat_list <-
  list(
    N = N ,
    K = 3 ,
    career = career ,
    family_income = family_income
  )

m11.14 <- stan(
  model_code = code_m11.14 ,
  data = dat_list ,
  chains = 4,
  cores = 4
)
precis(m11.14 , 2)


# UCB admissions data as multinomial in disguise as Poisson -------


library(rethinking)
data(UCBadmit)
d <- UCBadmit

# binomial model of overall admission probability
m_binom <- quap(alist(admit ~ dbinom(applications, p),
                      logit(p) <- a,
                      a ~ dnorm(0 , 1.5)), data = d)

# Poisson model of overall admission rate and rejection rate
# 'reject' is a reserved word in Stan, cannot use as variable name
dat <- list(admit = d$admit , rej = d$reject)
m_pois <- ulam(
  alist(
    admit ~ dpois(lambda1),
    rej ~ dpois(lambda2),
    log(lambda1) <- a1,
    log(lambda2) <- a2,
    c(a1, a2) ~ dnorm(0, 1.5)
  ),
  data = dat ,
  chains = 3 ,
  cores = 3
)

# inferred binomial prob of admission 
inv_logit(coef(m_binom))

# Poisson model implied prob of admission 
k <- coef(m_pois)
a1 <- k['a1']
a2 <- k['a2']
exp(a1)/(exp(a1)+exp(a2))

# censoring and survival example-----

library(rethinking)
data(AustinCats)
d <- AustinCats
d$adopt <- ifelse(d$out_event == "Adoption" , 1L , 0L)
dat <- list(
  days_to_event = as.numeric(d$days_to_event),
  color_id = ifelse(d$color == "Black" , 1L , 2L) ,
  adopted = d$adopt
)

m11.15 <- ulam(
  alist(
    days_to_event | adopted == 1 ~ exponential(lambda),
    days_to_event |
      adopted == 0 ~ custom(exponential_lccdf(!Y | lambda)),
    lambda <- 1.0 / mu,
    log(mu) <- a[color_id],
    a[color_id] ~ normal(0, 1)
  ),
  data = dat ,
  chains = 4 ,
  cores = 4
)
precis(m11.15 , 2)

# avg time to adoption
post <- extract.samples(m11.15)
post$D <- exp(post$a)
precis(post , 2)

log(0.35/(1-0.35))
inv_logit(3.2)
exp(1.7)
