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
xtabs(~ treatment + prosoc_left + condition , d)

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
m11.2 <- quap(
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a + b[treatment] ,
    a ~ dnorm( 0 , 1.5 ),
    b[treatment] ~ dnorm( 0 , 10 )
  ) , data=d )

# sample priors
set.seed(1999)
prior <- extract.prior( m11.2 , n=1e4 )
# prob for pulling left for each treatment 
p <- sapply( 1:4 , function(k) inv_logit( prior$a + prior$b[,k] ) )
# plot priors-- looking at absolute prior difference between first two treatments
dens( abs( p[,1] - p[,2] ) , adj=0.1 )

# try again, prior with SD of 0.5 for beta
m11.3 <- quap(
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a + b[treatment] ,
    a ~ dnorm( 0 , 1.5 ),
    b[treatment] ~ dnorm( 0 , 0.5 )
  ) , data=d )
set.seed(1999)
prior <- extract.prior( m11.3 , n=1e4 )
p <- sapply( 1:4 , function(k) inv_logit( prior$a + prior$b[,k] ) )
# what is average prior difference?
mean( abs( p[,1] - p[,2] ) )
# 0.098 --> 10% avg diff 

# run model using MCMC 

# prior trimmed data list
dat_list <- list(
  pulled_left = d$pulled_left,
  actor = d$actor,
  treatment = as.integer(d$treatment) )

m11.4 <- ulam(
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a[actor] + b[treatment] ,
    a[actor] ~ dnorm( 0 , 1.5 ),
    b[treatment] ~ dnorm( 0 , 0.5 )
  ) , data=dat_list , chains=4 , log_lik=TRUE )
precis( m11.4 , depth=2 )

# look at the intercepts on the outcome scale 
# intercepts are unique to each chimp 
# represent the tendency of each individual to pull the left lever 
post <- extract.samples(m11.4)
p_left <- inv_logit( post$a )
plot( precis( as.data.frame(p_left) ) , xlim=c(0,1) )

# consider the treatment effects 
labs <- c("R/N","L/N","R/P","L/P") # R/N means food on right, no partner
plot( precis( m11.4 , depth=2 , pars="b" ) , labels=labs )

# calc diff between no-partner/partner 
# contrasts on the log odds scale 
diffs <- list(
  db13 = post$b[,1] - post$b[,3],
  db24 = post$b[,2] - post$b[,4] )
plot( precis(diffs) )
# no compelling evidence 

# posterior prediction check: summ the proportions of left pulls for each actor 
# in each treatment and then plot against the posterior predictions 

# calc the proportion in each combo of actor and treatment 
pl <- by( d$pulled_left , list( d$actor , d$treatment ) , mean )
# first actor results: pulling left for each treatment 
pl[1,]

# plot these values against the posterior predictions 
plot( NULL , xlim=c(1,28) , ylim=c(0,1) , xlab="" ,
      ylab="proportion left lever" , xaxt="n" , yaxt="n" )
axis( 2 , at=c(0,0.5,1) , labels=c(0,0.5,1) )
abline( h=0.5 , lty=2 )
for ( j in 1:7 ) abline( v=(j-1)*4+4.5 , lwd=0.5 )
for ( j in 1:7 ) text( (j-1)*4+2.5 , 1.1 , concat("actor ",j) , xpd=TRUE )
for ( j in (1:7)[-2] ) {
  lines( (j-1)*4+c(1,3) , pl[j,c(1,3)] , lwd=2 , col=rangi2 )
  lines( (j-1)*4+c(2,4) , pl[j,c(2,4)] , lwd=2 , col=rangi2 )
}
points( 1:28 , t(pl) , pch=16 , col="white" , cex=1.7 )
points( 1:28 , t(pl) , pch=c(1,1,16,16) , col=rangi2 , lwd=2 )
yoff <- 0.01
text( 1 , pl[1,1]-yoff , "R/N" , pos=1 , cex=0.8 )
text( 2 , pl[1,2]+yoff , "L/N" , pos=3 , cex=0.8 )
text( 3 , pl[1,3]-yoff , "R/P" , pos=1 , cex=0.8 )
text( 4 , pl[1,4]+yoff , "L/P" , pos=3 , cex=0.8 )
mtext( "observed proportions\n" )


# posterior predictions 
dat <- list( actor=rep(1:7,each=4) , treatment=rep(1:4,times=7) )
p_post <- link( m11.4 , data=dat )
p_mu <- apply( p_post , 2 , mean )
p_ci <- apply( p_post , 2 , PI )


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
  cond = d$cond )
m11.5 <- ulam(
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a[actor] + bs[side] + bc[cond] ,
    a[actor] ~ dnorm( 0 , 1.5 ),
    bs[side] ~ dnorm( 0 , 0.5 ),
    bc[cond] ~ dnorm( 0 , 0.5 )
  ) , data=dat_list2 , chains=4 , log_lik=TRUE )

compare( m11.5 , m11.4 , func=PSIS )
