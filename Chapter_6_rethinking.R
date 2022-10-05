# Chapter 6: The haunted DAG and the causal terror 

# Megan Verma 
# 10/6/2022

# multicollinear legs ex----

N <- 100 #number of individuals

set.seed(909)
height <- rnorm(N,10,2) # sim total height of each leg 
leg_prop <- runif(N,0.4,0.5) #as proportion of height
leg_left <- leg_prop*height + # sim left leg as proportion + error
  rnorm( N , 0 , 0.02 )

leg_right <- leg_prop*height + # sim right leg as proportion + error
  rnorm( N , 0 , 0.02 )
d <- data.frame(height,leg_left,leg_right)

# model height as a function of leg lengths (including BOTH legs)
m6.1 <- quap(
  alist(
    height ~ dnorm(mu , sigma) ,
    mu <- a + bl * leg_left + br * leg_right ,
    a ~ dnorm(10 , 100) ,
    bl ~ dnorm(2 , 10) ,
    br ~ dnorm(2 , 10) ,
    sigma ~ dexp(1)
  ),
  data = d
)
precis(m6.1)
plot(precis(m6.1))

# plot post dist for br 
post <- extract.samples(m6.1)
plot(bl ~ br , post , col = col.alpha(rangi2, 0.1) , pch = 16)

# plot post dist for bl + br 
sum_blbr <- post$bl + post$br
dens( sum_blbr , col=rangi2 , lwd=2 , xlab="sum of bl and br" )


# fit the regression with only1 of the leg lengths
m6.2 <- quap(alist(
  height ~ dnorm(mu , sigma) ,
  mu <- a + bl * leg_left,
  a ~ dnorm(10 , 100) ,
  bl ~ dnorm(2 , 10) ,
  sigma ~ dexp(1)
),
data = d)
precis(m6.2)








# multicollinear milk ex-----

data(milk)
d <- milk
d$K <- scale(d$kcal.per.g)
d$F <- scale(d$perc.fat)
d$L <- scale(d$perc.lactose)

# kcal.per.g regressed on perc.fat
m6.3 <- quap(alist(
  K ~ dnorm(mu , sigma) ,
  mu <- a + bF*F ,
  a ~ dnorm(0 , 0.2) ,
  bF ~ dnorm(0 , 0.5) ,
  sigma ~ dexp(1)
) ,
data = d)

# kcal.per.g regressed on perc.lactose
m6.4 <- quap(alist(
  K ~ dnorm(mu , sigma) ,
  mu <- a + bL*L ,
  a ~ dnorm(0 , 0.2) ,
  bL ~ dnorm(0 , 0.5) ,
  sigma ~ dexp(1)
) ,
data = d)

precis(m6.3) # the more fat, the more kcal in milk
precis(m6.4) # the more lactose, the less kcal in milk 

# regress kcal on both! 
m6.5 <- quap(alist(
  K ~ dnorm(mu , sigma) ,
  mu <- a + bF * F + bL * L ,
  a ~ dnorm(0 , 0.2) ,
  bF ~ dnorm(0 , 0.5) ,
  bL ~ dnorm(0 , 0.5) ,
  sigma ~ dexp(1)
),
data = d)

precis(m6.5) # post means of both bF and bL are closer to 0, SDs are 2x as large
# percent lactose & percent fat contain too much of the same info! 

# check correlation using a pairs plot
pairs (~kcal.per.g + perc.fat + perc.lactose, data=d, col=rangi2)
# perc fat & perc lactose are strongly neg correlated with each other! 



