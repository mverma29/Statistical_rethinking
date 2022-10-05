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








# multicollinear milk 
