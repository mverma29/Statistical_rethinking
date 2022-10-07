# Chapter 6: The haunted DAG and the causal terror 

# Megan Verma 
# 10/6/2022

library(rethinking)

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





# post-treatment effect: plant fungus example----

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

# DAG of these relationships
library(dagitty)
plant_dag <- dagitty("dag {
    H_0 -> H_1
    F -> H_1
    T -> F
}")

coordinates(plant_dag) <- list(x = c(
  H_0 = 0,
  T = 2,
  F = 1.5,
  H_1 = 1
),
y = c(
  H_0 = 0,
  T = 0,
  F = 0,
  H_1 = 0
))

drawdag(plant_dag)

impliedConditionalIndependencies(plant_dag)

# unobserved variable (moisture) example, influences final plant height & fungus
set.seed(71)
N <- 1000
h0 <- rnorm(N, 10, 2)
treatment <- rep(0:1 , each = N / 2)
M <- rbern(N)
fungus <- rbinom(N , size = 1 , prob = 0.5 - treatment * 0.4 + 0.4 * M)
h1 <- h0 + rnorm(N , 5 + 3 * M)
d2 <-
  data.frame(
    h0 = h0 ,
    h1 = h1 ,
    treatment = treatment ,
    fungus = fungus)

# re-run m6.7 and 6.8 now
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
  data = d2
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
  data = d2
)
precis(m6.8) # impact of treatment is now positive


# collider bias ex: happiness and age, conditioning on marriage----

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

# unknown collider bias ex: grandparents, parents, children-----

N <- 200 # number of grandparent-parent-child triads
b_GP <- 1 # direct effect of G on P
b_GC <- 0 # direct effect of G on C
b_PC <- 1 # direct effect of P on C
b_U <- 2 # direct effect of U on P and C

# draw random observations 
set.seed(1)
U <- 2*rbern( N , 0.5 ) - 1
G <- rnorm( N )
P <- rnorm( N , b_GP*G + b_U*U )
C <- rnorm( N , b_PC*P + b_GC*G + b_U*U )
d <- data.frame( C=C , P=P , G=G , U=U )

# model, controlling for parents
m6.11 <- quap(alist(
  C ~ dnorm(mu , sigma),
  mu <- a + b_PC * P + b_GC * G,
  a ~ dnorm(0 , 1),
  c(b_PC, b_GC) ~ dnorm(0 , 1),
  sigma ~ dexp(1)
),
data = d)
precis(m6.11) # effect of parents is too big (2x as big as it should)

# how do we fix this?
# must measure and condition on unknown (u)

m6.12 <- quap(
  alist(
    C ~ dnorm(mu , sigma),
    mu <- a + b_PC * P + b_GC * G + b_U * U,
    a ~ dnorm(0 , 1),
    c(b_PC, b_GC, b_U) ~ dnorm(0 , 1),
    sigma ~ dexp(1)
  ),
  data = d
)
precis(m6.12)

# DAG showing what to condition on 
library(dagitty)
dag_6.1 <- dagitty("dag {
U [unobserved]
X -> Y
X <- U <- A -> C -> Y
U -> B <- C
}")
adjustmentSets(dag_6.1 , exposure = "X" , outcome = "Y")
# condition on either A or C 


# waffle house divorce extended DAG 
dag_6.2 <- dagitty("dag {
A -> D
A -> M -> D
A <- S -> M
S -> W -> D
}")
adjustmentSets(dag_6.2 , exposure = "W" , outcome = "D")
# what are the conditional independencies?
impliedConditionalIndependencies(dag_6.2)

# DAG question 1 medium: showing what to condition on 
dag_6.3 <- dagitty("dag {
U [unobserved]
X -> Y
X <- U <- A -> C -> Y
U -> B <- C
C <- V -> Y
}")
adjustmentSets(dag_6.3 , exposure = "X" , outcome = "Y")
# condition on either C and V or A 


