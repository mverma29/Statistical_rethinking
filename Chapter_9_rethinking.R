# Chp 9 Exercises

# Megan Verma
# 10/11/2022

library(rethinking)

# Markov's island hopping ex:----
num_weeks <- 1e5
positions <- rep(0, num_weeks)
current <- 10
for (i in 1:num_weeks) {
  # record current position
  positions[i] <- current
  # flip coin to generate proposal
  proposal <- current + sample(c(-1, 1) , size = 1)
  # now make sure he loops around the archipelago
  if (proposal < 1)
    proposal <- 10
  if (proposal > 10)
    proposal <- 1
  # move?
  prob_move <- proposal / current
  current <- ifelse(runif(1) < prob_move , proposal , current)
}

plot(1:100 , positions[1:100])
plot(table(positions))

# sample randomly from a high dimension dist (01 dimensions is enough)
# plot radial dist of the points
D <- 10 
T <- 1e3
Y <- rmvnorm(T, rep(0, D), diag(D))
rad_dist <- function(Y)
  sqrt(sum(Y ^ 2))
Rd <- sapply(1:T , function(i)
  rad_dist(Y[i, ]))
dens(Rd)
