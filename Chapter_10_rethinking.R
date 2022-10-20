# Chp 10 Exercises

# Megan Verma
# 10/13/2022

library(rethinking)

# maximum entropy ex:----

# put each dist of pebbles in a list
p <- list()
p$A <- c(0,0,10,0,0)
p$B <- c(0,1,8,1,0)
p$C <- c(0,2,6,2,0)
p$D <- c(1,2,4,2,1)
p$E <- c(2,2,2,2,2)


# normalize each such that it is a probability dist 
# divide each count of pebbles by the total # of pebbles
p_norm <- lapply(p , function(q)
  q / sum(q))

# compute the information entropy of each
(H <- sapply(p_norm , function(q)
  - sum(ifelse(q == 0, 0, q * log(
    q
  )))))

# count up ways each dist can be realized 
ways <- c(1,90,1260,37800,113400)
# divide that logarithm by 10 (# of pebbles) to get log ways er pebble
logwayspp <- log(ways)/10

# binomial dist max entropy example------

# marbles example:
# verify that none of the other dist are binomial by putting them in a list
# and passing each ot an expected value formula

# build list of the candidate distributions
p <- list()
p[[1]] <- c(1 / 4, 1 / 4, 1 / 4, 1 / 4)
p[[2]] <- c(2 / 6, 1 / 6, 1 / 6, 2 / 6)
p[[3]] <- c(1 / 6, 2 / 6, 2 / 6, 1 / 6)
p[[4]] <- c(1 / 8, 4 / 8, 2 / 8, 1 / 8)
# compute expected value of each
sapply(p , function(p)
  sum(p * c(0, 1, 1, 2)))

# compute entropy of each distribution
sapply(p , function(p)
  - sum(p * log(p)))
# dist A (first one) has max entropy!


# less special example where p=0.7
# binomial dist of this expected value is :
p <- 0.7
( A <- c( (1-p)^2 , p*(1-p) , (1-p)*p , p^2 ) )











