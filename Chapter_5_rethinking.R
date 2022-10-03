# Chapter 5: The many variables & the spurious waffles

# Megan Verma 
# 10/3/2022

# marriage rate data-----

# load data and copy
library(rethinking)
data(WaffleDivorce)
d   <- WaffleDivorce

# standardize variables
d$A <- scale(d$MedianAgeMarriage)
d$D <- scale(d$Divorce)