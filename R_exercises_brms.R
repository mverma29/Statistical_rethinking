# Statistical rethinking

# setup ----

packages <- c("ape", "bayesplot", "brms", "broom", "dagitty", "devtools", "flextable",
              "GGally", "ggdag", "ggdark", "ggmcmc", "ggrepel", "ggthemes", "ggtree", 
              "ghibli", "gtools", "loo", "patchwork", "psych", "rcartocolor", "Rcpp", 
              "remotes", "rstan", "StanHeaders", "statebins", "tidybayes", "tidyverse", 
              "viridis", "viridisLite", "wesanderson")

install.packages(packages, dependencies = T)

# packages that are not in CRAN
devtools::install_github("stan-dev/cmdstanr")
devtools::install_github("EdwinTh/dutchmasters")
devtools::install_github("gadenbuie/ggpomological")
devtools::install_github("rmcelreath/rethinking")
 
devtools::install_github("UrbanInstitute/urbnmapr")
remotes::install_github("stan-dev/posterior")

# Session 2.2: Building a model----

# save globe-tossing data in a tibble: 
# library(tidyverse)

(d <- tibble(toss = c("w", "l", "w", "w", "w", "l", "w", "l", "w")))

# Bayesian updating 

# add the cumulative number of trials, n_trials, and the cumulative number of 
# successes, n_successes (i.e., toss == "w"), to the data.

(
  d <-
    d %>% 
    mutate(n_trials  = 1:9,
           n_success = cumsum(toss == "w"))
)

sequence_length <- 50

d %>% 
  expand(nesting(n_trials, toss, n_success), 
         p_water = seq(from = 0, to = 1, length.out = sequence_length)) %>% 
  group_by(p_water) %>% 
  # you can learn more about lagging here: https://www.rdocumentation.org/packages/stats/versions/3.5.1/topics/lag or here: https://dplyr.tidyverse.org/reference/lead-lag.html
  mutate(lagged_n_trials  = lag(n_trials, k = 1),
         lagged_n_success = lag(n_success, k = 1)) %>% 
  ungroup() %>% 
  mutate(prior      = ifelse(n_trials == 1, .5,
                             dbinom(x    = lagged_n_success, 
                                    size = lagged_n_trials, 
                                    prob = p_water)),
         likelihood = dbinom(x    = n_success, 
                             size = n_trials, 
                             prob = p_water),
         strip      = str_c("n = ", n_trials)) %>% 
  # the next three lines allow us to normalize the prior and the likelihood, 
  # putting them both in a probability metric 
  group_by(n_trials) %>% 
  mutate(prior      = prior / sum(prior),
         likelihood = likelihood / sum(likelihood)) %>%   
  
  # plot!
  ggplot(aes(x = p_water)) +
  geom_line(aes(y = prior), 
            linetype = 2) +
  geom_line(aes(y = likelihood)) +
  scale_x_continuous("proportion water", breaks = c(0, .5, 1)) +
  scale_y_continuous("plausibility", breaks = NULL) +
  theme(panel.grid = element_blank()) +
  facet_wrap(~ strip, scales = "free_y")

# dashed curves are normalized prior densities 
# solid curves are normalized likelihoods 


# Session 2.3: Components of the model----

# 3 components:
# (1) a likelihood function: “the number of ways each conjecture could produce an observation,”
# (2) one or more parameters: “the accumulated number of ways each conjecture cold produce the entire data,”
# (3) a prior: “the initial plausibility of each conjectured cause of the data”

# unobserved variables are called parameters

# distribution function assigned to the observed variable is called a likelihood 

# determine the likelihood of 6 out of 9 tosses coming out water 
dbinom(x = 6, size = 9, prob = .5)

# change the values of prob over the parameter space, [0,1]

tibble(prob = seq(from = 0, to = 1, by = .01)) %>% 
  ggplot(aes(x = prob, y = dbinom(x = 6, size = 9, prob = prob))) +
  geom_line() +
  labs(x = "probability",
       y = "binomial likelihood") +
  theme(panel.grid = element_blank())


# Session 2.4.3: Grid approximation----

# compute the posterior over 1,000 evenly-spaced points on the probability space
# here we do 20 points

(
  d <-
    tibble(p_grid = seq(from = 0, to = 1, length.out = 20),      # define grid
           prior  = 1) %>%                                       # define prior
    mutate(likelihood = dbinom(6, size = 9, prob = p_grid)) %>%  # compute likelihood at each value in grid
    mutate(unstd_posterior = likelihood * prior) %>%             # compute product of likelihood and prior
    mutate(posterior = unstd_posterior / sum(unstd_posterior))   # standardize the posterior, so it sums to 1
)


# show how more points gives better curve in figure 
p1 <-
  d %>% 
  ggplot(aes(x = p_grid, y = posterior)) +
  geom_point() +
  geom_line() +
  labs(subtitle = "20 points",
       x = "probability of water",
       y = "posterior probability") +
  theme(panel.grid = element_blank())

p2 <-
  tibble(p_grid = seq(from = 0, to = 1, length.out = 5),
         prior  = 1) %>%
  mutate(likelihood = dbinom(6, size = 9, prob = p_grid)) %>%
  mutate(unstd_posterior = likelihood * prior) %>%
  mutate(posterior = unstd_posterior / sum(unstd_posterior)) %>% 
  
  ggplot(aes(x = p_grid, y = posterior)) +
  geom_point() +
  geom_line() +
  labs(subtitle = "5 points",
       x = "probability of water",
       y = "posterior probability") +
  theme(panel.grid = element_blank())

p1
p2

# Session 2.4.3: Quadratic approximation----

# aka Laplace 
# posterior dist can be usefully approximated using the Gaussian dist
# the log of a Gaussian dist forms a parabola -- hence "quadratic"

library(rethinking)

globe_qa <- quap(
  alist(
    W ~ dbinom(W + L, p),  # binomial likelihood 
    p ~ dunif(0, 1)        # uniform prior
  ), 
  data = list(W = 6, L = 3)
)

# display summary of quadratic approximation 
precis(globe_qa)

# Session 2.4.3: Markov chain Monte Carlo approximation----

# brms pacakge uses version of MCMC to fit Bayesian models 

library(brms)


# refit a model we previously fit, where w=24 and n=36 

b2.1 <-
  brm(data = list(w = 24), 
      family = binomial(link = "identity"),
      w | trials(36) ~ 0 + Intercept,
      prior(beta(1, 1), class = b, lb = 0, ub = 1),
      seed = 2,
      file = "/Users/meganverma/Desktop/Statistical_Rethinking/outputs/b02.01")

print(b2.1)
