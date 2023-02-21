
# Generating 95%CIs for PAFs using Monte Carlo simulation

# Set working directory
setwd("~/.")

# Load packages
pacman::p_load(dplyr, readr, mc2d, skimr)

# 10000 simulations
B <- 1e4

## opioid use data ---

d <- read.csv("simulation_input_git1.csv") # opioid inputs

# rate ratios
d$logrr <- log(d$rr)
d$se <- (d$logrr - log(d$rr_lci)) / qnorm(0.975)
rrs <- mapply(rnorm, mean = d$logrr, sd = d$se, n = list(B), SIMPLIFY = T)
rrs <- t(exp(rrs))

# populations - using modified PERT distribution
pops <- mapply(rpert, min = d$pop_lower, mode = d$pop_mean, max = d$pop_upper, n = list(B), SIMPLIFY = T)
pops <- t(pops)

# calculate PAFs
pafs <- pops * (rrs-1) / (pops * (rrs-1)+1)
paf_results <- t(apply(pafs, 1, quantile, probs = c(0.5, 0.025, 0.975))) * 100
paf_results <- as.data.frame(paf_results)

## all-cause IH data ---

d2 <- read.csv("simulation_input_git2.csv") # all-cause IH inputs

# rate ratios
d2$logrr <- log(d2$rr)
d2$se <- (d2$logrr - log(d2$rr_lci)) / qnorm(0.975)
rrs2 <- mapply(rnorm, mean = d2$logrr, sd = d2$se, n = list(B), SIMPLIFY = T)
rrs2 <- t(exp(rrs2))

# populations 
d2$logp <- log(d2$p)
d2$pse <- (d2$logp - log(d2$p_lci)) / qnorm(0.975)
pops2 <- mapply(rnorm, mean = d2$logp, sd = d2$pse, n = list(B), SIMPLIFY = T)
pops2 <- t(exp(pops2))

# calculate pafs2
pafs2 <- pops2 * (rrs2-1) / (pops2 * (rrs2-1)+1)
paf_results2 <- t(apply(pafs2, 1, quantile, probs = c(0.5, 0.025, 0.975))) * 100
paf_results2 <- as.data.frame(paf_results2)

