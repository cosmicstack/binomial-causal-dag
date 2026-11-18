library(tidyverse)
library(rethinking)
library(dagitty)
library(ggdag)

# Draw Causal DAG
dag <- dagitty("dag{N -> V; N -> AN -> AC; E -> P -> D; P -> AC; N -> E -> V; V -> S; V -> D -> S; D -> AC -> S; T -> S; U -> S; V -> AC}")
coordinates(dag) <- list(
  x = c(S = 2, T = 2, U = 2, D = 1.375, V = 1, AC = 1, E = 0, P = 0.25, N = 0, AN = 0.5),
  y = c(S = 0, T = 2, U = -2, D = 0, V = 2, AC = -2, E = 2, P = 0, N = -2, AN = -2)
)
ggdag(dag)

impliedConditionalIndependencies(dag)

# Data simulation
set.seed(1994)
n <- 1000
V <- rbeta(n, 3.5, 5) * 0.5
D <- 10*V + rbinom(n, 15, 0.3)
AC <- 0.1*D + 0.5*V + rbeta(n, 0.5, 0.5)
p <- inv_logit(0.3 + 1.1*V - 0.17*D + 1.4*AC + rnorm(n))
N <- rpois(n, 500)
S <- rbinom(n, N, p)

dat_sim <- data.frame(
  V, D, AC, p, N, S
)

# First Bayesian Model
## General starter
dat0 <- list(
  N = N,
  S = S
)

m0 <- ulam(
  alist(
    S ~ dbinom(N, p),
    logit(p) <- a,
    a ~ dnorm(0, 1.5)
  ),
  data = dat0,
  chains = 4,
  cores = 4
)
trankplot(m0)
precis(m0)

# Posterior Analysis
m0.p <- link(m0)
m0.p.mean <- apply(m0.p, 2, mean)

data.frame(p.true = dat_sim$p, p.gen = m0.p.mean) %>%
  pivot_longer(cols = c("p.true", "p.gen"), names_to = "type", values_to = "prob") %>%
  ggplot() +
  geom_density(aes(prob, fill = type), alpha = 0.3)

m0.p.diff <- dat_sim$p - m0.p.mean
ggplot() +
  geom_density(aes(m0.p.diff), fill = "cornsilk") +
  geom_vline(xintercept = 0, linetype = 2, color = "salmon") +
  theme_minimal()

## Total effect of V; no adjustment needed

dat1 <- list(
  N = N,
  S = S,
  V = V
)

m1 <- ulam(
  alist(
    S ~ dbinom(N, p),
    logit(p) <- a + bV*V,
    a ~ dnorm(0, 1.5),
    bV ~ dnorm(0, 0.5)
  ),
  data = dat1,
  chains = 4,
  cores = 4
)
trankplot(m1)
precis(m1)

# Prior Checks
m1.prior <- extract.prior(m1)
dens(m1.prior$a)
dens(m1.prior$bV)

m1.prior.p <- inv_logit(m1.prior$a + m1.prior$bV*V)
dens(m1.prior.p)
m1.prior.S <- rbinom(n, N, m1.prior.p)
hist(S)
hist(m1.prior.S)

# Posterior Analysis
postcheck(m1, window = 50, n=25)
m1.post <- extract.samples(m1, n=1e3)
m1.p <- link(m1)
m1.p.mean <- apply(m1.p, 2, mean)
m1.p.PI <- apply(m1.p, 2, PI)

data.frame(p.true = dat_sim$p, p.gen = m1.p.mean) %>%
  pivot_longer(cols = c("p.true", "p.gen"), names_to = "type", values_to = "prob") %>%
  ggplot() +
  geom_density(aes(prob, fill = type), alpha = 0.3)

m1.p.diff <- dat_sim$p - m1.p.mean
ggplot() +
  geom_density(aes(m1.p.diff), fill = "cornsilk") +
  geom_vline(xintercept = 0, linetype = 2, color = "salmon") +
  theme_minimal()


## Total effect of V; no adjustment needed
## But over-dispersed model

m2 <- ulam(
  alist(
    S ~ dbetabinom(N, pbar, theta),
    logit(pbar) <- a + bV*V,
    a ~ dnorm(0, 1.5),
    bV ~ dnorm(0, 0.5),
    transpars> theta <<- phi + 2.0,
    phi ~ dexp(1)
  ),
  data = dat1,
  chains = 4,
  cores = 4
)
trankplot(m2)
precis(m2)

# Posterior Analysis
postcheck(m2)
m2.post <- extract.samples(m1, n=1e3)
m2.p <- link(m2)
m2.p.mean <- apply(m2.p, 2, mean)
m2.p.PI <- apply(m2.p, 2, PI)

data.frame(p.true = dat_sim$p, p.gen = m2.p.mean) %>%
  pivot_longer(cols = c("p.true", "p.gen"), names_to = "type", values_to = "prob") %>%
  ggplot() +
  geom_density(aes(prob, fill = type), alpha = 0.3)

m2.p.diff <- dat_sim$p - m2.p.mean
ggplot() +
  geom_density(aes(m2.p.diff), fill = "cornsilk") +
  geom_vline(xintercept = 0, linetype = 2, color = "salmon") +
  theme_minimal()

# m2 base, but direct effect of V
dat1.1 <- list(
  N = N,
  S = S,
  V = V,
  AC = AC,
  D = D
)

m2.1 <- ulam(
  alist(
    S ~ dbetabinom(N, pbar, theta),
    logit(pbar) <- a + bV*V + bAC*AC + bD*D,
    a ~ dnorm(0, 1.5),
    c(bV, bAC, bD) ~ dnorm(0, 0.5),
    transpars> theta <<- phi + 2.0,
    phi ~ dexp(1)
  ),
  data = dat1.1,
  chains = 4,
  cores = 4
)
trankplot(m2.1)
precis(m2.1)

# for a std 5 node DAG, with 1 agent node
# AC = 0.1; D = 5
dat_cf.1 <- data.frame(V=c(0.2, 0.5, 0.8), AC = 0.1, D = 5)
m2.1.p <- link(m2.1, data = dat_cf.1)
m2.1.sim <- sim(m2.1, data = dat_cf.1)

m2.1.diff1 <- m2.1.p[, 2] - m2.1.p[, 1]
dens(m2.1.diff1)

m2.1.diff2 <- m2.1.p[, 3] - m2.1.p[, 2]
dens(m2.1.diff2)

m2.1.diff3 <- m2.1.p[, 3] - m2.1.p[, 1]
dens(m2.1.diff3)
