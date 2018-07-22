library(tidyverse)
library(rstan)

reprex_data <- readRDS("reprex_data.RDS")

chains <- 2

refresh <- 25

iter <- 6000

warmup <-  4000

seed <- 42

cores <- 2

adapt_delta <- 0.8

max_treedepth <- 12

scrooge_file <- "reprex"

inits <-
  map(
    1:chains,
    ~ list(
      log_initial_effort = jitter(log(reprex_data$observed_effort[1])),
      sigma_r = jitter(1e-3)
    )
  )

fit <-
  rstan::stan(
    file = here::here("src", paste0(scrooge_file, ".stan")),
    data = reprex_data,
    chains = chains,
    refresh = refresh,
    cores = cores,
    iter = iter,
    warmup = warmup,
    control = list(adapt_delta = adapt_delta,
                   max_treedepth = max_treedepth),
    init = inits,
    seed = seed
  )


