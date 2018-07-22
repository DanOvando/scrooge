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

selectivity <- rstan::extract(fit, c("length_50_sel"), inc_warmup = TRUE, permuted = FALSE) %>%
  drop() %>%
  as_data_frame() %>%
  gather(chain, value) %>%
  group_by(chain) %>%
  mutate(iteration = 1:length(value),
        variable =  "length_50_sel") %>%
  ungroup()

sigma_r <- rstan::extract(fit, c("sigma_r"), inc_warmup = TRUE, permuted = FALSE) %>%
  drop() %>%
  as_data_frame() %>%
  gather(chain, value) %>%
  group_by(chain) %>%
  mutate(iteration = 1:length(value),
         variable =  "sigma_r") %>%
  ungroup()

compare <- selectivity %>%
  bind_rows(sigma_r) %>%
  mutate(warmup = iteration <= warmup) %>%
  filter(warmup == T) %>%
  ggplot(aes(iteration,value, color = chain)) +
  geom_line() +
  facet_grid(variable~chain, scales = "free_y")

compare

