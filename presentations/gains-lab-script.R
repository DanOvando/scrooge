# library(rstan)
# library(spasm)
# library(FishLife)
# library(ggridges)
# library(scales)
# library(patchwork)
# library(LBSPR)
library(ggridges)
library(wesanderson)
library(scales)
# library(tidyverse)

rstan::rstan_options(auto_write = TRUE)


# Recruitment and Fishing Pulses

set.seed(42)
fish <-
  create_fish(
    query_fishlife = T,
    mat_mode = "length",
    linf = 100,
    max_age = 30,
    time_step = 1,
    sigma_r = 0.4,
    rec_ac = 0.4,
    price = 10,
    price_cv = 0,
    price_ac = 0,
    steepness = 0.6
  )



fleet <- create_fleet(
  fish = fish,
  cost = 3,
  cost_cv =  0,
  cost_ac = 0,
  q_cv = 0,
  q_ac = 0,
  profit_lags = 4,
  length_50_sel = 50,
  theta = 0.1,
  fleet_model = "constant-effort",
  initial_effort = 100
)

fleet$e_msy <-  NA

fleet$p_msy <-  NA

sim <- sim_fishery(
  fish = fish,
  fleet = fleet,
  manager = create_manager(mpa_size = 0),
  num_patches = 1,
  sim_years = 100,
  burn_year = 50,
  time_step = fish$time_step,
  est_msy = F,
  tune_costs = F,
  b_v_bmsy_oa = 1.5
)


length_comps <- sim %>%
  select(year, patch, age, numbers, numbers_caught) %>%
  nest(-year,-patch, .key = n_at_age) %>%
  mutate(catch_length = map(
    n_at_age,
    ~ spasm::sample_lengths(
      n_at_age = .x,
      cv = fish$cv_len,
      k = fish$vbk,
      linf = fish$linf,
      t0 = fish$t0,
      sample_type = "catch",
      percent_sampled = 1,
      time_step = fish$time_step,
      linf_buffer = 1.25
    )
  )) %>%
  select(year, catch_length) %>%
  unnest() %>%
  mutate(numbers = round(numbers)) %>%
  mutate(year = year - min(year) + 1) %>%
  group_by(year) %>%
  mutate(scaled_numbers = numbers / sum(numbers))

r_lengths <- length_comps %>%
  filter(year > 75) %>%
  ggplot(aes(length_bin, year, height = scaled_numbers, group = year)) +
  geom_density_ridges(stat = "identity") +
  labs(x = "Length (cm)", title = "Proportional Length Distribution")


# fishing pulses

set.seed(123)
fish <-
  create_fish(
    query_fishlife = T,
    mat_mode = "length",
    linf = 100,
    max_age = 30,
    time_step = 1,
    sigma_r = 0,
    rec_ac = 0,
    price = 10,
    price_cv = 0,
    price_ac = 0,
    steepness = 0.6
  )



fleet <- create_fleet(
  fish = fish,
  cost = 3,
  cost_cv =  0,
  cost_ac = 0,
  q_cv = 0,
  q_ac = 0,
  profit_lags = 4,
  length_50_sel = 50,
  theta = 0.1,
  fleet_model = "constant-effort",
  sigma_effort = 0.2,
  effort_ac = 0.25,
  initial_effort = 25
)

fleet$e_msy <-  NA

fleet$p_msy <-  NA

sim <- sim_fishery(
  fish = fish,
  fleet = fleet,
  manager = create_manager(mpa_size = 0),
  num_patches = 1,
  sim_years = 100,
  burn_year = 50,
  time_step = fish$time_step,
  est_msy = F,
  tune_costs = F,
  b_v_bmsy_oa = 1.5
)


length_comps <- sim %>%
  select(year, patch, age, numbers, numbers_caught) %>%
  nest(-year,-patch, .key = n_at_age) %>%
  mutate(catch_length = map(
    n_at_age,
    ~ spasm::sample_lengths(
      n_at_age = .x,
      cv = fish$cv_len,
      k = fish$vbk,
      linf = fish$linf,
      t0 = fish$t0,
      sample_type = "catch",
      percent_sampled = 1,
      time_step = fish$time_step,
      linf_buffer = 1.25
    )
  )) %>%
  select(year, catch_length) %>%
  unnest() %>%
  mutate(numbers = round(numbers)) %>%
  mutate(year = year - min(year) + 1) %>%
  group_by(year) %>%
  mutate(scaled_numbers = numbers / sum(numbers))

sim %>%
  group_by(year) %>%
  summarise(effort = mean(effort)) %>%
  ggplot(aes(year, effort)) +
  geom_point()

f_lengths <- length_comps %>%
  filter(year > 75) %>%
  ggplot(aes(length_bin, year, height = scaled_numbers, group = year)) +
  geom_density_ridges(stat = "identity") +
  labs(x = "Length (cm)", title = "Proportional Length Distribution")



# Ideal conditions

easy <- fisheries_sandbox %>%
  filter(fleet_model == "open-access",
         sigma_r == min(sigma_r),
         sigma_effort == min(sigma_effort),
         price_cv == min(price_cv),
         cost_cv == min(cost_cv),
         steepness == min(steepness),
         obs_error == min(obs_error),
         b_v_bmsy_oa == min(b_v_bmsy_oa)) %>%
  slice(1)

easy <- easy %>%
  mutate(prepped_fishery = map(prepped_fishery, subsample_data, window = 10, period = "end")) %>%
  mutate(scrooge_fit = pmap(
    list(
      data = map(prepped_fishery, "scrooge_data"),
      fish = map(prepped_fishery, "fish"),
      fleet = map(prepped_fishery, "fleet")),
    fit_scrooge,
    iter = 4000,
    warmup = 2000,
    adapt_delta = 0.8,
    economic_model = 1,
    effort_data_weight = 1,
    scrooge_file = "scrooge",
    in_clouds = F,
    experiment = "pfo",
    max_f_v_fmsy_increase = 0.5,
    chains = 1,
    cv_effort = 0.5,
    max_expansion = 1.5
  ))

# easy <- easy %>%
#   mutate(lime_fit = pmap(list(
#     data = map(prepped_fishery, "scrooge_data"),
#     fish = map(prepped_fishery, "fish"),
#     fleet = map(prepped_fishery, "fleet")
#   ), fit_lime))
#
# easy <- easy %>%
#   mutate(processed_lime = map2(
#     lime_fit,
#     map(prepped_fishery, "sampled_years"),
#     safely(process_lime)
#   ))


easy <- easy %>%
  mutate(processed_scrooge = map2(
    scrooge_fit,
    map(prepped_fishery, "sampled_years"),
    process_scrooge
  )) %>%
  mutate(observed = map(prepped_fishery, "simed_fishery")) %>%
  mutate(scrooge_performance = pmap(
    list(observed = observed,
         predicted = processed_scrooge),
    judge_performance
  )) %>%
  mutate(
    scrooge_rec_performance = pmap(
      list(observed = observed,
           predicted = processed_scrooge),
      judge_performance,
      observed_variable = rec_dev,
      predicted_variable = "log_rec_dev_t"
    )
  )  %>%
  mutate(lcomps = map2(prepped_fishery, processed_scrooge, ~process_lcomps(.x, .y$n_tl))) %>%
  mutate(
    rmse = map_dbl(scrooge_performance, ~ .x$comparison_summary$rmse),
    bias =  map_dbl(scrooge_performance, ~ .x$comparison_summary$bias)
  ) %>%
  arrange(rmse)

# assumptions stressed

medium <- fisheries_sandbox %>%
  filter(fleet_model == "open-access",
         sigma_r == max(sigma_r),
         sigma_effort == max(sigma_effort),
         price_cv == max(price_cv),
         cost_cv == max(cost_cv),
         steepness == min(steepness),
         obs_error == min(obs_error),
         b_v_bmsy_oa == min(b_v_bmsy_oa)) %>%
  slice(1)

medium <- medium %>%
  mutate(prepped_fishery = map(prepped_fishery, subsample_data, window = 10, period = "end")) %>%
  mutate(scrooge_fit = pmap(
    list(
      data = map(prepped_fishery, "scrooge_data"),
      fish = map(prepped_fishery, "fish"),
      fleet = map(prepped_fishery, "fleet")),
    fit_scrooge,
    iter = 4000,
    warmup = 2000,
    adapt_delta = 0.8,
    economic_model = 1,
    effort_data_weight = 1,
    scrooge_file = "scrooge",
    in_clouds = F,
    experiment = "pfo",
    max_f_v_fmsy_increase = 0.5,
    chains = 1,
    cv_effort = 0.5,
    max_expansion = 1.5
  ))

# medium <- medium %>%
#   mutate(lime_fit = pmap(list(
#     data = map(prepped_fishery, "scrooge_data"),
#     fish = map(prepped_fishery, "fish"),
#     fleet = map(prepped_fishery, "fleet")
#   ), fit_lime))
#
# medium <- medium %>%
#   mutate(processed_lime = map2(
#     lime_fit,
#     map(prepped_fishery, "sampled_years"),
#     safely(process_lime)
#   ))


medium <- medium %>%
  mutate(processed_scrooge = map2(
    scrooge_fit,
    map(prepped_fishery, "sampled_years"),
    process_scrooge
  )) %>%
  mutate(observed = map(prepped_fishery, "simed_fishery")) %>%
  mutate(scrooge_performance = pmap(
    list(observed = observed,
         predicted = processed_scrooge),
    judge_performance
  )) %>%
  mutate(
    scrooge_rec_performance = pmap(
      list(observed = observed,
           predicted = processed_scrooge),
      judge_performance,
      observed_variable = rec_dev,
      predicted_variable = "log_rec_dev_t"
    )
  )  %>%
  mutate(lcomps = map2(prepped_fishery, processed_scrooge, ~process_lcomps(.x, .y$n_tl))) %>%
  mutate(
    rmse = map_dbl(scrooge_performance, ~ .x$comparison_summary$rmse),
    bias =  map_dbl(scrooge_performance, ~ .x$comparison_summary$bias)
  ) %>%
  arrange(rmse)


# Totally Wrong


toughest <- fisheries_sandbox %>%
  filter(fleet_model == "supplied-catch",
         sigma_r == max(sigma_r),
         sigma_effort == max(sigma_effort),
         price_cv == max(price_cv),
         cost_cv == max(cost_cv),
         steepness == min(steepness),
         obs_error == min(obs_error),
         b_v_bmsy_oa == min(b_v_bmsy_oa)) %>%
  slice(1)

toughest <- toughest %>%
  mutate(prepped_fishery = map(prepped_fishery, subsample_data, window = 10, period = "end")) %>%
  mutate(scrooge_fit = pmap(
    list(
      data = map(prepped_fishery, "scrooge_data"),
      fish = map(prepped_fishery, "fish"),
      fleet = map(prepped_fishery, "fleet")),
    fit_scrooge,
    iter = 4000,
    warmup = 2000,
    adapt_delta = 0.8,
    economic_model = 1,
    effort_data_weight = 1,
    scrooge_file = "scrooge",
    in_clouds = F,
    experiment = "pfo",
    max_f_v_fmsy_increase = 0.5,
    chains = 1,
    cv_effort = 0.5,
    max_expansion = 1.5
  ))

# toughest <- toughest %>%
#   mutate(lime_fit = pmap(list(
#     data = map(prepped_fishery, "scrooge_data"),
#     fish = map(prepped_fishery, "fish"),
#     fleet = map(prepped_fishery, "fleet")
#   ), fit_lime))
#
# toughest <- toughest %>%
#   mutate(processed_lime = map2(
#     lime_fit,
#     map(prepped_fishery, "sampled_years"),
#     safely(process_lime)
#   ))


toughest <- toughest %>%
  mutate(processed_scrooge = map2(
    scrooge_fit,
    map(prepped_fishery, "sampled_years"),
    process_scrooge
  )) %>%
  mutate(observed = map(prepped_fishery, "simed_fishery")) %>%
  mutate(scrooge_performance = pmap(
    list(observed = observed,
         predicted = processed_scrooge),
    judge_performance
  )) %>%
  mutate(
    scrooge_rec_performance = pmap(
      list(observed = observed,
           predicted = processed_scrooge),
      judge_performance,
      observed_variable = rec_dev,
      predicted_variable = "log_rec_dev_t"
    )
  )  %>%
  mutate(lcomps = map2(prepped_fishery, processed_scrooge, ~process_lcomps(.x, .y$n_tl))) %>%
  mutate(
    rmse = map_dbl(scrooge_performance, ~ .x$comparison_summary$rmse),
    bias =  map_dbl(scrooge_performance, ~ .x$comparison_summary$bias)
  ) %>%
  arrange(rmse)



# how much does it help ---------------------------------------------------

helps <- fisheries_sandbox %>%
  filter(fleet_model == "open-access",
         sigma_r == min(sigma_r),
         sigma_effort == min(sigma_effort),
         price_cv == max(price_cv),
         cost_cv == max(cost_cv),
         steepness == min(steepness),
         obs_error == min(obs_error),
         b_v_bmsy_oa == min(b_v_bmsy_oa),
         q_cv == min(q_cv),
         q_ac == min(q_ac)) %>%
  slice(1)

helps$summary_plot

econ_helps <- helps %>%
  mutate(prepped_fishery = map(prepped_fishery, subsample_data, window = 10, period = "end")) %>%
  mutate(scrooge_fit = pmap(
    list(
      data = map(prepped_fishery, "scrooge_data"),
      fish = map(prepped_fishery, "fish"),
      fleet = map(prepped_fishery, "fleet")),
    fit_scrooge,
    iter = 4000,
    warmup = 2000,
    adapt_delta = 0.8,
    economic_model = 1,
    effort_data_weight = 1,
    scrooge_file = "scrooge",
    in_clouds = F,
    experiment = "pfo",
    max_f_v_fmsy_increase = 0.1,
    chains = 1,
    cv_effort = 0.1,
    max_expansion = 1.5
  ))

econ_helps <- econ_helps %>%
  mutate(processed_scrooge = map2(
    scrooge_fit,
    map(prepped_fishery, "sampled_years"),
    process_scrooge
  )) %>%
  mutate(observed = map(prepped_fishery, "simed_fishery")) %>%
  mutate(scrooge_performance = pmap(
    list(observed = observed,
         predicted = processed_scrooge),
    judge_performance
  )) %>%
  mutate(
    scrooge_rec_performance = pmap(
      list(observed = observed,
           predicted = processed_scrooge),
      judge_performance,
      observed_variable = rec_dev,
      predicted_variable = "log_rec_dev_t"
    )
  )  %>%
  mutate(lcomps = map2(prepped_fishery, processed_scrooge, ~process_lcomps(.x, .y$n_tl))) %>%
  mutate(
    rmse = map_dbl(scrooge_performance, ~ .x$comparison_summary$rmse),
    bias =  map_dbl(scrooge_performance, ~ .x$comparison_summary$bias)
  ) %>%
  arrange(rmse)

noecon_helps <- helps %>%
  mutate(prepped_fishery = map(prepped_fishery, subsample_data, window = 10, period = "end")) %>%
  mutate(scrooge_fit = pmap(
    list(
      data = map(prepped_fishery, "scrooge_data"),
      fish = map(prepped_fishery, "fish"),
      fleet = map(prepped_fishery, "fleet")),
    fit_scrooge,
    iter = 4000,
    warmup = 2000,
    adapt_delta = 0.8,
    economic_model = 0,
    effort_data_weight = 0,
    scrooge_file = "scrooge",
    in_clouds = F,
    experiment = "pfo",
    max_f_v_fmsy_increase = 0.1,
    chains = 1,
    cv_effort = 0.1,
    max_expansion = 1.5
  ))

noecon_helps <- noecon_helps %>%
  mutate(processed_scrooge = map2(
    scrooge_fit,
    map(prepped_fishery, "sampled_years"),
    process_scrooge
  )) %>%
  mutate(observed = map(prepped_fishery, "simed_fishery")) %>%
  mutate(scrooge_performance = pmap(
    list(observed = observed,
         predicted = processed_scrooge),
    judge_performance
  )) %>%
  mutate(
    scrooge_rec_performance = pmap(
      list(observed = observed,
           predicted = processed_scrooge),
      judge_performance,
      observed_variable = rec_dev,
      predicted_variable = "log_rec_dev_t"
    )
  )  %>%
  mutate(lcomps = map2(prepped_fishery, processed_scrooge, ~process_lcomps(.x, .y$n_tl))) %>%
  mutate(
    rmse = map_dbl(scrooge_performance, ~ .x$comparison_summary$rmse),
    bias =  map_dbl(scrooge_performance, ~ .x$comparison_summary$bias)
  ) %>%
  arrange(rmse)

true_f <- helps$prepped_fishery[[1]]$simed_fishery %>%
  group_by(year) %>%
  summarise(f = mean(f)) %>%
  mutate(year = year - min(year) + 1)


compare_econ <- econ_helps$processed_scrooge[[1]]$f_t %>%
  mutate(model = "With Economic Data") %>%
  bind_rows( noecon_helps$processed_scrooge[[1]]$f_t %>% mutate(model = "W/O Economic Data")) %>%
  rename(f_hat = value) %>%
  left_join(true_f, by = "year") %>%
  ungroup()

trend_plot <- compare_econ %>%
  group_by(year, model) %>%
  mutate(rank_f = percent_rank(f_hat)) %>%
  summarise(mean_f_hat = median(f_hat),
            upper_75 = max(f_hat[rank_f <= 0.95]),
            lower_75 = min(f_hat[rank_f >= 0.05]),
            f = mean(f)) %>%
  ungroup() %>%
  ggplot() +
  geom_ribbon(aes(
    x = year,
    ymax = upper_75,
    ymin = lower_75,
    fill = model
  ), alpha = 0.25,
  color = "grey",
  show.legend = F) +
  geom_line(aes(year, mean_f_hat, color = model), linetype = 2, size = 1.25) +
  geom_line(aes(year, f), size = 1.25) +
  geom_point(aes(year, f), size = 4) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(y = "Fishing Mortality",
       caption = "Ribbon is 90% Credible Interval") +
  scale_color_manual(values = wes_palette("Royal1")) +
  scale_fill_manual(values = wes_palette("Royal1"))


bias_plot <- compare_econ %>%
  group_by(iteration, model) %>%
  summarise(rmse = sqrt(mean((f_hat - f)^2)),
            bias = median((f_hat - f) / f)) %>%
  ungroup() %>%
  group_by(model) %>%
  mutate(mean_bias = mean(bias),
         mean_rmse = mean(rmse)) %>%
  ggplot() +
  geom_vline(aes(xintercept = 0), size = 2) +
  geom_density(aes(bias, fill = model), alpha = 0.5, color = NA) +
  geom_vline(aes(xintercept = mean_bias, color = model)) +
  scale_x_continuous(labels = percent) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title = "Bias") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = wes_palette("Royal1")) +
  scale_fill_manual(values = wes_palette("Royal1"))

rmse_plot <- compare_econ %>%
  group_by(iteration, model) %>%
  summarise(rmse = sqrt(mean((f_hat - f)^2)),
            bias = median((f_hat - f) / f)) %>%
  ungroup() %>%
  group_by(model) %>%
  mutate(mean_bias = mean(bias),
         mean_rmse = mean(rmse)) %>%
  ggplot() +
  geom_density(aes(rmse, fill = model), alpha = 0.5, color = NA, show.legend = F) +
  geom_vline(aes(xintercept = mean_rmse, color = model), show.legend = F) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title = "RMSE") +
  scale_color_manual(values = wes_palette("Royal1")) +
  scale_fill_manual(values = wes_palette("Royal1"))

helps_plot <- ((rmse_plot / bias_plot) | trend_plot) + plot_layout(ncol = 2, widths = c(1, 2))


# hurts -------------------------------------------------------------------

hurts <- fisheries_sandbox %>%
  filter(fleet_model == "constant-effort",
         sigma_r == min(sigma_r),
         sigma_effort == min(sigma_effort),
         price_cv == min(price_cv),
         cost_cv == min(cost_cv),
         steepness == max(steepness),
         obs_error == min(obs_error),
         b_v_bmsy_oa == min(b_v_bmsy_oa),
         q_cv == min(q_cv),
         q_ac == min(q_ac)) %>%
  slice(1)

econ_hurts <- hurts %>%
  mutate(prepped_fishery = map(prepped_fishery, subsample_data, window = 10, period = "end")) %>%
  mutate(scrooge_fit = pmap(
    list(
      data = map(prepped_fishery, "scrooge_data"),
      fish = map(prepped_fishery, "fish"),
      fleet = map(prepped_fishery, "fleet")),
    fit_scrooge,
    iter = 4000,
    warmup = 2000,
    adapt_delta = 0.8,
    economic_model = 1,
    effort_data_weight = 1,
    scrooge_file = "scrooge",
    in_clouds = F,
    experiment = "pfo",
    max_f_v_fmsy_increase = 0.5,
    chains = 1,
    cv_effort = 0.1,
    max_expansion = 1.5
  ))

econ_hurts <- econ_hurts %>%
  mutate(processed_scrooge = map2(
    scrooge_fit,
    map(prepped_fishery, "sampled_years"),
    process_scrooge
  )) %>%
  mutate(observed = map(prepped_fishery, "simed_fishery")) %>%
  mutate(scrooge_performance = pmap(
    list(observed = observed,
         predicted = processed_scrooge),
    judge_performance
  )) %>%
  mutate(
    scrooge_rec_performance = pmap(
      list(observed = observed,
           predicted = processed_scrooge),
      judge_performance,
      observed_variable = rec_dev,
      predicted_variable = "log_rec_dev_t"
    )
  )  %>%
  mutate(lcomps = map2(prepped_fishery, processed_scrooge, ~process_lcomps(.x, .y$n_tl))) %>%
  mutate(
    rmse = map_dbl(scrooge_performance, ~ .x$comparison_summary$rmse),
    bias =  map_dbl(scrooge_performance, ~ .x$comparison_summary$bias)
  ) %>%
  arrange(rmse)

noecon_hurts <- hurts %>%
  mutate(prepped_fishery = map(prepped_fishery, subsample_data, window = 10, period = "end")) %>%
  mutate(scrooge_fit = pmap(
    list(
      data = map(prepped_fishery, "scrooge_data"),
      fish = map(prepped_fishery, "fish"),
      fleet = map(prepped_fishery, "fleet")),
    fit_scrooge,
    iter = 4000,
    warmup = 2000,
    adapt_delta = 0.8,
    economic_model = 0,
    effort_data_weight = 0,
    scrooge_file = "scrooge",
    in_clouds = F,
    experiment = "pfo",
    max_f_v_fmsy_increase = 0.5,
    chains = 1,
    cv_effort = 0.5,
    max_expansion = 1.5
  ))

noecon_hurts <- noecon_hurts %>%
  mutate(processed_scrooge = map2(
    scrooge_fit,
    map(prepped_fishery, "sampled_years"),
    process_scrooge
  )) %>%
  mutate(observed = map(prepped_fishery, "simed_fishery")) %>%
  mutate(scrooge_performance = pmap(
    list(observed = observed,
         predicted = processed_scrooge),
    judge_performance
  )) %>%
  mutate(
    scrooge_rec_performance = pmap(
      list(observed = observed,
           predicted = processed_scrooge),
      judge_performance,
      observed_variable = rec_dev,
      predicted_variable = "log_rec_dev_t"
    )
  )  %>%
  mutate(lcomps = map2(prepped_fishery, processed_scrooge, ~process_lcomps(.x, .y$n_tl))) %>%
  mutate(
    rmse = map_dbl(scrooge_performance, ~ .x$comparison_summary$rmse),
    bias =  map_dbl(scrooge_performance, ~ .x$comparison_summary$bias)
  ) %>%
  arrange(rmse)

true_f <- hurts$prepped_fishery[[1]]$simed_fishery %>%
  group_by(year) %>%
  summarise(f = mean(f)) %>%
  mutate(year = year - min(year) + 1)


compare_econ <- econ_hurts$processed_scrooge[[1]]$f_t %>%
  mutate(model = "With Economic Data") %>%
  bind_rows( noecon_hurts$processed_scrooge[[1]]$f_t %>% mutate(model = "W/O Economic Data")) %>%
  rename(f_hat = value) %>%
  left_join(true_f, by = "year") %>%
  ungroup()

trend_plot <- compare_econ %>%
  group_by(year, model) %>%
  mutate(rank_f = percent_rank(f_hat)) %>%
  summarise(mean_f_hat = median(f_hat),
            upper_75 = max(f_hat[rank_f <= 0.95]),
            lower_75 = min(f_hat[rank_f >= 0.05]),
            f = mean(f)) %>%
  ungroup() %>%
  ggplot() +
  geom_ribbon(aes(
    x = year,
    ymax = upper_75,
    ymin = lower_75,
    fill = model
  ), alpha = 0.25,
  color = "grey",
  show.legend = F) +
  geom_line(aes(year, mean_f_hat, color = model), linetype = 2, size = 1.25) +
  geom_line(aes(year, f), size = 1.25) +
  geom_point(aes(year, f), size = 4) +
  scale_color_manual(values = wes_palette("Royal1")) +
  scale_fill_manual(values = wes_palette("Royal1")) +
  labs(y = "Fishing Mortality",
       caption = "Ribbon is 90% Credible Interval")


bias_plot <- compare_econ %>%
  group_by(iteration, model) %>%
  summarise(rmse = sqrt(mean((f_hat - f)^2)),
            bias = median((f_hat - f) / f)) %>%
  ungroup() %>%
  group_by(model) %>%
  mutate(mean_bias = mean(bias),
         mean_rmse = mean(rmse)) %>%
  ggplot() +
  geom_vline(aes(xintercept = 0), size = 2) +
  geom_density(aes(bias, fill = model), alpha = 0.5, color = "grey", show.legend = F) +
  geom_vline(aes(xintercept = mean_bias, color = model), show.legend = F) +
  scale_x_continuous(labels = percent) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title = "Bias") +
  scale_color_manual(values = wes_palette("Royal1")) +
  scale_fill_manual(values = wes_palette("Royal1"))

rmse_plot <- compare_econ %>%
  group_by(iteration, model) %>%
  summarise(rmse = sqrt(mean((f_hat - f)^2)),
            bias = median((f_hat - f) / f)) %>%
  ungroup() %>%
  group_by(model) %>%
  mutate(mean_bias = mean(bias),
         mean_rmse = mean(rmse)) %>%
  ggplot() +
  geom_density(aes(rmse, fill = model), alpha = 0.5, color = NA, show.legend = F) +
  geom_vline(aes(xintercept = mean_rmse, color = model), show.legend = F) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title = "RMSE") +
  scale_color_manual(values = wes_palette("Royal1")) +
  scale_fill_manual(values = wes_palette("Royal1"))

hurts_plot <- ((rmse_plot / bias_plot) | trend_plot) + plot_layout(ncol = 2, widths = c(1, 2))

# how does it compare? ----------------------------------------------------


compare_sim <- fisheries_sandbox %>%
  filter(fleet_model == "open-access",
         sigma_r == max(sigma_r),
         sigma_effort == min(sigma_effort),
         price_cv == max(price_cv),
         cost_cv == max(cost_cv),
         steepness == min(steepness),
         obs_error == min(obs_error),
         b_v_bmsy_oa == min(b_v_bmsy_oa),
         q_cv == min(q_cv)) %>%
  slice(1)

compare_sim$summary_plot

compare_sim <- compare_sim %>%
  mutate(prepped_fishery = map(prepped_fishery, subsample_data, window = 20, period = "end")) %>%
  mutate(scrooge_fit = pmap(
    list(
      data = map(prepped_fishery, "scrooge_data"),
      fish = map(prepped_fishery, "fish"),
      fleet = map(prepped_fishery, "fleet")),
    fit_scrooge,
    iter = 4000,
    warmup = 2000,
    adapt_delta = 0.8,
    economic_model = 1,
    effort_data_weight = 1,
    scrooge_file = "scrooge",
    in_clouds = F,
    experiment = "pfo",
    max_f_v_fmsy_increase = 0.5,
    chains = 1,
    cv_effort = 0.5,
    max_expansion = 1.5
  ))

compare_sim <- compare_sim %>%
  mutate(processed_scrooge = map2(
    scrooge_fit,
    map(prepped_fishery, "sampled_years"),
    process_scrooge
  )) %>%
  mutate(observed = map(prepped_fishery, "simed_fishery")) %>%
  mutate(scrooge_performance = pmap(
    list(observed = observed,
         predicted = processed_scrooge),
    judge_performance
  )) %>%
  mutate(
    scrooge_rec_performance = pmap(
      list(observed = observed,
           predicted = processed_scrooge),
      judge_performance,
      observed_variable = rec_dev,
      predicted_variable = "log_rec_dev_t"
    )
  )  %>%
  mutate(lcomps = map2(prepped_fishery, processed_scrooge, ~process_lcomps(.x, .y$n_tl))) %>%
  mutate(
    rmse = map_dbl(scrooge_performance, ~ .x$comparison_summary$rmse),
    bias =  map_dbl(scrooge_performance, ~ .x$comparison_summary$bias)
  ) %>%
  arrange(rmse)


compare_sim <- compare_sim %>%
  mutate(lime_fit = pmap(list(
    data = map(prepped_fishery, "scrooge_data"),
    fish = map(prepped_fishery, "fish"),
    fleet = map(prepped_fishery, "fleet")
  ), fit_lime))

compare_sim <- compare_sim %>%
  mutate(processed_lime = map2(
    lime_fit,
    map(prepped_fishery, "sampled_years"),
    safely(process_lime)
  ))


# fit lbspr

compare_sim <- compare_sim %>%
  mutate(lbspr_fit = pmap(list(
    data = map(prepped_fishery, "scrooge_data"),
    fish = map(prepped_fishery, "fish"),
    fleet = map(prepped_fishery, "fleet")
  ), fit_lbspr))

compare_sim <- compare_sim %>%
  mutate(processed_lbspr = map2(
    lbspr_fit,
    map(prepped_fishery, "sampled_years"),
    process_lbspr
  ))



true_f <- compare_sim$prepped_fishery[[1]]$simed_fishery %>%
  group_by(year) %>%
  summarise(f = unique(f)) %>%
  mutate(year = year - min(year) + 1)


lime_fits <- compare_sim$processed_lime[[1]]$result %>%
  select(year, estimate,lower, upper) %>%
  mutate(model = "lime")

lbspr_fits <- compare_sim$processed_lbspr[[1]] %>%
  select(year, fm) %>%
  mutate(lower = NA, upper = NA, model = "lbspr") %>%
  mutate(estimate = fm * compare_sim$prepped_fishery[[1]]$fish$m) %>%
  select(year, estimate, lower, upper, model)


scrooge_fits <- compare_sim$processed_scrooge[[1]]$f_t %>%
  group_by(year) %>%
  mutate(rank_f = percent_rank(value)) %>%
  summarise(estimate = mean(value),
  upper = max(value[rank_f <= 0.975]),
            lower = min(value[rank_f >= 0.025])) %>%
  mutate(model = "scrooge") %>%
  select(year, estimate, lower, upper, model)


model_fits <- lbspr_fits %>%
  bind_rows(lime_fits) %>%
  bind_rows(scrooge_fits) %>%
  left_join(true_f, by = "year") %>%
  mutate(upper = pmin(upper, 2))

model_trends_plot <- model_fits %>%
  ggplot() +
  geom_ribbon(aes(year, ymin = lower, ymax = upper, fill = model),
              alpha = 0.5, color = "grey") +
  geom_line(aes(year, estimate, color = model), linetype = 2, size = 1.25) +
  geom_line(aes(year, f, color = "true"), color = "black") +
  geom_point(aes(year, f), size = 4, alpha = 0.8) +
  scale_fill_manual(values = wes_palette("Zissou1")) +
  labs(y = "Fishing Mortality", x = "Year", caption = "Black points are true values")



save(file = here::here("presentations","gaines-lab.Rdata"), easy, medium, toughest, f_lengths, r_lengths,model_trends_plot,
     helps_plot, hurts_plot,model_trends_plot,model_fits)

