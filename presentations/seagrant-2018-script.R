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
library(tweenr)
library(caret)
library(recipes)
library(gganimate)
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
    rec_ac = 0.6,
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
  filter(year < 40, length_bin > 30) %>%
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
  cost = 5,
  cost_cv =  0,
  cost_ac = 0,
  q_cv = 0,
  q_ac = 0,
  profit_lags = 4,
  length_50_sel = 50,
  theta = 0.1,
  fleet_model = "constant-effort",
  sigma_effort = 0.4,
  effort_ac = 0.8,
  initial_effort = 40
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
  summarise(
    effort = mean(effort),
    biomass = sum(biomass),
    ncaught = sum(numbers_caught)
  ) %>%
  ggplot(aes(year, biomass)) +
  geom_point()

f_lengths <- length_comps %>%
  filter(year < 40, length_bin > 30) %>%
  ggplot(aes(length_bin, year, height = scaled_numbers, group = year)) +
  geom_density_ridges(stat = "identity") +
  labs(x = "Length (cm)", title = "Proportional Length Distribution")

# what's the evidence?

set.seed(42)
fish <-
  create_fish(
    query_fishlife = T,
    scientific_name = "Scomber japonicus",
    mat_mode = "length",
    linf = 100,
    max_age = 30,
    time_step = 1,
    sigma_r = 0,
    rec_ac = 0,
    price = 10,
    price_cv = 0,
    price_ac = 0,
    steepness = 0.8
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
  fleet_model = "open-access",
  initial_effort = 100
)

fleet$e_msy <-  NA

fleet$p_msy <-  NA

sim <- sim_fishery(
  fish = fish,
  fleet = fleet,
  manager = create_manager(mpa_size = 0),
  num_patches = 1,
  sim_years = 200,
  burn_year = 50,
  time_step = fish$time_step,
  est_msy = T,
  tune_costs = T,
  b_v_bmsy_oa = 1.5
)

oa <- sim %>%
  group_by(year) %>%
  summarise(
    effort = mean(effort),
    biomass = sum(biomass),
    profits = sum(profits)
  ) %>%
  mutate(ease = "linear",
         id = 1) %>%
  ungroup()

oa_tween <- oa %>%
  tween_elements(., "year", "id", "ease", nframes = 500)


oa_plot <- oa_tween %>%
  filter(year < 200) %>%
  ggplot(aes(biomass, effort, frame = .frame)) +
  geom_path(
    aes(cumulative = T),
    size = 1,
    show.legend = F,
    color = "steelblue"
  ) +
  labs(x = "Biomass", y = "Fishing Effort")


gganimate::gganimate(
  oa_plot,
  filename = "presentations/oa.gif",
  interval = 0.05,
  title_frame = FALSE
)


gganimate::gganimate(
  oa_plot,
  filename = "presentations/small-oa.gif",
  interval = 0.05,
  title_frame = FALSE,
  ani.width = 200, ani.height = 300
)


# how much does it help ---------------------------------------------------

helps <- fisheries_sandbox %>%
  filter(
    fleet_model == "open-access",
    sigma_r == max(sigma_r),
    sigma_effort == min(sigma_effort),
    price_cv == max(price_cv),
    cost_cv == max(cost_cv),
    steepness == min(steepness),
    obs_error == min(obs_error),
    b_v_bmsy_oa == min(b_v_bmsy_oa),
    q_cv == max(q_cv),
    q_ac == max(q_ac)
  ) %>%
  slice(1)

helps$summary_plot

econ_helps <- helps %>%
  mutate(prepped_fishery = map(
    prepped_fishery,
    subsample_data,
    window = 20,
    period = "middle"
  )) %>%
  mutate(
    scrooge_fit = pmap(
      list(
        data = map(prepped_fishery, "scrooge_data"),
        fish = map(prepped_fishery, "fish"),
        fleet = map(prepped_fishery, "fleet")
      ),
      fit_scrooge,
      iter = 4000,
      warmup = 2000,
      adapt_delta = 0.8,
      economic_model = 1,
      effort_data_weight = 0,
      scrooge_file = "scrooge",
      in_clouds = F,
      experiment = "pfo",
      max_f_v_fmsy_increase = 0.5,
      chains = 1,
      cv_effort = 0.25,
      max_expansion = 1.5
    )
  )

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
  mutate(lcomps = map2(
    prepped_fishery,
    processed_scrooge,
    ~ process_lcomps(.x, .y$n_tl)
  )) %>%
  mutate(
    rmse = map_dbl(scrooge_performance, ~ .x$comparison_summary$rmse),
    bias =  map_dbl(scrooge_performance, ~ .x$comparison_summary$bias)
  ) %>%
  arrange(rmse)

noecon_helps <- helps %>%
  mutate(prepped_fishery = map(
    prepped_fishery,
    subsample_data,
    window = 20,
    period = "middle"
  )) %>%
  mutate(
    scrooge_fit = pmap(
      list(
        data = map(prepped_fishery, "scrooge_data"),
        fish = map(prepped_fishery, "fish"),
        fleet = map(prepped_fishery, "fleet")
      ),
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
      cv_effort = 0.25,
      max_expansion = 1.5
    )
  )

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
  mutate(lcomps = map2(
    prepped_fishery,
    processed_scrooge,
    ~ process_lcomps(.x, .y$n_tl)
  )) %>%
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
  bind_rows(noecon_helps$processed_scrooge[[1]]$f_t %>% mutate(model = "W/O Economic Data")) %>%
  rename(f_hat = value) %>%
  left_join(true_f, by = "year") %>%
  ungroup()

trend_plot <- compare_econ %>%
  group_by(year, model) %>%
  mutate(rank_f = percent_rank(f_hat)) %>%
  summarise(
    mean_f_hat = median(f_hat),
    upper_75 = max(f_hat[rank_f <= 0.95]),
    lower_75 = min(f_hat[rank_f >= 0.05]),
    f = mean(f)
  ) %>%
  ungroup() %>%
  ggplot() +
  geom_ribbon(
    aes(
      x = year,
      ymax = upper_75,
      ymin = lower_75,
      fill = model
    ),
    alpha = 0.25,
    color = "grey",
    show.legend = F
  ) +
  geom_line(aes(year, mean_f_hat, color = model),
            linetype = 2,
            size = 1.25) +
  geom_line(aes(year, f), size = 1.25) +
  geom_point(aes(year, f), size = 4) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(y = "Fishing Mortality",
       caption = "Ribbon is 90% Credible Interval", x = "Year") +
  scale_color_manual(name = "", values = wes_palette("Royal1")) +
  scale_fill_manual(name = "", values = wes_palette("Royal1")) +
  theme(plot.margin = unit(c(.1, .1, .1, .1), "lines"),
        legend.position = "top")


bias_plot <- compare_econ %>%
  group_by(iteration, model) %>%
  summarise(rmse = sqrt(mean((f_hat - f) ^ 2)),
            bias = pmin(4, median((f_hat - f) / f))) %>%
  ungroup() %>%
  group_by(model) %>%
  mutate(mean_bias = mean(bias),
         mean_rmse = mean(rmse)) %>%
  ggplot() +
  geom_vline(aes(xintercept = 0), size = 1, linetype = 2) +
  geom_density(
    aes(bias, fill = model),
    alpha = 0.5,
    color = NA,
    show.legend = F
  ) +
  # geom_vline(aes(xintercept = mean_bias, color = model)) +
  scale_x_continuous(labels = percent) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Bias") +
  theme(
    legend.position = "bottom",
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(.1, .1, .1, .1), "lines")
  ) +
  scale_color_manual(values = wes_palette("Royal1")) +
  scale_fill_manual(values = wes_palette("Royal1"))

rmse_plot <- compare_econ %>%
  group_by(iteration, model) %>%
  summarise(rmse = sqrt(mean((f_hat - f) ^ 2)),
            bias = median((f_hat - f) / f)) %>%
  ungroup() %>%
  group_by(model) %>%
  mutate(mean_bias = mean(bias),
         mean_rmse = mean(rmse)) %>%
  ggplot() +
  geom_density(
    aes(rmse, fill = model),
    alpha = 0.5,
    color = NA,
    show.legend = F
  ) +
  geom_vline(aes(xintercept = 0), size = 1, linetype = 2) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "RMSE") +
  scale_color_manual(values = wes_palette("Royal1")) +
  scale_fill_manual(values = wes_palette("Royal1")) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(.1, .1, .1, .1), "lines")
  )

helps_plot <-
  (trend_plot |
     (rmse_plot / bias_plot)) + plot_layout(ncol = 2, widths = c(2, 1))


# hurts -------------------------------------------------------------------

hurts <- fisheries_sandbox %>%
  filter(
    fleet_model == "constant-effort",
    sigma_r == max(sigma_r),
    sigma_effort == min(sigma_effort),
    price_cv == max(price_cv),
    cost_cv == max(cost_cv),
    steepness == max(steepness),
    obs_error == min(obs_error),
    b_v_bmsy_oa == min(b_v_bmsy_oa),
    q_cv == max(q_cv),
    q_ac == max(q_ac)
  ) %>%
  slice(1)

econ_hurts <- hurts %>%
  mutate(prepped_fishery = map(
    prepped_fishery,
    subsample_data,
    window = 10,
    period = "end"
  )) %>%
  mutate(
    scrooge_fit = pmap(
      list(
        data = map(prepped_fishery, "scrooge_data"),
        fish = map(prepped_fishery, "fish"),
        fleet = map(prepped_fishery, "fleet")
      ),
      fit_scrooge,
      iter = 4000,
      warmup = 2000,
      adapt_delta = 0.8,
      economic_model = 1,
      effort_data_weight = 0,
      scrooge_file = "scrooge",
      in_clouds = F,
      experiment = "pfo",
      max_f_v_fmsy_increase = 0.5,
      chains = 1,
      cv_effort = 0.25,
      max_expansion = 1.5
    )
  )

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
  mutate(lcomps = map2(
    prepped_fishery,
    processed_scrooge,
    ~ process_lcomps(.x, .y$n_tl)
  )) %>%
  mutate(
    rmse = map_dbl(scrooge_performance, ~ .x$comparison_summary$rmse),
    bias =  map_dbl(scrooge_performance, ~ .x$comparison_summary$bias)
  ) %>%
  arrange(rmse)

noecon_hurts <- hurts %>%
  mutate(prepped_fishery = map(
    prepped_fishery,
    subsample_data,
    window = 10,
    period = "end"
  )) %>%
  mutate(
    scrooge_fit = pmap(
      list(
        data = map(prepped_fishery, "scrooge_data"),
        fish = map(prepped_fishery, "fish"),
        fleet = map(prepped_fishery, "fleet")
      ),
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
      cv_effort = 0.25,
      max_expansion = 1.5
    )
  )

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
  mutate(lcomps = map2(
    prepped_fishery,
    processed_scrooge,
    ~ process_lcomps(.x, .y$n_tl)
  )) %>%
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
  bind_rows(noecon_hurts$processed_scrooge[[1]]$f_t %>% mutate(model = "W/O Economic Data")) %>%
  rename(f_hat = value) %>%
  left_join(true_f, by = "year") %>%
  ungroup()

trend_plot <- compare_econ %>%
  group_by(year, model) %>%
  mutate(rank_f = percent_rank(f_hat)) %>%
  summarise(
    mean_f_hat = median(f_hat),
    upper_75 = max(f_hat[rank_f <= 0.95]),
    lower_75 = min(f_hat[rank_f >= 0.05]),
    f = mean(f)
  ) %>%
  ungroup() %>%
  ggplot() +
  geom_ribbon(
    aes(
      x = year,
      ymax = upper_75,
      ymin = lower_75,
      fill = model
    ),
    alpha = 0.25,
    color = "grey",
    show.legend = F
  ) +
  geom_line(aes(year, mean_f_hat, color = model),
            linetype = 2,
            size = 1.25) +
  geom_line(aes(year, f), size = 1.25) +
  geom_point(aes(year, f), size = 4) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(y = "Fishing Mortality",
       caption = "Ribbon is 90% Credible Interval", x = "Year") +
  scale_color_manual(name = "", values = wes_palette("Royal1")) +
  scale_fill_manual(name = "", values = wes_palette("Royal1")) +
  theme(plot.margin = unit(c(.1, .1, .1, .1), "lines"),
        legend.position = "top")


bias_plot <- compare_econ %>%
  group_by(iteration, model) %>%
  summarise(rmse = sqrt(mean((f_hat - f) ^ 2)),
            bias = median((f_hat - f) / f)) %>%
  ungroup() %>%
  group_by(model) %>%
  mutate(mean_bias = mean(bias),
         mean_rmse = mean(rmse)) %>%
  ggplot() +
  geom_vline(aes(xintercept = 0), size = 1, linetype = 2) +
  geom_density(
    aes(bias, fill = model),
    alpha = 0.5,
    color = NA,
    show.legend = F
  ) +
  # geom_vline(aes(xintercept = mean_bias, color = model)) +
  scale_x_continuous(labels = percent) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Bias") +
  theme(
    legend.position = "bottom",
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(.1, .1, .1, .1), "lines")
  ) +
  scale_color_manual(values = wes_palette("Royal1")) +
  scale_fill_manual(values = wes_palette("Royal1"))

rmse_plot <- compare_econ %>%
  group_by(iteration, model) %>%
  summarise(rmse = sqrt(mean((f_hat - f) ^ 2)),
            bias = median((f_hat - f) / f)) %>%
  ungroup() %>%
  group_by(model) %>%
  mutate(mean_bias = mean(bias),
         mean_rmse = mean(rmse)) %>%
  ggplot() +
  geom_density(
    aes(rmse, fill = model),
    alpha = 0.5,
    color = NA,
    show.legend = F
  ) +
  geom_vline(aes(xintercept = 0), size = 1, linetype = 2) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "RMSE") +
  scale_color_manual(values = wes_palette("Royal1")) +
  scale_fill_manual(values = wes_palette("Royal1")) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(.1, .1, .1, .1), "lines")
  )

hurts_plot <-
  (trend_plot |
     (rmse_plot / bias_plot)) + plot_layout(ncol = 2, widths = c(2, 1))


# how does it compare? ----------------------------------------------------


compare_sim <- fisheries_sandbox %>%
  filter(
    fleet_model == "open-access",
    sigma_r == max(sigma_r),
    sigma_effort == min(sigma_effort),
    price_cv == max(price_cv),
    cost_cv == max(cost_cv),
    steepness == min(steepness),
    obs_error == min(obs_error),
    b_v_bmsy_oa == min(b_v_bmsy_oa),
    q_cv == max(q_cv),
    q_ac == max(q_ac)
  ) %>%
  slice(1)

compare_sim$summary_plot

compare_sim <- compare_sim %>%
  mutate(prepped_fishery = map(
    prepped_fishery,
    subsample_data,
    window = 20,
    period = "middle"
  )) %>%
  mutate(
    scrooge_fit = pmap(
      list(
        data = map(prepped_fishery, "scrooge_data"),
        fish = map(prepped_fishery, "fish"),
        fleet = map(prepped_fishery, "fleet")
      ),
      fit_scrooge,
      iter = 4000,
      warmup = 2000,
      adapt_delta = 0.8,
      economic_model = 1,
      effort_data_weight = 0,
      scrooge_file = "scrooge",
      in_clouds = F,
      experiment = "pfo",
      max_f_v_fmsy_increase = 0.5,
      chains = 1,
      cv_effort = 0.25,
      max_expansion = 1.5
    )
  )

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
  mutate(lcomps = map2(
    prepped_fishery,
    processed_scrooge,
    ~ process_lcomps(.x, .y$n_tl)
  )) %>%
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
  select(year, estimate, lower, upper) %>%
  mutate(model = "lime") %>%
  mutate(upper = pmin(upper, 1))

lbspr_fits <- compare_sim$processed_lbspr[[1]] %>%
  select(year, fm) %>%
  mutate(lower = NA,
         upper = NA,
         model = "lbspr") %>%
  mutate(estimate = fm * compare_sim$prepped_fishery[[1]]$fish$m) %>%
  select(year, estimate, lower, upper, model)


scrooge_fits <- compare_sim$processed_scrooge[[1]]$f_t %>%
  group_by(year) %>%
  mutate(rank_f = percent_rank(value)) %>%
  summarise(
    estimate = mean(value),
    upper = max(value[rank_f <= 0.975]),
    lower = min(value[rank_f >= 0.025])
  ) %>%
  mutate(model = "scrooge") %>%
  select(year, estimate, lower, upper, model)


model_fits <- lbspr_fits %>%
  bind_rows(lime_fits) %>%
  bind_rows(scrooge_fits) %>%
  left_join(true_f, by = "year") %>%
  mutate(upper = pmin(upper, 2))

model_trends_plot <- model_fits %>%
  ggplot() +
  geom_ribbon(aes(
    year,
    ymin = lower,
    ymax = upper,
    fill = model
  ),
  alpha = 0.5,
  color = "grey") +
  geom_line(aes(year, estimate, color = model),
            linetype = 2,
            size = 1.25) +
  geom_line(aes(year, f, color = "true"), color = "black") +
  geom_point(aes(year, f), size = 4, alpha = 0.8) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  lims(y = c(0, 1)) +
  # scale_fill_manual(values = wes_palette("Zissou1")) +
  labs(y = "Fishing Mortality", x = "Year", caption = "Black points are true values")



# when does it work? ------------------------------------------------------
run_name <- "v1.0"
in_clouds <-  T
if (in_clouds == T) {
  system("umount results/scrooge-results")

  system("rm -r results/scrooge-results")

  if (dir.exists("results/scroote-results") == F) {
    system("mkdir results/scrooge-results")

  }

  system("gcsfuse scrooge-results results/scrooge-results")

  cloud_dir <- here::here("results", "scrooge-results", run_name)

}

experiment_files <- list.files(cloud_dir)

experiment_files <-
  experiment_files[!str_detect(experiment_files, "description")]

experiment_numbers <-
  str_replace_all(experiment_files, "\\D", "") %>% as.numeric()

experiments <-
  expand.grid(
    period = c("middle", "end"),
    window = c(2, 5, 10),
    economic_model = c(0, 1),
    effort_data_weight = c(0, 1),
    fishery = 1:nrow(fisheries_sandbox),
    stringsAsFactors = F
  )

fisheries_sandbox <- fisheries_sandbox %>%
  select(-economic_model) %>%
  mutate(fishery = 1:nrow(.)) %>%
  left_join(experiments, by = "fishery") %>%
  mutate(prepped_fishery = pmap(
    list(
      prepped_fishery = prepped_fishery,
      window = window,
      period = period
    ),
    subsample_data
  ))

fisheries_sandbox$experiment <- 1:nrow(fisheries_sandbox)

# map(prepped_fishery, "sampled_years"),


experiment_sandbox <- fisheries_sandbox %>%
  filter(experiment %in% experiment_numbers) %>%
  arrange(experiment) %>%
  mutate(
    performance = map2(experiment, prepped_fishery, process_experiment, cloud_dir = cloud_dir)
  )


experiment_summary <- experiment_sandbox %>%
  select(-fleet_params,-prepped_fishery,-summary_plot) %>%
  unnest() %>%
  mutate(emodel = ifelse(economic_model == 1, "Uses Economic Data", "Ignores Economic Data")) %>%
  mutate(eweight = ifelse(effort_data_weight == 1, "Effort Priors", "Incentive Priors"))


rmse_plot <- experiment_summary %>%
  ggplot(aes(pmin(1, rmse), fill = eweight)) +
  geom_vline(aes(xintercept = 0), linetype = 2) +
  geom_density(alpha = 0.5) +
  facet_grid(emodel ~ fleet_model, scales = "free") +
  scale_fill_viridis_d(name = "") +
  labs(x = 'RMSE', caption = "Note: Random sample of experiments") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank())

bias_plot <- experiment_summary %>%
  ggplot(aes(pmin(1, bias), fill = eweight)) +
  geom_vline(aes(xintercept = 0), linetype = 2) +
  geom_density(alpha = 0.5) +
  facet_grid(emodel ~ fleet_model, scales = "free") +
  scale_x_continuous(labels = percent, name = "Bias") +
  scale_fill_viridis_d(name = "") +
  labs(caption = "Note: Random sample of experiments") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank())


# machine learning --------------------------------------------------------

linf <- fisheries_sandbox %>%
  select(sci_name, prepped_fishery) %>%
  mutate(linf = map_dbl(prepped_fishery, c("fish", "linf"))) %>%
  select(-prepped_fishery) %>%
  group_by(sci_name) %>%
  summarise(linf = mean(linf))

experiment_frame <- experiment_summary %>%
  left_join(linf, by = "sci_name") %>%
  select(
    rmse,
    experiment,
    economic_model,
    effort_data_weight,
    sigma_r,
    linf,
    steepness,
    q_cv,
    price_cv,
    cost_cv
  ) %>%
  filter(!(economic_model == 0 & effort_data_weight == 1)) %>%
  mutate(economic_model = ifelse(economic_model == 1, "Uses Economic Data", "Ignores Economic Data")) %>%
  mutate(effort_data_weight = ifelse(effort_data_weight == 1, "Effort Priors", "Incentive Priors")) %>%
  mutate(model = glue::glue("{economic_model}-{effort_data_weight}")) %>%
  select(-economic_model,-effort_data_weight) %>%
  filter(rmse < 1)

scrooge_split <-
  rsample::initial_split(experiment_frame, strata = "model")

scrooge_train <- rsample::training(scrooge_split)

scrooge_test <- rsample::testing(scrooge_split)


fitfoo <- function(data) {
  data <- data %>%
    select(rmse, sigma_r, cost_cv, price_cv, q_cv, linf)

  scrooge_recipe <- recipe(rmse ~ ., data = data)  %>%
    step_num2factor(all_predictors())

  prepped_scrooge_recipe <-
    prep(scrooge_recipe, data = data, retain = TRUE)

  model <- caret::train(
    scrooge_recipe,
    data = data,
    method = "ranger",
    importance = "impurity_corrected"
  )

}

scrooge_train <- scrooge_train %>%
  nest(-model) %>%
  mutate(scrooge_fit = map(data, fitfoo))

get_varimp <- function(model) {
  varimp <-  caret::varImp(model)$importance %>%
    as.data.frame() %>%
    mutate(variable = rownames(.)) %>%
    rename(importance = Overall)


}


scrooge_train$varimp <- map(scrooge_train$scrooge_fit, get_varimp)

varimp_plot <- scrooge_train %>%
  select(model, varimp) %>%
  unnest() %>%
  ggplot(aes(variable, importance, fill = model)) +
  geom_col(position = "dodge", color = "grey") +
  scale_fill_viridis_d(name = "") +
  labs(y = "Importance Score") +
  theme(axis.title.y = element_blank()) +
  # facet_wrap(~model) +
  coord_flip()



add_preds <- function(model, newdata) {
  newdata <- newdata %>%
    select(rmse, experiment, sigma_r, cost_cv, price_cv, q_cv, linf, model)

  preds <- predict(model, newdata = newdata)

  newdata$pred_rmse <- preds

  return(newdata)

}

scrooge_preds <- scrooge_train %>%
  mutate(pred_test = map(scrooge_fit, add_preds, scrooge_test)) %>%
  rename(fitted_model = model) %>%
  select(-data,-scrooge_fit, -varimp) %>%
  unnest() %>%
  arrange(rmse) %>%
  mutate(fishery = factor(rmse))


scrooge_preds <- scrooge_preds %>%
  mutate(fishery = as.numeric(factor(fishery)))

scrooge_pred_rmse_plot <- scrooge_preds %>%
  ggplot(aes(factor(fishery), pred_rmse, fill = fitted_model), alpha = 0.75) +
  geom_point(shape = 21) +
  coord_flip() +
  labs(x = "Fishery", y = "Predicted RMSE") +
  scale_fill_viridis_d(name = "")


save(
  file = here::here("presentations", "seagrant-presentation.Rdata"),
  f_lengths,
  r_lengths,
  model_trends_plot,
  helps_plot,
  hurts_plot,
  model_fits,
  rmse_plot,
  bias_plot,
  experiment_sandbox,
  varimp_plot,
  scrooge_pred_rmse_plot
)
