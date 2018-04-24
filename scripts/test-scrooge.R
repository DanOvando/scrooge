fisheries_sandbox <-
  purrr::cross_df(
    list(
      sci_name = c("Lutjanus campechanus"),
      fleet_model = c("open-access"),
      sigma_r = c(0),
      sigma_effort = c(0),
      price_cv = c(0),
      cost_cv = c(0),
      price_ac = 0,
      cost_ac = 0,
      q_cv = 0,
      q_ac = 0,
      economic_model = c(1),
      steepness = c(0.6),
      obs_error = c(0),
      b_v_bmsy_oa = c(0.5)
    )
  )

fleet_model_params <- data_frame(
  fleet_model = c(
    "constant-catch",
    "constant-effort",
    "supplied-catch",
    "open-access"
  ),
  fleet_params = list(
    list(target_catch = 10000),
    list(initial_effort = 200),
    list(catches = cdfw_catches$catch[cdfw_catches$sci_name == "semicossyphus pulcher"]),
    list(theta = 0.5, initial_effort = 10)
  )
)

fisheries_sandbox <- fisheries_sandbox %>%
  left_join(fleet_model_params, by = "fleet_model") %>%
  mutate(
    prepped_fishery = pmap(
      list(
        sci_name = sci_name,
        fleet_model = fleet_model,
        fleet_params = fleet_params,
        sigma_r = sigma_r,
        sigma_effort = sigma_effort,
        price_cv = price_cv,
        cost_cv = cost_cv,
        price_ac = price_ac,
        cost_ac = cost_ac,
        steepness = steepness,
        obs_error = obs_error,
        b_v_bmsy_oa = b_v_bmsy_oa,
        q_cv = q_cv,
        q_ac = q_ac
      ),
      safely(prepare_fishery),
      sim_years = 200,
      burn_years = 50,
      price = 10,
      cost = 5,
      profit_lags = 1,
      query_price = F,
      use_effort_data = 1
    )
  )

no_error <- map(fisheries_sandbox$prepped_fishery, 'error') %>%
  map_lgl(is.null)

fisheries_sandbox <- fisheries_sandbox %>%
  filter(no_error == T)

fisheries_sandbox <- fisheries_sandbox %>%
  mutate(prepped_fishery = map(prepped_fishery, "result")) %>%
  mutate(summary_plot = map(prepped_fishery, plot_simmed_fishery))


tester <- fisheries_sandbox %>%
  mutate(prepped_fishery = map(
    prepped_fishery,
    subsample_data,
    window = 20,
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
      use_effort_data = 0,
      scrooge_file = "new_scrooge",
      in_clouds = F,
      experiment = "pfo",
      max_f_v_fmsy_increase = 0.05
    )
  )

tester <- tester %>%
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



scrooge_fit <- tester$scrooge_fit[[1]]

burn <- rstan::extract(scrooge_fit, "n_a_init")$n_a_init %>%
  purrr::array_branch(1) %>%
  map(as_data_frame) %>%
  bind_rows(.id = "iteration") %>%
  group_by(iteration) %>%
  mutate(year = 1:length(iteration)) %>%
  gather(age,numbers, -iteration,-year) %>%
  ungroup() %>%
  mutate(age = str_replace(age,"\\D","") %>% as.numeric()) %>%
  group_by(year, age) %>%
  summarise(numbers = mean(numbers))

burn %>%
  group_by(year) %>%
  filter(age == min(age)) %>%
  summarise(tn = sum(numbers)) %>%
  ggplot(aes(year, tn)) +
  geom_point()

r_lengths <- burn %>%
  filter(year > 90) %>%
  ggplot(aes(age, year, height = numbers, group = year)) +
  geom_density_ridges(stat = "identity") +
  labs(x = "Length (cm)", title = "Proportional Length Distribution")



tester$scrooge_performance[[1]]$comparison_plot +
  lims(x = c(35, NA))

tester$scrooge_rec_performance[[1]]$comparison_plot+
  lims(x = c(35, NA))


tester <- tester %>%
  mutate(lime_fit = pmap(list(
    data = map(prepped_fishery, "scrooge_data"),
    fish = map(prepped_fishery, "fish"),
    fleet = map(prepped_fishery, "fleet")
  ), fit_lime))

tester <- tester %>%
  mutate(processed_lime = map2(
    lime_fit,
    map(prepped_fishery, "sampled_years"),
    safely(process_lime)
  ))


# fit lbspr

tester <- tester %>%
  mutate(lbspr_fit = pmap(list(
    data = map(prepped_fishery, "scrooge_data"),
    fish = map(prepped_fishery, "fish"),
    fleet = map(prepped_fishery, "fleet")
  ), fit_lbspr))

tester <- tester %>%
  mutate(processed_lbspr = map2(
    lbspr_fit,
    map(prepped_fishery, "sampled_years"),
    process_lbspr
  ))



true_f <- tester$prepped_fishery[[1]]$simed_fishery %>%
  group_by(year) %>%
  summarise(f = unique(f)) %>%
  mutate(year = year - min(year) + 1)


lime_fits <- tester$processed_lime[[1]]$result %>%
  select(year, estimate, lower, upper) %>%
  mutate(model = "lime")

lbspr_fits <- tester$processed_lbspr[[1]] %>%
  select(year, fm) %>%
  mutate(lower = NA,
         upper = NA,
         model = "lbspr") %>%
  mutate(estimate = fm * compare_sim$prepped_fishery[[1]]$fish$m) %>%
  select(year, estimate, lower, upper, model)


scrooge_fits <- tester$processed_scrooge[[1]]$f_t %>%
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
  scale_fill_manual(values = wes_palette("Zissou1")) +
  labs(y = "Fishing Mortality", x = "Year", caption = "Black points are true values")
