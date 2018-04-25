set.seed(42)
rstan::rstan_options(auto_write = TRUE)

custom_sandbox <- F

if (custom_sandbox == T){
fisheries_sandbox <-
  purrr::cross_df(
    list(
      sci_name = c("Lutjanus campechanus"),
      fleet_model = c("open-access"),
      sigma_r = c(0.4),
      sigma_effort = c(0.1),
      price_cv = c(0.25),
      cost_cv = c(0.5),
      price_ac = 0.5,
      cost_ac = 0.5,
      q_cv = .2,
      q_ac = 0.5,
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
      profit_lags = 3,
      query_price = F,
      use_effort_data = 1,
      cv_len = 0.1,
      cv_effort = 1.6
    )
  )

no_error <- map(fisheries_sandbox$prepped_fishery, 'error') %>%
  map_lgl(is.null)

fisheries_sandbox <- fisheries_sandbox %>%
  filter(no_error == T)

fisheries_sandbox <- fisheries_sandbox %>%
  mutate(prepped_fishery = map(prepped_fishery, "result")) %>%
  mutate(summary_plot = map(prepped_fishery, plot_simmed_fishery))
} else{

  load(file = here::here("processed_data", "fisheries_sandbox.Rdata"))

}


test_set <- fisheries_sandbox %>%
  dplyr::filter(
    fleet_model == "open-access",
    sigma_r == max(sigma_r),
    sigma_effort == min(sigma_effort),
    price_cv == max(price_cv),
    price_ac == max(price_ac),
    cost_ac == max(cost_ac),
    cost_cv == max(cost_cv),
    q_cv == max(q_cv),
    q_ac == max(q_ac),
    obs_error == min(obs_error),
    b_v_bmsy_oa == 0.5
  ) %>%
  slice(2)


tester <- test_set %>%
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
      scrooge_file = "scrooge",
      in_clouds = F,
      experiment = "pfo",
      max_f_v_fmsy_increase = 0.5,
      chains = 1,
      cv_effort = 0.5,
      max_expansion = 1.5
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

tester$scrooge_performance[[1]]$comparison_plot

tester$scrooge_rec_performance[[1]]$comparison_plot


scrooge_fit <- tester$scrooge_fit[[1]]


profits <- rstan::extract(scrooge_fit, "profit_t")$profit_t %>%
  as_data_frame() %>%
  mutate(iteration = 1:nrow(.)) %>%
  gather(year, profits,-iteration) %>%
  mutate(year = str_replace(year,"\\D","") %>% as.numeric()) %>%
  group_by(year) %>%
  summarise(mean_profits = mean(profits)) %>%
  ggplot(aes(year, mean_profits)) +
  geom_point()



numbers <- rstan::extract(scrooge_fit, "n_ta")$n_ta %>%
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


numbers %>%
  group_by(year) %>%
  filter(age == 1) %>%
  summarise(tn = sum(numbers)) %>%
  ggplot(aes(year, tn)) +
  geom_point()

selectivity <- rstan::extract(scrooge_fit, "mean_selectivity_at_age")$mean_selectivity_at_age %>%
  purrr::array_branch(1) %>%
  map(~as_data_frame(t(.x))) %>%
  bind_rows(.id = "iteration") %>%
  group_by(iteration) %>%
  mutate(year = 1:length(iteration)) %>%
  gather(age,selectivity, -iteration,-year) %>%
  ungroup() %>%
  mutate(age = str_replace(age,"\\D","") %>% as.numeric()) %>%
  group_by(year, age) %>%
  summarise(selectivity = mean(selectivity)) %>%
  mutate(true_selectivity = tester$prepped_fishery[[1]]$fleet$sel_at_age,
         mean_length = tester$prepped_fishery[[1]]$fish$length_at_age)


sel_50_hat <- rstan::extract(scrooge_fit, "length_50_sel")$length_50_sel %>%
  mean()

sel_95_hat <- sel_50_hat + 2

sel_50_true <- tester$prepped_fishery[[1]]$fleet$length_50_sel

manual <- 1.0 / (1 + exp(-log(19) * ((tester$prepped_fishery[[1]]$fish$length_at_age
 - sel_50_hat) / 2)))

1 / (1 + exp(-log(19)*(tester$prepped_fishery[[1]]$scrooge_data$mean_length_at_age - sel_95_hat)/(2)))

1 / (1 + exp(-log(19)*(tester$prepped_fishery[[1]]$fish$length_at_age - sel_50_true)/(2)))

tester$prepped_fishery[[1]]$fleet$sel_at_age

selectivity$manual <- manual


tester$prepped_fishery[[1]]$scrooge_data$mean_length_at_age

tester$prepped_fishery[[1]]$fish$length_at_age

selectivity %>%
  ggplot() +
  geom_line(aes(mean_length, selectivity)) +
  geom_line(aes(mean_length, manual), color = "blue") +
  geom_point(aes(mean_length, true_selectivity)) +
  geom_vline(aes(xintercept = sel_50_hat), alpha = 0.5) +
  geom_vline(aes(xintercept = sel_50_true), color = "red", alpha = 0.5)


tester$prepped_fishery[[1]]$fleet$length_50_sel

tester$prepped_fishery[[1]]$fleet$length_95_sel

numbers %>%
  group_by(year) %>%
  filter(age == 1) %>%
  summarise(tn = sum(numbers)) %>%
  ggplot(aes(year, tn)) +
  geom_point()

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
  summarise(tn = sum(numbers)) %>%
  ggplot(aes(year, tn)) +
  geom_point()

r_lengths <- burn %>%
  ggplot(aes(age, year, height = numbers, group = year)) +
  geom_density_ridges(stat = "identity") +
  labs(x = "Length (cm)", title = "Proportional Length Distribution")



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
  mutate(estimate = fm * tester$prepped_fishery[[1]]$fish$m) %>%
  select(year, estimate, lower, upper, model)


scrooge_fits <- tester$processed_scrooge[[1]]$f_t %>%
  group_by(year) %>%
  mutate(rank_f = percent_rank(value)) %>%
  summarise(
    estimate = median(value),
    upper = max(value[rank_f <= 0.975]),
    lower = min(value[rank_f >= 0.025])
  ) %>%
  mutate(model = "scrooge") %>%
  select(year, estimate, lower, upper, model)


model_fits <- lbspr_fits %>%
  bind_rows(lime_fits) %>%
  bind_rows(scrooge_fits) %>%
  left_join(true_f, by = "year")

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
  geom_line(aes(year, f), color = "black", alpha = 0.75) +
  geom_point(aes(year, f), size = 4, alpha = 0.75) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  labs(y = "Fishing Mortality", x = "Year", caption = "Black points are true values")

