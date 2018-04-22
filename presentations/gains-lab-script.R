library(rstan)
library(spasm)
library(FishLife)
library(tidyverse)


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
    iter = 8000,
    warmup = 4000,
    adapt_delta = 0.8,
    economic_model = 1,
    use_effort_data = 0,
    scrooge_file = "scrooge",
    in_clouds = F,
    experiment = "pfo",
    max_f_v_fmsy_increase = 0.5
  )
  )

easy <- easy %>%
  mutate(lime_fit = pmap(list(
    data = map(prepped_fishery, "scrooge_data"),
    fish = map(prepped_fishery, "fish"),
    fleet = map(prepped_fishery, "fleet")
  ), fit_lime))

easy <- easy %>%
  mutate(processed_lime = map2(
    lime_fit,
    map(prepped_fishery, "sampled_years"),
    safely(process_lime)
  ))


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
      predicted_variable = "rec_dev_t"
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
         obs_error == max(obs_error),
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
    iter = 2000,
    warmup = 1000,
    adapt_delta = 0.8,
    economic_model = 1,
    use_effort_data = 1,
    scrooge_file = "scrooge",
    in_clouds = F,
    experiment = "pfo",
    max_f_v_fmsy_increase = 0.5
  )
  )

medium <- medium %>%
  mutate(lime_fit = pmap(list(
    data = map(prepped_fishery, "scrooge_data"),
    fish = map(prepped_fishery, "fish"),
    fleet = map(prepped_fishery, "fleet")
  ), fit_lime))

medium <- medium %>%
  mutate(processed_lime = map2(
    lime_fit,
    map(prepped_fishery, "sampled_years"),
    safely(process_lime)
  ))


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
      predicted_variable = "rec_dev_t"
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
         obs_error == max(obs_error),
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
    iter = 2000,
    warmup = 1000,
    adapt_delta = 0.8,
    economic_model = 1,
    use_effort_data = 0,
    scrooge_file = "scrooge",
    in_clouds = F,
    experiment = "pfo",
    max_f_v_fmsy_increase = 0.5
  )
  )

toughest <- toughest %>%
  mutate(lime_fit = pmap(list(
    data = map(prepped_fishery, "scrooge_data"),
    fish = map(prepped_fishery, "fish"),
    fleet = map(prepped_fishery, "fleet")
  ), fit_lime))

toughest <- toughest %>%
  mutate(processed_lime = map2(
    lime_fit,
    map(prepped_fishery, "sampled_years"),
    safely(process_lime)
  ))


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
      predicted_variable = "rec_dev_t"
    )
  )  %>%
  mutate(lcomps = map2(prepped_fishery, processed_scrooge, ~process_lcomps(.x, .y$n_tl))) %>%
  mutate(
    rmse = map_dbl(scrooge_performance, ~ .x$comparison_summary$rmse),
    bias =  map_dbl(scrooge_performance, ~ .x$comparison_summary$bias)
  ) %>%
  arrange(rmse)

save(file = here::here("presentations","gaines-lab.Rdata"), easy, medium, toughest)

