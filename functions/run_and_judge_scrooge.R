run_and_judge_scrooge <- function(fishery,
                                  total_iterations = 4000,
                                  warmup = 2000,
                                  adapt_delta = 0.8,
                                  chains = 1,
                                  economic_model = 1,
                                  effort_data_weight = 0,
                                  scrooge_file = "scrooge",
                                  experiment_name = "blah",
                                  max_f_v_fmsy_increase = 0.5,
                                  cv_effort = 0.25,
                                  max_expansion = 1.5,
                                  window = 10,
                                  period = "end",
                                  prop_years_lcomp_data = 1
                                  )
{

  assessed_fishery <- fishery %>%
  mutate(prepped_fishery = map(
    prepped_fishery,
    subsample_data,
    window = window,
    period = period,
    prop_years_lcomp_data = prop_years_lcomp_data
  )) %>%
  mutate(
    scrooge_fit = pmap(
      list(
        data = map(prepped_fishery, "scrooge_data"),
        fish = map(prepped_fishery, "fish"),
        fleet = map(prepped_fishery, "fleet")
      ),
      fit_scrooge,
      iter = total_iterations,
      warmup = warmup,
      adapt_delta = adapt_delta,
      economic_model = economic_model,
      effort_data_weight = effort_data_weight,
      scrooge_file = scrooge_file,
      in_clouds = F,
      experiment = experiment_name,
      max_f_v_fmsy_increase = max_f_v_fmsy_increase,
      chains = chains,
      cv_effort = cv_effort,
      max_expansion = max_expansion
    )
  )

  assessed_fishery <- assessed_fishery %>%
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

return(assessed_fishery)
}