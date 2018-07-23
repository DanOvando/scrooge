

experiments <- expand.grid(period = c("end"), window = c(20),
                           economic_model = c(0),
                           fishery = 1:nrow(fisheries_sandbox), stringsAsFactors = F)

fisheries_sandbox <- fisheries_sandbox %>%
  mutate(fishery = 1:nrow(.)) %>%
  slice(3) %>%
  left_join(experiments, by = "fishery") %>%
  mutate(prepped_fishery = pmap(list(
    prepped_fishery = prepped_fishery,
    window = window,
    period = period,
    prop_years_lcomp_data = 1), subsample_data))

fisheries_sandbox$experiment <- 1:nrow(fisheries_sandbox)

fisheries_sandbox <- fisheries_sandbox %>% filter(economic_model == 0)

sfs <- safely(fit_scrooge)

doParallel::registerDoParallel(cores = n_cores)

foreach::getDoParWorkers()

set.seed(42)
fits <- foreach::foreach(i = 1:nrow(fisheries_sandbox)) %dopar% {
  out <- sfs(
    data = fisheries_sandbox$prepped_fishery[[i]]$scrooge_data,
    fish = fisheries_sandbox$prepped_fishery[[i]]$fish,
    fleet = fisheries_sandbox$prepped_fishery[[i]]$fleet,
    experiment = fisheries_sandbox$experiment[i],
    economic_model = 3,
    scrooge_file = "scrooge",
    iter = 6000,
    warmup = 4000,
    adapt_delta = 0.8,
    max_treedepth = 12,
    max_perc_change_f = 0.2,
    in_clouds = in_clouds,
    cloud_dir = cloud_dir,
    chains = 1,
    cv_effort = 0.5,
    q_guess = mean(possible_q),
    r0 = mean(possible_r0),
    sd_sigma_r = .4,
    cores = 1
  )

} # close fitting loop

#was getting stuck in the low thousand range
rstanarm::launch_shinystan(fits[[1]]$result)


fits[[1]]$result -> a

rstan::check_energy(a)

ppue <- tidybayes::spread_samples(a, ppue_hat[year])

true_ppue <- data_frame(year = 1:length(fisheries_sandbox$prepped_fishery[[1]]$scrooge_data$ppue_t
), ppue = fisheries_sandbox$prepped_fishery[[1]]$scrooge_data$ppue_t
)


ppue %>%
  left_join(true_ppue, by = "year") %>%
  ggplot() +
  geom_line(aes(year, ppue_hat, group = .iteration),alpha = 0.25) +
  geom_point(data = true_ppue, aes(year, ppue), color = "red")

ppue %>%
  left_join(true_ppue, by = "year") %>%
  ggplot() +
  geom_point(aes(ppue, ppue_hat, group = .iteration),alpha = 0.25)


wtf <- rstan::extract(a, "sigma_r", inc_warmup = TRUE, permuted = FALSE)

plot(wtf)

catch <- tidybayes::spread_samples(a, c_ta[year,age])

catch <- catch %>%
  group_by(year, age) %>%
  summarise(mean_catch = mean(c_ta))


 catch %>%
  ggplot(aes(age, mean_catch, fill = factor(year))) +
   geom_col(position = "dodge", show.legend = F)


rec_devs <- tidybayes::spread_samples(a, rec_dev_t[year])

rec_devs %>%
  ggplot(aes(year,rec_dev_t, group = .iteration)) +
  geom_line(alpha = 0.25)


lcomps <- tidybayes::spread_samples(a, p_lbin_sampled[year,lbin])

lcomps %>%
  group_by(year, lbin) %>%
  summarise(mean_n = mean(p_lbin_sampled)) %>%
  group_by(year) %>%
  mutate(mean_n = mean_n / sum(mean_n)) %>%
  ggplot(aes(lbin, mean_n, color = year, group = year)) +
  geom_line(show.legend = T)

psel_50_selectivity <- tidybayes::spread_samples(a, p_length_50_sel)


sel_50_selectivity <- tidybayes::spread_samples(a, length_50_sel)


selectivity <- tidybayes::spread_samples(a, mean_selectivity_at_age[age])

selectivity %>%
  ggplot(aes(age,mean_selectivity_at_age, group = .iteration)) +
  geom_line(alpha = 0.25)



selectivity %>%
  ggplot(aes(age,mean_selectivity_at_age, group = .iteration)) +
  geom_line(alpha = 0.25)


n_pop <- tidybayes::spread_samples(a, n_ta[year,age])

n_pop %>%
  group_by(year, age) %>%
  summarise(mean_n = mean(n_ta)) %>%
  group_by(year) %>%
  mutate(mean_n = mean_n / sum(mean_n)) %>%
  ggplot(aes(age, mean_n, color = year, group = year)) +
  geom_line(show.legend = T)

ssb_pop <- tidybayes::spread_samples(a, ssb_ta[year,age])


annual_n <- n_pop %>%
  group_by(year, .chain, .iteration) %>%
  summarise(n = sum(n_ta),
            recruits = n_ta[age == min(age)])

annual_ssb <- ssb_pop %>%
  group_by(year, .chain, .iteration) %>%
  summarise(ssb = sum(ssb_ta))

pop <- annual_n %>%
  left_join(annual_ssb, by = c("year", ".chain", ".iteration")) %>%
  group_by(.iteration) %>%
  arrange(year) %>%
  mutate(lr = lead(recruits)) %>%
  ungroup()

wtf <- pop %>%
  mutate(lr = lead)
  filter(ssb > 7.5e7, lead(recruits) < 22500000)

pop %>%
  ggplot(aes(year, n, group = .iteration)) +
  geom_line(alpha = 0.5) +
  geom_point()

pop %>%
  ggplot(aes(year, recruits, group = .iteration)) +
  geom_point(alpha = 0.5)

fisheries_sandbox$prepped_fishery[[1]]$simed_fishery %>%
  group_by(year) %>%
  summarise(recruits = numbers[age == min(age)],
            ssb = sum(ssb),
            tn = sum(numbers)) %>%
  ungroup() %>%
  ggplot(aes(ssb, lead(recruits))) +
  geom_point()

pairs(fits[[1]]$result, pars = c("sigma_r", "sigma_effort"))


pairs(fits[[1]]$result, pars = c("sigma_r", "p_length_50_sel"))

pairs(fits[[1]]$result, pars = c("p_length_50_sel","initial_f"))

pairs(fits[[1]]$result, pars = c("sigma_r", "p_length_50_sel","burn_f","initial_f"))


pairs(fits[[1]]$result, pars = c("sigma_r", "exp_rec_dev_t"))

pairs(fits[[1]]$result, pars = c("sigma_effort", "oc_effort_t"))

fisheries_sandbox <- fisheries_sandbox %>%
  slice(1)

fisheries_sandbox$scrooge_fit <- fits

scrooge_worked <-
  map(fisheries_sandbox$scrooge_fit, 'error') %>% map_lgl(is.null)


stan_worked <-
  map_lgl(map(fisheries_sandbox$scrooge_fit,"result"),
          ~ !(nrow(.x) %>% is.null()))


fisheries_sandbox <-
  fisheries_sandbox %>%
  filter(scrooge_worked) %>%  #,stan_worked, lime_worked, period != "middle")
  mutate(scrooge_fit = map(scrooge_fit, "result"))

scrooge_converged <-
  fisheries_sandbox$scrooge_fit %>%
  map_lgl( ~ class(rstan::get_logposterior(.x)) == "list")

fisheries_sandbox <-
  fisheries_sandbox %>%
  filter(scrooge_converged)

if (in_clouds == T) {
  loadfoo <- function(experiment, cloud_dir) {
    readRDS(glue::glue("{cloud_dir}/{experiment}"))
  }

  fisheries_sandbox <-
    fisheries_sandbox %>%
    mutate(scrooge_fit = map(scrooge_fit, safely(loadfoo), cloud_dir = cloud_dir))

}


processed_sandbox <-
  fisheries_sandbox %>%
  mutate(processed_scrooge = map2(
    scrooge_fit,
    map(prepped_fishery, "sampled_years"),
    process_scrooge,
    to_tidy = c("f_t", "n_tl","rec_dev_t")
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
  # mutate(lime_performance = pmap(
  #   list(observed = observed,
  #        predicted = processed_lime),
  #   judge_lime
  # )) %>%
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
