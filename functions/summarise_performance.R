summarise_performance <- function(experiment, prepped_fishery, cloud_dir) {

  fit <- readRDS(glue::glue("{cloud_dir}/{experiment}"))

  sampled_years <- prepped_fishery$sampled_years

  sampled_years <- data_frame(year =1:length(sampled_years),
                              sampled_year = sampled_years)

  f_t <- tidybayes::spread_samples(fit, f_t[year]) %>%
    ungroup() %>%
    left_join(sampled_years, by = "year") %>%
    mutate(year = sampled_year) %>%
    select(-sampled_year)

  true_f <- prepped_fishery$simed_fishery %>%
    group_by(year) %>%
    summarise(true_f = unique(f))

  f_pref <- f_t %>%
    left_join(true_f, by = "year") %>%
    mutate(resid = f_t - true_f) %>%
    rename(predicted = f_t,
           observed = true_f) %>%
    mutate(sq_er = resid^2,
           variable = "f")

  return(f_pref)

}