assess_fits <- function(prepped_fishery, fit){

  sampled_years <- prepped_fishery$sampled_years

  sampled_years <- data_frame(year =1:length(sampled_years),
                              sampled_year = sampled_years)

  rec_devs <- tidybayes::spread_samples(fit, rec_dev_t[year])

  true_recdevs <- prepped_fishery$simed_fishery %>%
    group_by(year) %>%
    summarise(true_rec_dev = exp(unique(rec_dev))) %>%
    ungroup() %>%
    mutate(year = 1:nrow(.))

  rec_pref <- rec_devs %>%
    left_join(true_recdevs, by = "year") %>%
    mutate(resid = rec_dev_t - true_rec_dev) %>%
    mutate(sq_er = resid^2,
           variable = "recruitment_deviates") %>%
    rename(predicted = rec_dev_t,
           observed = true_rec_dev)

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

# ppue

ppue_hat_t <- tidybayes::spread_samples(fit, ppue_hat_t[year]) %>%
    ungroup() %>%
    left_join(sampled_years, by = "year") %>%
    mutate(year = sampled_year) %>%
    select(-sampled_year)

  true_ppue <- prepped_fishery$simed_fishery %>%
    group_by(year) %>%
    summarise(profits = sum(profits),
              effort = unique(effort)) %>%
    ungroup() %>%
    mutate(ppue = profits / effort)

  ppue_pref <- ppue_hat_t %>%
    left_join(true_ppue, by = "year") %>%
    mutate(resid = ppue_hat_t - ppue) %>%
    rename(predicted = ppue_hat_t,
           observed = ppue) %>%
    mutate(sq_er = resid^2,
           variable = "ppue") %>%
    mutate(predicted = predicted / max(predicted),
           observed = observed / max(observed))

  out <- f_pref %>%
    bind_rows(rec_pref) %>%
    bind_rows(ppue_pref)


return(out)

}