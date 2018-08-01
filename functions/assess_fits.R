assess_fits <- function(prepped_fishery, fit){

  sampled_years <- prepped_fishery$sampled_years

  sampled_years <- data_frame(year =1:length(sampled_years),
                              sampled_year = sampled_years)



  length_comp_sampled_years <- data_frame(year = 1:nrow(prepped_fishery$scrooge_data$length_comps),
                                          sampled_year = prepped_fishery$scrooge_data$length_comps_years)

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

ppue_hat <- tidybayes::spread_samples(fit, ppue_hat[year]) %>%
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

  ppue_pref <- ppue_hat %>%
    left_join(true_ppue, by = "year") %>%
    mutate(resid = ppue_hat - ppue) %>%
    rename(predicted = ppue_hat,
           observed = ppue) %>%
    mutate(sq_er = resid^2,
           variable = "ppue") %>%
    mutate(predicted = predicted / max(predicted),
           observed = observed / max(observed))

  #  percent change in effort

  perc_change_effort_hat <- tidybayes::spread_samples(fit, perc_change_effort_hat[year]) %>%
    ungroup() %>%
    left_join(sampled_years, by = "year") %>%
    mutate(year = sampled_year) %>%
    select(-sampled_year)

  true_perc_change_effort_hat <- prepped_fishery$simed_fishery %>%
    group_by(year) %>%
    summarise(effort = unique(effort)) %>%
    ungroup() %>%
    arrange(year) %>%
    mutate(perc_change_effort = lead(effort) / effort)

  perc_change_effort_pref <- perc_change_effort_hat %>%
    left_join(true_perc_change_effort_hat, by = "year") %>%
    mutate(resid = perc_change_effort_hat - perc_change_effort) %>%
    rename(predicted = perc_change_effort_hat,
           observed = perc_change_effort) %>%
    mutate(sq_er = resid^2,
           variable = "percent_change_effort") %>%
    mutate(predicted = predicted / max(predicted),
           observed = observed / max(observed))

  # length comps

  observed_lcomps <- prepped_fishery$scrooge_data$length_comps %>%
    as_data_frame() %>%
    mutate(year = prepped_fishery$scrooge_data$length_comps_years) %>%
    gather(lbin, numbers, -year) %>%
    mutate(lbin = as.numeric(lbin)) %>%
    group_by(year) %>%
    mutate(numbers = numbers / sum(numbers))

  pp_length_comps <- tidybayes::spread_samples(fit, n_tl[year,length_bin]) %>%
    ungroup() %>%
    left_join(observed_lcomps, by = c("year", "length_bin" = "lbin")) %>%
    mutate(source = "posterior_predictive")  %>%
    rename(observed = numbers, predicted = n_tl)

  fitted_length_comps <- tidybayes::spread_samples(fit, p_lbin_sampled[year,lbin]) %>%
    ungroup() %>%
    left_join(observed_lcomps, by = c("year", "lbin" = "lbin")) %>%
    mutate(source = "fitted")  %>%
    rename(observed = numbers, predicted = p_lbin_sampled) %>%
    rename(length_bin = lbin)

  length_comps <- pp_length_comps %>%
    bind_rows(fitted_length_comps)


  out <- list(length_comps = length_comps,
              others =     f_pref %>%
                bind_rows(rec_pref) %>%
                bind_rows(ppue_pref) %>%
                bind_rows(perc_change_effort_pref))

return(out)

}