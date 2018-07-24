plot_fits <- function(prepped_fishery, fit){

# check <- case_studies %>%
#   filter(economic_model == 3, case_study == "realistic",
#          period == "end", window == "10")
#

  check <- case_studies %>% slice(5)

  prepped_fishery <- check$prepped_fishery[[1]]

  fit <- check$scrooge_fit[[1]]$result

  sampled_years <- prepped_fishery$sampled_years

  sampled_years <- data_frame(year =1:length(sampled_years),
                              sampled_year = sampled_years)

  rec_devs <- tidybayes::spread_samples(fit, rec_dev_t[year])

  true_recdevs <- prepped_fishery$simed_fishery %>%
    group_by(year) %>%
    summarise(true_rec_dev = exp(unique(rec_dev))) %>%
    ungroup() %>%
    mutate(year = 1:nrow(.))

  rec_devs %>%
    left_join(true_recdevs, by = "year") %>%
    ggplot(aes(year,rec_dev_t, group = .iteration)) +
    geom_line(alpha = 0.25) +
    geom_point(aes(year,true_rec_dev), color = "red") +
    geom_hline(aes(yintercept = 1), color = "red")


  f_t <- tidybayes::spread_samples(fit, f_t[year]) %>%
    ungroup() %>%
    left_join(sampled_years, by = "year") %>%
    mutate(year = sampled_year) %>%
    select(-sampled_year)

  true_f <- check$prepped_fishery[[1]]$simed_fishery %>%
    group_by(year) %>%
    summarise(true_f = unique(f)) %>%
    mutate(year = year - min(year) + 1)

  f_plot <- f_t %>%
    left_join(true_f, by = "year") %>%
    group_by(year) %>%
    summarise(
      lower_90 = quantile(f_t, 0.05),
      upper_90 = quantile(f_t, 0.95),
      lower_50 = quantile(f_t, 0.25),
      upper_50 = quantile(f_t, 0.75),
      mean_f = mean(f_t),
              true_f = mean(true_f)) %>%
    ggplot() +
    geom_ribbon(aes(year, ymin = lower_90, ymax = upper_90), fill = "lightgrey") +
    geom_ribbon(aes(year, ymin = lower_50, ymax = upper_50), fill = "darkgrey") +
    geom_line(aes(year,mean_f), color = "steelblue") +
    geom_point(aes(year, true_f), fill = "tomato", size = 4, shape = 21) +
    labs(y = "F", x = "Year")

  ppue <- tidybayes::spread_samples(fit, ppue_hat[year]) %>%
    ungroup() %>%
    left_join(sampled_years, by = "year") %>%
    mutate(year = sampled_year) %>%
    select(-sampled_year)

  true_ppue <-
    data_frame(year = sampled_years$sampled_year,
               true_ppue =  prepped_fishery$scrooge_data$ppue_t) %>%
    ungroup()


  ppue_plot <- ppue %>%
    left_join(true_ppue, by = "year") %>%
    mutate(resid = ppue_hat - true_ppue) %>%
    group_by(year) %>%
    summarise(
      lower_90 = quantile(ppue_hat, 0.05),
      upper_90 = quantile(ppue_hat, 0.95),
      lower_50 = quantile(ppue_hat, 0.25),
      upper_50 = quantile(ppue_hat, 0.75),
      mean_val = mean(ppue_hat),
      true_val = mean(true_ppue)) %>%
    ggplot() +
    geom_ribbon(aes(year, ymin = lower_90, ymax = upper_90), fill = "lightgrey") +
    geom_ribbon(aes(year, ymin = lower_50, ymax = upper_50), fill = "darkgrey") +
    geom_line(aes(year,mean_val), color = "steelblue") +
    geom_point(aes(year, true_val), fill = "tomato", size = 4, shape = 21) +
    labs(y = "PPUE", x = "Year")

  observed_lcomps <- prepped_fishery$scrooge_data$length_comps %>%
    as_data_frame() %>%
    mutate(year = sampled_years$sampled_year) %>%
    gather(lbin, numbers, -year) %>%
    mutate(lbin = as.numeric(lbin)) %>%
    group_by(year) %>%
    mutate(numbers = numbers / sum(numbers))

  pp_n_tl <- tidybayes::spread_samples(fit, n_tl[year,length_bin]) %>%
    ungroup() %>%
    left_join(sampled_years, by = "year") %>%
    mutate(year = sampled_year) %>%
    select(-sampled_year)

  pp_length_plot <- pp_n_tl %>%
    group_by(year, .chain,.iteration) %>%
    mutate(p_n_tl = n_tl / sum(n_tl)) %>%
    group_by(year, .chain, length_bin) %>%
    summarise(lower_90 = quantile(p_n_tl,0.05),
              upper_90 = quantile(p_n_tl,0.95)) %>%
    ggplot() +
    geom_ribbon(aes(x = length_bin, ymin = lower_90, ymax = upper_90), fill = "lightgrey") +
    facet_wrap(~year) +
    theme_minimal() +
    geom_point(data = observed_lcomps, aes(lbin, numbers), size = .5, alpha = 0.5, color = "red")


  lcomps <- tidybayes::spread_samples(fit, p_lbin_sampled[year,lbin]) %>%
    ungroup() %>%
    left_join(sampled_years, by = "year") %>%
    mutate(year = sampled_year) %>%
    select(-sampled_year)


 length_plot <-  lcomps %>%
    group_by(year, lbin) %>%
    summarise(mean_n = mean(p_lbin_sampled)) %>%
    group_by(year) %>%
    mutate(mean_n = mean_n / sum(mean_n)) %>%
    ggplot() +
    geom_line(aes(lbin, mean_n, color = year, group = year),show.legend = F) +
    geom_point(data = observed_lcomps, aes(lbin, numbers), size = .5, alpha = 0.5) +
    facet_wrap(~year) +
    theme_classic()

}