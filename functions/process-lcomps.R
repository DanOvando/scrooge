process_lcomps <- function(observed,predicted, sampled_years){

  sampled_years <- observed$sampled_years

  observed <- observed$length_comps


  observed <- observed %>%
    gather(length_bin, numbers,-year) %>%
    mutate(length_bin = as.numeric(length_bin)) %>%
    group_by(year) %>%
    mutate(p_numbers = numbers / sum(numbers)) %>%
    mutate(source = "observed")

  predicted <- predicted %>%
    gather(length_bin, numbers,-year,-iteration) %>%
    mutate(length_bin = as.numeric(str_replace_all(length_bin,"V",""))) %>%
    group_by(year, length_bin) %>%
    summarise(numbers = mean(numbers)) %>%
    mutate(p_numbers = numbers / sum(numbers)) %>%
    mutate(source = "predicted")


  length_fit_plot <-  ggplot() +
    geom_point(
      data = observed %>% filter(year %in% sampled_years),
      aes(length_bin, p_numbers),
      color = "red",
      alpha = 0.5
    ) +
    geom_line(data = predicted %>% filter(year %in% sampled_years),
              aes(length_bin, p_numbers), size = 1.25, alpha = 0.9) +
    facet_wrap( ~ year) +
    labs(x = "Lenth Bin (cm)", y = 'Proportional Numbers',
         caption = "Red points = observed, Black line = predicted") +
    theme(panel.spacing = unit(0.5, "lines"))

   out <- list(observed = observed,
               predicted = predicted,
               length_fit_plot = length_fit_plot)

}