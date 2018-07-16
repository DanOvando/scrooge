plot_simmed_fishery <- function(prepped_fishery){

  sim <- prepped_fishery$simed_fishery

  current_theme <- theme_get()

  theme_set(theme_classic())


  biomass_plot <- sim %>%
    group_by(year) %>%
    summarise(biomass = sum(biomass)) %>%
    ggplot(aes(year, biomass)) +
    geom_line() +
    geom_point()


 catch_plot <- sim %>%
    group_by(year) %>%
    summarise(catch = sum(biomass_caught)) %>%
    ggplot(aes(year, catch)) +
    geom_line() +
    geom_point()


 f_plot <- sim %>%
   group_by(year) %>%
   summarise(f = unique(f)) %>%
   ggplot(aes(year, f)) +
   geom_line() +
   geom_point()


 length_plot <- prepped_fishery$scrooge_data$length_comps %>%
   gather(length_bin, numbers, -year) %>%
   mutate(length_bin = as.numeric(length_bin)) %>%
   ungroup() %>%
   filter(year %in% c(min(year),max(year))) %>%
   group_by(year) %>%
   mutate(numbers = numbers / sum(numbers)) %>%
   ungroup() %>%
   ggplot(aes(length_bin, numbers, fill = factor(year))) +
   geom_col(position = "dodge") +
   scale_fill_discrete(name = 'Year')

 outplot <- biomass_plot + catch_plot + f_plot + length_plot + plot_layout(ncol = 2, nrow = 2)


  theme_set(current_theme)

  return(outplot)

}