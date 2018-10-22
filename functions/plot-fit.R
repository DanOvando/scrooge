plot_fit <- function(performance,lcomp_years, plot_vars = c("f","recruitment_deviates")){

  test <- case_studies %>%
    slice(5)

  performance <- test$performance[[1]]

  lcomp_years <- test$prepped_fishery[[1]]$scrooge_data$length_comps_years


  if (plot_var == "all"){

    plot_vars <- unique(performance$others$variable)

  }

to_plot = performance$others %>%
  filter(variable %in% plot_vars)


ci90 <- to_plot %>%
  group_by(year, variable) %>%
  mutate(predicted_percentile = dplyr::percent_rank(predicted)) %>%
  filter(predicted_percentile < 0.95 & predicted_percentile > 0.05) %>%
  ungroup()

summaries <- to_plot %>%
  group_by(year, variable) %>%
  summarise(mp = median(predicted),
            observed = unique(observed))

# fill_vector <- rep("darkgrey", length (unique(summaries$year)))
#
# fill_vector[lcomp_years] <- "tomato"
#
# fill_vector <- fill_vector

ci90 %>%
  ggplot() +
  stat_bin2d(aes(year, predicted), binwidth = c(1, .01)) +
  geom_line(data = summaries,
            aes(year, mp),
            size = 2,
            color = "lightgrey") +
  geom_point(
    data = summaries,
    aes(year, observed),
    size = 4,
    fill = "darkgrey",
    shape = 21
  ) +
  scale_fill_viridis(option = "A",
                     guide = guide_colorbar(frame.colour = "black",
                                            barheight = 20)) +
  labs(x = "Year", y =  plot_var) +
  facet_wrap(~variable, scales = "free")

}