plot_case_studies <- function(cs1,
         cs2,
         performance_stats,
         case_studies,
         cs_name = "realistic",
         lcomp_years = "max",
         plot_variable = "f",
         plot_period = "middle",
         modelo_name) {

  if (lcomp_years == "max"){

    lyears <- max(performance_stats$prop_years_lcomp_data)
  } else if (lcomp_years == "min"){

    lyears <- min(performance_stats$prop_years_lcomp_data)

  }


cs_summary <- performance_stats %>%
  ungroup() %>%
  filter(case_study == case_study,
         prop_years_lcomp_data == lyears ,
         variable == plot_variable,
         model == cs1 | model == cs2)

# cs_summary$model <- fct_recode(as.factor(cs_summary$model),  "Lengths Only" = cs1,
#                                 "Lengths + Economics" = cs2 )
cs_summ_plot <- cs_summary %>%
  gather(metric, value, recent_rmse, recent_median_bias) %>%
  ggplot(aes(value, fill = model)) +
  geom_vline(aes(xintercept = 0), linetype = 2) +
  geom_density_ridges(aes(x = value, y = metric, fill = model),
                      alpha = 0.3,
                      color = "darkgrey",
                      show.legend = F) +
  scale_y_discrete(labels = c("bias", "rmse"), position = "right", name = element_blank())+
  theme(
    plot.margin = unit(c(0, 0, 0, 1), units = "lines"),
    axis.title.y = element_blank(),
    axis.title.x = element_blank()
  ) +
  theme(legend.title = element_blank(), axis.text.x = element_text(size = 10)) +
  ggsci::scale_fill_aaas() +
  labs(caption = "Calculated over last 5 years", title = "B")


lcomp_years <- case_studies %>%
  filter(case_study == cs_name,
         prop_years_lcomp_data == lyears,
         period == plot_period)

lcomp_years <- lcomp_years$prepped_fishery[[1]]$scrooge_data$length_comps_years


cs <- perf_summaries %>%
  ungroup() %>%
  filter(case_study == cs_name,
         prop_years_lcomp_data == lyears,
         model == cs1 | model == cs2)

fill_vector <- rep("transparent", length (unique(cs$year[cs$variable == "f"])))

fill_vector[lcomp_years] <- "tomato"

fill_vector <- c(fill_vector, fill_vector)

cs_plot <- cs %>%
  filter(variable == plot_variable) %>%
  mutate(year = year - min(year) + 1) %>%
  mutate(lcomp_year = year %in% lcomp_years) %>%
  arrange(experiment) %>%
  ggplot() +
  geom_point(aes(year, observed), size = 3, alpha = 0.75, shape = 21, fill =fill_vector) +
  geom_ribbon(aes(year, ymin = lower_90, ymax = upper_90, fill = model), alpha = 0.2, show.legend = FALSE) +
  geom_line(aes(year,mean_predicted, color = model), show.legend = FALSE) +
  labs(title = "A",x = "Year", y = "Fishing Mortality", caption = "Points = Data, Shaded = Predictions") +
  ggsci::scale_color_aaas() +
  ggsci::scale_fill_aaas() +
  facet_wrap(~model, labeller = labeller(model = modelo_name)) +
  theme(plot.margin = unit(c(0,0,0,0), units = "lines"),
        panel.spacing = unit(0.1, units = "lines"),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 14))


out <- cs_plot + cs_summ_plot + plot_layout(nrow = 1, ncol = 2, widths = c(3,1))

}
