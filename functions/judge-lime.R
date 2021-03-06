judge_lime <-
  function(observed,
           predicted,
           observed_variable = f,
           predicted_variable = "f_t",
           group_level = year) {

    group_level = enquo(group_level)

    observed_variable = enquo(observed_variable)

    predicted_values <- predicted %>%
      filter(variable == predicted_variable)

    observed_values <- observed %>%
      dplyr::select(!!group_level,!!observed_variable) %>%
      group_by(!!group_level) %>%
      summarise(observed = mean(!!observed_variable)) %>%
      mutate(year = year - min(year) + 1)

      comparison <- observed_values %>%
      left_join(predicted_values, by = "year") %>%
      ungroup() ## fix this later for full flexibility

    comparison_summary <- comparison %>%
      filter(is.na(estimate) == F) %>%
      summarise(rmse = sqrt(mean((estimate - observed) ^ 2)),
                bias = median((estimate - observed) / observed)) #%>%


    comparison_plot <- comparison %>%
      ggplot() +
      geom_pointrange(aes(year, estimate, ymin = lower, ymax = upper), alpha = 0.75) +
      geom_line(aes(year, observed), color = "red") +
      geom_point(aes(year, observed), color = "red") +
      labs(
        caption = glue::glue(
          "RMSE = {prettyNum(comparison_summary$rmse, digits = 2)} ; Bias = {prettyNum(comparison_summary$bias, digits = 2)}"
        )
      )

    return(list(comparison_plot = comparison_plot,
                comparison_summary = comparison_summary))
  }