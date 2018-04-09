judge_performance <-
  function(observed,
           predicted,
           observed_variable = f,
           predicted_variable = "f_t",
           group_level = year) {
    group_level = enquo(group_level)

    observed_variable = enquo(observed_variable)

    predicted_values <- purrr::pluck(predicted, predicted_variable)

    observed_values <- observed %>%
      select(!!group_level,!!observed_variable) %>%
      group_by(!!group_level) %>%
      summarise(observed = mean(!!observed_variable)) %>%
      mutate(year = year - min(year) + 1)


      comparison <- observed_values %>%
      left_join(predicted_values, by = "year") %>%
      ungroup() ## fix this later for full flexibility

    comparison_summary <- comparison %>%
      filter(is.na(value) == F) %>%
      summarise(rmse = sqrt(mean((value - observed) ^ 2)),
                bias = median((value - observed) / observed)) #%>%


    comparison_plot <- comparison %>%
      ggplot() +
      geom_boxplot(aes(year, value, group = year), alpha = 0.75) +
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