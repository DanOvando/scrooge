process_lime <-
  function(fit,
           sampled_years,
           vars = c("log_F_t_input")) {
    report_names <- fit$Sdreport %>%
      summary() %>%
      rownames()


    sd_report <- fit$Sdreport %>%
      summary() %>%
      as_data_frame() %>%
      mutate(variable = report_names) %>%
      filter(variable %in% vars) %>%
      mutate(year = sampled_years) %>%
      mutate(
        estimate = exp(Estimate),
        lower = exp(Estimate - 1.96 * `Std. Error`),
        upper = exp(Estimate + 1.96 * `Std. Error`)
      )
    out <- sd_report %>%
      select(variable, year, estimate, lower, upper) %>%
      mutate(model = "lime") %>%
      mutate(variable = "f_t")

  }