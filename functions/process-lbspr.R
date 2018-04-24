process_lbspr <- function(fit, sampled_years){



  out <- fit@Ests %>%
    as_data_frame() %>%
    set_names(tolower) %>%
    mutate(year = sampled_years)

}