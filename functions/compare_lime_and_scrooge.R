compare_lime_and_scrooge <- function(fitted_lime, processed_scrooge){

  comparison <- fitted_scrooge %>%
    mutate(year = (year - min(year) + 1)) %>%
    group_by(year) %>%
    summarise(scrooge_pred = mean(predicted),
              observed = mean(observed)) %>%
    ungroup() %>%
    left_join(fitted_lime %>% select(year, predicted) %>% rename(lime_pred = predicted), by = "year")

}