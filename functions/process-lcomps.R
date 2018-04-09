process_lcomps <- function(observed,predicted){
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
    geom_point(data = observed, aes(length_bin, p_numbers), color = "red")+
    geom_line(data = predicted,aes(length_bin, p_numbers)) +
    facet_wrap(~year) +
    labs(x = "Lenth Bin (cm)", y= 'Proportional Numbers',
         caption = "Red points = observed, Black line = predicted")

   out <- list(observed = observed,
               predicted = predicted,
               length_fit_plot = length_fit_plot)

}