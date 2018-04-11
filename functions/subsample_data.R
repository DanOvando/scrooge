subsample_data <- function(prepped_fishery, period = "beginning", window = 5){

  years <- prepped_fishery$length_comps$year %>% unique()

  if (period == "beginning"){


    sampled_years <- 1:window

  } else if (period == "middle"){

    sampled_years <-
      pmax(0, median(years) - floor(window / 2)):pmin(max(years), median(years) + floor(window/2))

  } else if (period == "end"){


    sampled_years <- (max(years) - window) : max(years)

  }

  prepped_fishery$sampled_years <- sampled_years

  prepped_fishery$scrooge_data$length_comps <- prepped_fishery$scrooge_data$length_comps %>%
    slice(sampled_years)

  prepped_fishery$scrooge_data$length_comps_years <- 1:length(sampled_years)

  prepped_fishery$scrooge_data$n_lcomps <- ifelse(length(sampled_years) > 1, length(sampled_years), array(1, dim = 1))

  prepped_fishery$scrooge_data$nt <- ifelse(length(sampled_years) > 1, length(sampled_years), array(1, dim = 1))

  prepped_fishery$scrooge_data$price_t <- prepped_fishery$scrooge_data$price_t %>%
    slice(sampled_years)

  prepped_fishery$scrooge_data$q_t <- prepped_fishery$scrooge_data$q_t %>%
 slice(sampled_years)

  prepped_fishery$scrooge_data$cost_t <- prepped_fishery$scrooge_data$cost_t %>%
    slice(sampled_years)

  return(prepped_fishery)


}