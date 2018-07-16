subsample_data <-
  function(prepped_fishery,
           period = "beginning",
           window = 5,
           prop_years_lcomp_data = 1) {


    if (prop_years_lcomp_data > 1){
      stop("prop_years_lcomp_data must be <= 1")
    }

    years <- prepped_fishery$scrooge_data$length_comps$year %>% unique()

    if (period == "beginning") {
      sampled_years <- 1:window

    } else if (period == "middle") {
      sampled_years <-
        floor(pmax(0, median(years) - floor(window / 2)):pmin(max(years), median(years) + floor(window /
                                                                                                  2)))

    } else if (period == "end") {
      sampled_years <- (max(years) - window):max(years)

    }

    length_years <-
      sampled_years[-(1:(
        length(sampled_years) - ceiling(length(sampled_years) * prop_years_lcomp_data)
      ))]

    prepped_fishery$sampled_years <- sampled_years

    prepped_fishery$scrooge_data$perc_change_effort <-
      prepped_fishery$scrooge_data$perc_change_effort[sampled_years]


    prepped_fishery$scrooge_data$observed_effort <-
      prepped_fishery$scrooge_data$observed_effort[sampled_years]

    prepped_fishery$scrooge_data$length_comps <-
      prepped_fishery$scrooge_data$length_comps %>%
      slice(length_years)

    prepped_fishery$scrooge_data$length_comps_years <-
      which(sampled_years %in% length_years)

    prepped_fishery$scrooge_data$n_lcomps <-
      ifelse(length(length_years) > 1,
             length(length_years),
             array(1, dim = 1))

    prepped_fishery$scrooge_data$nt <-
      ifelse(length(sampled_years) > 1,
             length(sampled_years),
             array(1, dim = 1))

    prepped_fishery$scrooge_data$price_t <-
      prepped_fishery$scrooge_data$price_t[sampled_years]

    prepped_fishery$scrooge_data$q_t <-
      prepped_fishery$scrooge_data$q_t[sampled_years]

    prepped_fishery$scrooge_data$cost_t <-
      prepped_fishery$scrooge_data$cost_t[sampled_years]


    return(prepped_fishery)


  }