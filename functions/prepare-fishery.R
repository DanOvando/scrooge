prepare_fishery <-
  function(sci_name,
           fleet_model,
           fleet_params,
           sigma_r,
           sigma_effort,
           price_cv,
           cost_cv,
           price_ac,
           cost_ac,
           q_cv = 0,
           q_ac = 0,
           time_step = 1,
           price = 1,
           initial_effort = 1000,
           cost = 1,
           percnt_loo_selected = 0.25,
           sim_years = 15,
           burn_years = 10,
           linf_buffer = 1.5,
           num_patches = 1,
           sample_type = "catch",
           percent_sampled = 1,
           economic_model = 1
           ) {



    fish <-
      create_fish(
        scientific_name = sci_name,
        query_fishlife = T,
        mat_mode = "length",
        time_step = time_step,
        sigma_r = sigma_r,
        price = price,
        price_cv = price_cv,
        price_ac = price_ac,
        r0 = 100000
      )

    if (fleet_model == "constant-catch"){
    fleet <- create_fleet(
      fish = fish,
      cost_cv =  cost_cv,
      cost_ac = cost_ac,
      q_cv = q_cv,
      q_ac = q_cv,
      fleet_model = fleet_model,
      target_catch = fleet_params$target_catch,
      cost = cost,
      sigma_effort = sigma_effort,
      length_50_sel = percnt_loo_selected * fish$linf
    )
    }

    if (fleet_model == "constant-effort"){
      fleet <- create_fleet(
        fish = fish,
        cost_cv =  cost_cv,
        cost_ac = cost_ac,
        q_cv = q_cv,
        q_ac = q_cv,
        fleet_model = fleet_model,
        initial_effort = fleet_params$initial_effort,
        cost = cost,
        sigma_effort = sigma_effort,
        length_50_sel = percnt_loo_selected * fish$linf
      )
    }

    if (fleet_model == "supplied-catch"){
      fleet <- create_fleet(
        fish = fish,
        cost_cv =  cost_cv,
        cost_ac = cost_ac,
        q_cv = q_cv,
        q_ac = q_cv,
        fleet_model = fleet_model,
        catches = fleet_params$catches,
        cost = cost,
        sigma_effort = sigma_effort,
        length_50_sel = percnt_loo_selected * fish$linf
      )

      fish$r0 <- max(fleet_params$catches)

      sim_years <- length(fleet$catches)

    }

    if (fleet_model == "open-access"){
      fleet <- create_fleet(
        fish = fish,
        cost_cv =  cost_cv,
        cost_ac = cost_ac,
        q_cv = q_cv,
        q_ac = q_cv,
        fleet_model = fleet_model,
        theta = fleet_params$theta,
        theta_tuner = fleet_params$theta_tuner,
        cost = cost,
        sigma_effort = sigma_effort,
        length_50_sel = percnt_loo_selected * fish$linf,
        initial_effort = fleet_params$initial_effort
      )
    }

    # fleet <-
    #   spasm::update_fleet(
    #     fleet = purrr::list_modify(
    #       fleet,
    #       cost = cost,
    #       price = price,
    #       sigma_effort = sigma_effort,
    #       length_50_sel = percnt_loo_selected * fish$linf
    #     ),
    #     fish = fish
    #   )

    sim <- spasm::sim_fishery(
      fish = fish,
      fleet = fleet,
      manager = create_manager(mpa_size = 0),
      num_patches = num_patches,
      sim_years = sim_years,
      burn_year = burn_years,
      time_step = fish$time_step
    )

    linf_buffer <- linf_buffer

    length_at_age_key <- generate_length_at_age_key(
      min_age = fish$min_age,
      max_age = fish$max_age,
      cv = fish$cv_len,
      linf = fish$linf,
      k = fish$vbk,
      t0 = fish$t0,
      time_step = fish$time_step,
      linf_buffer = linf_buffer
    ) %>%
      ungroup() %>%
      select(age, length_bin, p_bin) %>%
      spread(length_bin, p_bin) %>%
      select(-age)


    length_comps <- sim %>%
      select(year, patch, age, numbers, numbers_caught) %>%
      nest(-year,-patch, .key = n_at_age) %>%
      mutate(catch_length = map(
        n_at_age,
        ~ sample_lengths(
          n_at_age = .x,
          cv = fish$cv_len,
          k = fish$vbk,
          linf = fish$linf,
          t0 = fish$t0,
          sample_type = sample_type,
          percent_sampled = percent_sampled,
          time_step = fish$time_step,
          linf_buffer = linf_buffer
        )
      )) %>%
      select(year, catch_length) %>%
      unnest() %>%
      mutate(numbers = round(numbers)) %>%
      spread(length_bin, numbers) %>%
      mutate(year = year - min(year) + 1)

    price_and_cost_history <- sim %>%
      group_by(year) %>%
      summarise(price = unique(price),
                cost = unique(cost),
                q = fleet$q) %>%
      gather(variable, value,-year) %>%
      group_by(variable) %>%
      mutate(lag_value = lag(value, 1)) %>%
      ungroup() %>%
      mutate(lag_value = ifelse(is.na(lag_value), value, lag_value))

    price_t <- price_and_cost_history %>%
      filter(variable == "price")

    cost_t <- price_and_cost_history %>%
      filter(variable == "cost")

    q_t <- price_and_cost_history %>%
      filter(variable == "q")


    scrooge_data <- list(
      economic_model = economic_model,
      estimate_recruits = 1,
      length_comps = length_comps %>% select(-year),
      length_comps_years  = length_comps$year,
      price_t = price_t %>% select(value, lag_value),
      cost_t = cost_t %>% select(value, lag_value),
      q_t = q_t %>% select(value, lag_value),
      beta = 1.3,
      # base_effort = fish$m / mean(q_t$value),
      length_50_sel_guess = fish$linf / 2,
      delta_guess = 2,
      n_lcomps = nrow(length_comps),
      nt = length(length_comps$year),
      n_ages = fish$max_age + 1,
      n_lbins = ncol(length_at_age_key),
      ages = 1:(fish$max_age + 1),
      mean_length_at_age = fish$length_at_age,
      mean_weight_at_age = fish$weight_at_age,
      mean_maturity_at_age = fish$maturity_at_age,
      m = fish$m,
      h = fish$steepness,
      r0 = fish$r0,
      k = fish$vbk,
      loo = fish$linf,
      t0 = fish$t0,
      length_at_age_key = as.matrix(length_at_age_key)
    )

    out <- list(simed_fishery = sim,
                length_comps = length_comps,
                scrooge_data = scrooge_data,
                fish = fish,
                fleet = fleet)

    return(out)


  }