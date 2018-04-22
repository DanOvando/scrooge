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
           b_v_bmsy_oa = 0.5,
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
           economic_model = 1,
           steepness = 0.8,
           profit_lags = 1,
           rec_ac = 0.25,
           query_price = T,
           bias = 0,
           obs_error = 0,
           tune_costs = F,
           est_msy = F,
           mey_buffer = 2,
           use_effort_data = 0
           ) {


    if (query_price == T) {
      prices <-
        read_csv(file = here::here("data", "Exvessel Price Database.csv"))

      this_price <- prices %>%
        filter(scientific_name %in% sci_name, Year > 2000)

      if (nrow(this_price == 0)) {
        price = 4

      } else{
        price <- mean(this_price$exvessel, na.rm = T) * .001
      }
    }


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
        steepness = steepness,
        r0 = 1000,
        rec_ac = rec_ac
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
      length_50_sel = percnt_loo_selected * fish$linf,
      profit_lags =  profit_lags
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
        cost = cost,
        sigma_effort = sigma_effort,
        length_50_sel = percnt_loo_selected * fish$linf,
        initial_effort = fleet_params$initial_effort,
        profit_lags =  profit_lags,
        mey_buffer = mey_buffer
      )

      tune_costs <- T

      est_msy <- T
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
    #


    sim <- sim_fishery(
      fish = fish,
      fleet = fleet,
      manager = create_manager(mpa_size = 0),
      num_patches = num_patches,
      sim_years = sim_years,
      burn_year = burn_years,
      time_step = fish$time_step,
      est_msy = est_msy,
      tune_costs = tune_costs,
      b_v_bmsy_oa = b_v_bmsy_oa
    )

    linf_buffer <- linf_buffer

    # add in observation error
    true_min <- fish$min_age

    true_max <- fish$max_age

    time_step <- fish$time_step

    fish <- map_if(fish, is.numeric, ~.x * exp(rnorm(1, rnorm(1,0, bias), obs_error)))

    fish$min_age <- true_min

    fish$max_age <- true_max

    fish$time_step <- time_step

    fish$length_at_age <-   fish$linf * (1 - exp(- fish$vbk * (seq(fish$min_age,fish$max_age, by = fish$time_step) -  fish$t0)))

    fish$weight_at_age <-  fish$weight_a * fish$length_at_age ^ fish$weight_b

    fish$maturity_at_age <-
      ((1 / (1 + exp(-log(
        19
      ) * ((seq(fish$min_age, fish$max_age, by = fish$time_step) - fish$age_50_mature) / pmax(0.01,fish$age_95_mature - fish$age_50_mature)
      )))))

    fleet <- map_if(fleet, is.numeric, ~.x * exp(rnorm(1, rnorm(1,0, bias), obs_error)))

    fleet$length_95_sel <- fleet$length_50_sel + fleet$delta

    fleet$sel_at_age <-
      ((1 / (1 + exp(-log(
        19
      ) * ((fish$length_at_age - fleet$length_50_sel) / pmax(0.01,fleet$length_95_sel - fleet$length_50_sel)
      )))))


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
        ~ spasm::sample_lengths(
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
                q = 0.1) %>%
      gather(variable, value,-year) %>%
      ungroup() %>%
      mutate(value = value * exp(rnorm(nrow(.), rnorm(nrow(.),0, bias), obs_error))) %>%
      group_by(variable) %>%
      mutate(centered_value = value / mean(value)) %>%
      ungroup()

    price_t <- price_and_cost_history %>%
      filter(variable == "price")

    cost_t <- price_and_cost_history %>%
      filter(variable == "cost")

    q_t <- price_and_cost_history %>%
      filter(variable == "q")

    effort_t <- sim %>%
      group_by(year) %>%
      summarise(effort = mean(effort)) %>%
      ungroup() %>%
      mutate(relative_effort = effort/max(effort))

    scrooge_data <- list(
      economic_model = economic_model,
      estimate_recruits = 1,
      length_comps = length_comps %>% select(-year),
      relative_effort = effort_t$relative_effort,
      length_comps_years  = length_comps$year,
      price_t = price_t$centered_value,
      cost_t = cost_t$centered_value,
      q_t = q_t$centered_value,
      beta = 1.3,
      # base_effort = fish$m / mean(q_t$value),
      length_50_sel_guess = fleet$length_50_sel * exp(rnorm(1, rnorm(1,0, bias), obs_error)),
      delta_guess = 2,
      n_lcomps = nrow(length_comps),
      nt = length(length_comps$year),
      n_ages = fish$max_age + 1,
      n_lbins = ncol(length_at_age_key),
      ages = 1:(fish$max_age + 1),
      m = fish$m,
      h = fish$steepness,
      r0 = fish$r0,
      k = fish$vbk,
      loo = fish$linf,
      t0 = fish$t0,
      length_at_age_key = as.matrix(length_at_age_key),
      mean_length_at_age = fish$length_at_age,
      mean_weight_at_age = fish$weight_at_age,
      mean_maturity_at_age = fish$maturity_at_age,
      use_effort_data = use_effort_data
    )

    out <- list(simed_fishery = sim,
                length_comps = length_comps,
                scrooge_data = scrooge_data,
                fish = fish,
                fleet = fleet)

    return(out)


  }