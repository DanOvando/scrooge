prepare_fishery <-
  function(sci_name,
           fleet_model,
           fleet_params,
           sigma_r,
           sigma_effort,
           price_cv,
           cost_cv,
           q_cv = 0,
           price_ac,
           cost_ac,
           q_ac = 0,
           price_slope = 0,
           cost_slope = 0,
           q_slope = 0,
           beta = 1.3,
           max_cp_ratio = 0.75,
           time_step = 1,
           price = 1,
           initial_f = 0.1,
           cost = 1,
           q = 1e-3,
           percnt_loo_selected = 0.25,
           sim_years = 15,
           burn_years = 10,
           linf_buffer = 1.5,
           num_patches = 1,
           r0 = 10000,
           sample_type = "catch",
           percent_sampled = 500,
           economic_model = 1,
           steepness = 0.8,
           profit_lags = 1,
           max_perc_change_f = 1,
           rec_ac = 0.25,
           query_price = T,
           bias = 0,
           obs_error = 0,
           est_msy = F,
           use_effort_data = 0,
           cv_effort = 0.25,
           cv_len = 0.1,
           seed = 42
           ) {

    initial_effort <- initial_f / q

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
        price_slope = price_slope,
        steepness = steepness,
        r0 = r0,
        rec_ac = rec_ac,
        cv_len = cv_len,
        density_movement_modifier = 0
      )

    if (fleet_model == "constant-catch"){
    fleet <- create_fleet(
      fish = fish,
      cost_cv =  cost_cv,
      cost_ac = cost_ac,
      cost_slope = cost_slope,
      q = q,
      q_cv = q_cv,
      q_ac = q_ac,
      q_slope = q_slope,
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
        q = q,
        q_cv = q_cv,
        q_ac = q_ac,
        cost_slope = cost_slope,
        q_slope = q_slope,
        fleet_model = fleet_model,
        initial_effort = initial_effort,
        cost = cost,
        sigma_effort = sigma_effort,
        length_50_sel = percnt_loo_selected * fish$linf
      )
    }

    if (fleet_model == "random-walk"){

      fleet <- create_fleet(
        fish = fish,
        cost_cv =  cost_cv,
        cost_ac = cost_ac,
        q = q,
        q_cv = q_cv,
        q_ac = q_ac,
        cost_slope = cost_slope,
        q_slope = q_slope,
        fleet_model = fleet_model,
        initial_effort = initial_effort,
        cost = cost,
        sigma_effort = sigma_effort,
        effort_ac = 0.75,
        length_50_sel = percnt_loo_selected * fish$linf
      )
    }


    if (fleet_model == "supplied-catch"){
      fleet <- create_fleet(
        fish = fish,
        cost_cv =  cost_cv,
        cost_ac = cost_ac,
        q = q,
        q_cv = q_cv,
        q_ac = q_ac,
        cost_slope = cost_slope,
        q_slope = q_slope,
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
        q = q,
        q_cv = q_cv,
        q_ac = q_ac,
        cost_slope = cost_slope,
        q_slope = q_slope,
        fleet_model = fleet_model,
        max_perc_change_f = max_perc_change_f,
        max_cp_ratio = max_cp_ratio,
        sigma_effort = sigma_effort,
        length_50_sel = percnt_loo_selected * fish$linf,
        initial_effort = initial_effort,
        profit_lags =  profit_lags,
        beta = beta
      )
    }

    set.seed(seed)
    sim <- sim_fishery(
      fish = fish,
      fleet = fleet,
      manager = create_manager(mpa_size = 0),
      num_patches = num_patches,
      sim_years = sim_years,
      burn_year = burn_years,
      time_step = fish$time_step,
      est_msy = FALSE
    )

    sim$year <- sim$year - min(sim$year) + 1

    scrooge_data <- prepare_scrooge_data(fish = fish,
                                     fleet = fleet,
                                     bias= bias,
                                     obs_error = obs_error,
                                     linf_buffer = linf_buffer,
                                     sim = sim,
                                     sample_type = "catch",
                                     percent_sampled = percent_sampled,
                                     cv_effort = cv_effort,
                                     economic_model = economic_model)



    out <- list(simed_fishery = sim,
                scrooge_data = scrooge_data,
                fish = fish,
                fleet = fleet)

    return(out)


  }
