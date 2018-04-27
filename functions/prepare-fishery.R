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
           r0 = 10000,
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
           use_effort_data = 0,
           cv_effort = 0.25,
           cv_len = 0.1
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
        r0 = r0,
        rec_ac = rec_ac,
        cv_len = cv_len
      )

    if (fleet_model == "constant-catch"){
    fleet <- create_fleet(
      fish = fish,
      cost_cv =  cost_cv,
      cost_ac = cost_ac,
      q_cv = q_cv,
      q_ac = q_ac,
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
        q_ac = q_ac,
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
        q_ac = q_ac,
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
        q_ac = q_ac,
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
                length_comps = length_comps,
                scrooge_data = scrooge_data,
                fish = fish,
                fleet = fleet)

    return(out)


  }