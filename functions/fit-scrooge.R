fit_scrooge <-
  function(data,
           fish,
           fleet,
           experiment,
           scrooge_file = "scrooge",
           chains = 1,
           refresh = 25,
           cores = 1,
           iter = 1000,
           warmup = 2000,
           adapt_delta = 0.8,
           economic_model = 1,
           likelihood_model = 0,
           model_type = "scrooge",
           max_treedepth = 10,
           max_perc_change_f = 0.5,
           in_clouds = F,
           cloud_dir = "results/scrooge_results",
           price = 2,
           q_guess = 0.1,
           r0 = 1000,
           init_f_v_m = 0.8,
           cv_effort = 1.6,
           cp_guess = 0.75,
           h = 0.8,
           sd_sigma_r = 0.001,
           seed = 42,
           n_burn = 50,
           sigma_effort = 0.2
           ) {

    data$age_sel <- floor((log(1-pmin(data$length_50_sel_guess, data$loo*.99)/data$loo)/-data$k)+data$t0)

    data$bin_mids <- as.numeric(colnames(data$length_at_age_key))

    data$bin_mids <- data$bin_mids + ((data$bin_mids[2] - data$bin_mids[1])/2)

    data$n_burn <-  n_burn

    data$sigma_effort <- sigma_effort

    data$length_comps <- data$length_comps %>%
      select(-year)

    data$economic_model <- economic_model

    data$sigma_r_guess <- 0

    data$r0 <- r0

    data$h <- fish$steepness

    data$f_init_guess <- fish$m * init_f_v_m

    data$sd_sigma_r <- sd_sigma_r

    if (economic_model == 0){

      data$q_t <- rep(mean(data$q_t), length(data$q_t))

    }

    if (is.na(cv_effort) == F){
    data$cv_effort <-  cv_effort
    }

    f_msy <-
      nlminb(
        data$m,
        est_msy,
        data = data,
        time = 100,
        lower = 0,
        upper = 2,
        use = "fit"
      )

    if (f_msy$convergence != 0){
      stop("check msy convergence")
    }

    b_msy <- est_msy(f_msy$par,
                     data = data,
                     time = 100,
                     use = "else")

    data$q_t <- q_guess * data$q_t

    n0_at_age <-
      r0  * exp(-fish$m * seq(fish$min_age, fish$max_age, fish$time_step))

    n0_at_age[fish$max_age + 1] <-
      n0_at_age[fish$max_age + 1] / (1 - exp(-fish$m))

    b0_at_age <- n0_at_age * fish$weight_at_age

    hyp_f <- 3*fish$m #hypothetical f

    hyp_effort <- hyp_f / min(data$q_t)

    hyp_f_at_age <- hyp_f * fleet$sel_at_age

    hyp_b0_catch <- sum((hyp_f_at_age / (hyp_f_at_age + fish$m))  * b0_at_age * (1 - exp(-(hyp_f_at_age + fish$m))))

    b0_revenue <- max(data$price_t) * hyp_b0_catch

    hyp_profits_guess <- b0_revenue * (1 - cp_guess)

    cost_guess <-
      (b0_revenue - hyp_profits_guess) / hyp_effort ^ fleet$beta

    data$max_cost_guess <- cost_guess

    data$relative_cost_t <-  data$cost_t

    data$p_response_guess <- (max_perc_change_f * hyp_effort) / (hyp_profits_guess / hyp_effort)

    # delta_effort <- data$p_response_guess * (hyp_profits_guess / hyp_effort)

    # data$ppue_t <- data$ppue_t

  # if (max(data$ppue_t) > 0){
  # data$ppue_t <- data$ppue_t/max(data$ppue_t)
  # } else{
  #
  #   data$ppue_t <- data$ppue_t/mean(data$ppue_t)
  #
  # }

  inits <-
    map(
      1:chains,
      ~ list(
        sigma_r = jitter(.1),
        p_length_50_sel = jitter(0.25),
        initial_f = jitter(data$m),
        log_effort_dev_t = rep(0, data$nt)
      )
    )

  fit <-
      rstan::stan(
        file = here::here("src", paste0(scrooge_file, ".stan")),
        data = data,
        chains = chains,
        refresh = refresh,
        cores = cores,
        iter = iter,
        warmup = warmup,
        control = list(adapt_delta = adapt_delta,
                       max_treedepth = max_treedepth),
        init = inits,
        seed = seed
      )
    # init = inits


    # clean up old DLLS per https://github.com/stan-dev/rstan/issues/448

    loaded_dlls = getLoadedDLLs()

    loaded_dlls = loaded_dlls[str_detect(names(loaded_dlls), '^file')]
    if (length(loaded_dlls) > 25) {
      for (dll in head(loaded_dlls,-15)) {
        # message("Unloading DLL ", dll[['name']], ": ", dll[['path']])
        dyn.unload(dll[['path']])
      }
    }



    if (in_clouds == T) {
      filename <- glue::glue("experiment_{experiment}.rds")

      saveRDS(fit, file = glue::glue("{cloud_dir}/{filename}"))

      rm(fit)

      fit <-  filename

    }

    out <- fit


    return(out)

  }
