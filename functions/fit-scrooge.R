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
           model_type = "scrooge",
           max_treedepth = 10,
           max_f_v_fmsy_increase = 0.5,
           in_clouds = F,
           cloud_dir = "results/scrooge_results",
           price = 1,
           q = 0.01,
           r0 = 100,
           b_v_bmsy_oa = 0.5,
           init_f_v_m = 0.8,
           use_effort_data = 1
           ) {

    data$use_effort_data <- use_effort_data

    data$economic_model <- economic_model

    data$sigma_r_guess <- 0.4

    data$r0 <- 100

    data$f_init_guess <- fish$m * init_f_v_m

    p_response = max_f_v_fmsy_increase

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

    # fs <- seq(0,2, by = 0.02)
    #
    # test <- NA
    # for (i in seq_along(fs)){
    #   test[i] <- est_msy(fs[i], data = data, time = 100, use = "fit")
    #
    # }

    e_msy <- (f_msy$par / q)

    msy <- -f_msy$objective

    cost <-
      nlminb(
        0.2,
        tune_costs,
        data = data,
        time = 200,
        lower = 0,
        p_response = p_response,
        b_v_bmsy_oa = b_v_bmsy_oa,
        msy = msy,
        e_msy = e_msy,
        b_msy = b_msy,
        price = price,
        q = q
      )

    if (cost$convergence != 0){
      stop("check cost convergence")
    }

    # check <- tune_costs(cost$par, data = data,
    #                     time = 1000,
    #                     p_response = p_response,
    #                     b_v_bmsy_target = 0.5,
    #                     msy = msy,
    #                     e_msy = e_msy,
    #                     b_msy = b_msy,
    #                     price = price,
    #                     q = q,
    #                     use = "blah")

    data$cost_t <- cost$par * data$cost_t

    data$price_t <- price * data$price_t

    data$q_t <- q * data$q_t

    p_msy <- price*msy - cost$par * e_msy ^ data$beta

    data$p_msy <- p_msy

    data$e_msy <- e_msy

    data$p_response_guess <- p_response

    inits <-
      map(
        1:chains,
        ~ list(
          p_length_50_sel = 0.25 * exp(rnorm(1, 0, .1)),
          p_response = p_response * exp(rnorm(1, 0, .1))
        )
      )


    # inits <-
    #   map(
    #     1:chains,
    #     ~ list(
    #       base_effort = 1 *  exp(rnorm(1, 0, .1)),
    #       p_length_50_sel = 0.25 * exp(rnorm(1, 0, .1)),
    #       p_response = p_response * exp(rnorm(1, 0, .1))
    #     )
    #   )



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
        init = inits
      )



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
