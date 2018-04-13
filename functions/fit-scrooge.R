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
           price = 0.1) {

    data$sigma_r_guess <- 0.4

    data$r0 <- 100

    p_expansion = max_f_v_fmsy_increase

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


    e_msy <- (f_msy$par / mean(data$q_t$value))

    msy <- -f_msy$objective

    cost <-
      nlminb(
        0.2,
        tune_costs,
        data = data,
        time = 1000,
        lower = 0,
        p_expansion = p_expansion,
        b_v_bmsy_target = 0.5,
        msy = msy,
        e_msy = e_msy,
        b_msy = b_msy,
        price = price
      )


    # check <- tune_costs(cost$par, data = data,
    #                     time = 10000,
    #                     p_expansion = p_expansion,
    #                     b_v_bmsy_target = 0.5,
    #                     msy = msy,
    #                     e_msy = e_msy,
    #                     b_msy = b_msy,
    #                     price = price,
    #                     use = "blah")


    data$cost_t$value <- cost$par

    data$price_t$value <- price

    p_msy <- price*msy - cost$par * e_msy ^ data$beta

    data$p_msy <- p_msy

    data$e_msy <- e_msy

    data$p_expansion <- p_expansion

    inits <-
      map(
        1:chains,
        ~ list(
          base_effort = 1 *  exp(rnorm(1, 0, .1)),
          p_length_50_sel = 0.25 * exp(rnorm(1, 0, .1)))
      )

    # inits <-
    #   map(
    #     1:chains,
    #     ~ list(
    #       base_effort = 1 *  exp(rnorm(1, 0, .1)),
    #       p_length_50_sel = 0.25 * exp(rnorm(1, 0, .1)),
    #       p_expansion = p_expansion * exp(rnorm(1, 0, .1))
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
