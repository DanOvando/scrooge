fit_scrooge <- function(data,fish, fleet,scrooge_file = "scrooge",
                        chains = 1, refresh = 25, cores = 1,
                        iter = 1000,
                        warmup = 2000,
                        adapt_delta = 0.8, economic_model = 1,
                        model_type = "scrooge",
                        max_treedepth = 10,
                        pmsy_expansion = 0.5){


  data$sigma_r_guess <- 0.4
  fmsy <- nlminb(data$m, est_msy, data = data, time = 10, lower = 0, upper = 2)

  pmsy <- mean(data$price_t$value) * -fmsy$objective - mean(data$cost_t$value) * (fmsy$par / mean(data$q_t$value)) ^ data$beta

  p_expansion = ((fmsy$par / mean(data$q_t$value)) * pmsy_expansion) / pmsy

  inits <-
    map(1:chains,  ~ list(base_effort = 1 *  exp(rnorm(
      1, 0, .1
    )),
    p_length_50_sel = 0.25 *exp(rnorm(
      1, 0, .1
    )),
    p_expansion = p_expansion * exp(rnorm(1,0,.1))))


fit <-
    rstan::stan(
      file = here::here("src", paste0(scrooge_file,".stan")),
      data = data,
      chains = chains,
      refresh = refresh,
      cores = cores,
      iter = iter,
      warmup = warmup,
      control = list(
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth),
      init = inits
    )



  # clean up old DLLS per https://github.com/stan-dev/rstan/issues/448

    loaded_dlls = getLoadedDLLs()

    loaded_dlls = loaded_dlls[str_detect(names(loaded_dlls), '^file')]
    if (length(loaded_dlls) > 25) {
      for (dll in head(loaded_dlls, -15)) {
        # message("Unloading DLL ", dll[['name']], ": ", dll[['path']])
        dyn.unload(dll[['path']])
      }
    }

    out <- fit

  # out <-  list(scrooge_fit = fit,
  #      lime_fit = res)

  return(out)

}
