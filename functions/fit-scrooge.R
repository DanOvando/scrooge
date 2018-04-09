fit_scrooge <- function(data, scrooge_file = "scrooge",
                        chains = 1, refresh = 25, cores = 1,
                        iter = 1000,
                        warmup = 2000,
                        adapt_delta = 0.8, economic_model = NA,
                        model_type = "scrooge",
                        max_treedepth = 10,
                        pmsy_expansion = 0.5){


  if (model_type == "scrooge"){

  if (is.na(economic_model) == F){
  data <- purrr::list_modify(data,economic_model = economic_model)
  }
  data$sigma_r_guess <- 0.4

  fmsy <- nlminb(data$m, calc_msy, data = data, time = 10, lower = 0, upper = 2)

  pmsy <- mean(data$price_t$value) * -fmsy$objective - mean(data$cost_t$value) * (fmsy$par / mean(data$q_t$value)) ^ data$beta

  p_expansion = ((fmsy$par / mean(data$q_t$value)) * pmsy_expansion) / pmsy

  # data$p_expansion <- p_expansion

  inits <-
    map(1:chains,  ~ list(base_effort = 1 *  exp(rnorm(
      1, 0, .1
    )),
    p_length_50_sel = 0.25 *exp(rnorm(
      1, 0, .1
    ),
    p_expansion = p_expansion * exp(rnorm(1,0,.1)))))


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

# init = list(list(effort_t = true_effort_devs$effort_devs)


  # clean up old DLLS per https://github.com/stan-dev/rstan/issues/448
  loaded_dlls = getLoadedDLLs()
  loaded_dlls = loaded_dlls[str_detect(names(loaded_dlls), '^file')]
  if (length(loaded_dlls) > 10) {
    for (dll in head(loaded_dlls, -10)) {
      # message("Unloading DLL ", dll[['name']], ": ", dll[['path']])
      dyn.unload(dll[['path']])
    }
  }
  # message("DLL Count = ", length(getLoadedDLLs()), ": [", str_c(names(loaded_dlls), collapse = ","), "]")
  }

  if (model_type == "lime"){

    fish <- vfo$prepped_fishery[[1]]$fish

    fleet <- vfo$prepped_fishery[[1]]$fleet

    scrooge_lh <- create_lh_list(vbk= fish$vbk,
                                 linf= fish$linf,
                                 t0= fish$t0,
                                 lwa= fish$weight_a,
                                 lwb= fish$weight_b,
                                 S50= fleet$length_50_sel,
                                 S95=fleet$length_95_sel,
                                 selex_input="length",
                                 selex_type=c("logistic"),
                                 M50= fish$length_50_mature,
                                 M95= fish$length_95_mature,
                                 maturity_input="length",
                                 M= fish$m,
                                 binwidth=1,
                                 CVlen= fish$cv_len,
                                 SigmaR= fish$sigma_r + .001,
                                 SigmaF= fleet$sigma_effort + .001,
                                 SigmaC=0.2,
                                 SigmaI=0.2,
                                 R0= fish$r0,
                                 qcoef=1e-5,
                                 start_ages=0,
                                 rho=0,
                                 nseasons=1)

    true <- generate_data(modpath=NULL,
                          itervec=1,
                          Fdynamics="Endogenous",
                          Rdynamics="BH",
                          lh=scrooge_lh,
                          Nyears=20,
                          Nyears_comp=20,
                          comp_sample=200,
                          init_depl=0.5,
                          seed=123)

    temp_LF_matrix <- vfo$prepped_fishery[[1]]$length_comps


    LF_matrix <- temp_LF_matrix %>%
      as.matrix()

    rownames(LF_matrix) <- LF_matrix[,"year"]


    LF_matrix <- LF_matrix[, -c(1)]

    scrooge_data_LF <-
      list("years" = 1:nrow(LF_matrix), "LF" = LF_matrix)


    start <- Sys.time()
    res <- run_LIME(modpath=NULL,
                    lh=scrooge_lh,
                    input_data=scrooge_data_LF,
                    est_sigma="log_sigma_R",
                    data_avail="LC",
                    newtonsteps=3)
    end <- Sys.time() - start
  }

  return(fit)

}