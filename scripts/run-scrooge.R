  # run-scrooge -------
  # Author: Dan Ovando
  # Project: scrooge
  # Summary:
  # Scrooge is a model for integrating economic priors into length-based fisheries
  # stock assessment. This script applies scrooge to a range of simulated fisheries and compares
  # performance among scenarios and across alternative assessment methods


  # setup -------------------------------------------------------------------
  set.seed(42)
  library(rstan)
  library(FishLife)
  library(spasm) # personal fishery simulator, will post public github link
  library(patchwork)
  library(hrbrthemes)
  library(foreach)
  library(doParallel)
  library(extrafont)
  library(LIME)
  library(LBSPR)
  library(stringr)
  library(glue)
  library(broom)
  library(tidyverse)
  rstan::rstan_options(auto_write = TRUE)
  extrafont::loadfonts()
  functions <- list.files(here::here("functions"))

  walk(functions, ~ here::here("functions", .x) %>% source()) # load local functions

  in_clouds <-  F

  run_name <- "v1.0"

  local_data <- T

  new_run <- F

  if (in_clouds == F){

    run_dir <- here::here("results", run_name)

  } else{

    if (new_run == T){
      set.seed(sample(1:200,1))
      run_name <- sample(fruit,1)
      set.seed(42)
    }

    run_dir <- here::here("results","scrooge-results", run_name)

  }


  if (dir.exists(run_dir) == F) {
    dir.create(run_dir, recursive = T)
  }

  run_description <- "making sure everything still works"


  scrooge_theme <- theme_ipsum(base_size = 14, axis_title_size = 18)

  theme_set(scrooge_theme)

  # run options -------------------------------------------------------------

  sim_fisheries <- T

  fit_models <- F

  run_tests <- F

  n_cores <- 1

  if (in_clouds == T){


    system("umount results/scrooge-results")

    system("rm -r results/scrooge-results")

    if (dir.exists("results/scroote-results") == F){

       system("mkdir results/scrooge-results")

    }

    system("gcsfuse scrooge-results results/scrooge-results")

    system("umount data/scrooge-data")

    system("rm -r data/scrooge-data")

    if (dir.exists("results/scroote-data") == F){

      system("mkdir data/scrooge-data")

    }

    # system("mkdir data/scrooge-data")

    system("gcsfuse scrooge-data data/scrooge-data")

    cloud_dir <- here::here("results","scrooge-results",run_name)

    if (dir.exists(cloud_dir) == F){

    dir.create(cloud_dir)

    }

  }

  write(run_description,
        file = glue::glue("{run_dir}/description.txt"))

  # load data ---------------------------------------------------------------

  if (in_clouds == F | local_data == T){

    data_dir <- "data"
  } else{

    data_dir <- "data/scrooge-data"
  }

  # load some sample cdfw catches

  cdfw_data <- read_csv(file = file.path(data_dir, 'cfdw-catches.csv'))


  has_timeseries <- cdfw_data %>%
    group_by(sci_name) %>%
    summarise(has_catch = min(year) <= 2000 & max(year) == 2015) %>%
    filter(has_catch == T)

  fill_catches <- function(min_year, max_year, catches) {
    full_frame <- data_frame(year = min_year:max_year) %>%
      left_join(catches, by = 'year') %>%
      mutate(catch = zoo::na.approx(catch_lbs, rule = 2)) %>%
      select(-catch_lbs)

  }

  cdfw_catches <- cdfw_data %>%
    filter(sci_name %in% unique(has_timeseries$sci_name)) %>%
    group_by(sci_name, year) %>%
    summarise(catch_lbs = sum(pounds_caught)) %>%
    group_by(sci_name) %>%
    nest(-sci_name, .key = 'catches') %>%
    filter(is.na(sci_name) == F) %>%
    mutate(catches = map(
      catches,
      fill_catches,
      min_year = 2000,
      max_year = 2015
    )) %>%
    unnest()

  # load in prices


  fish_prices <-
    readr::read_csv(file.path(data_dir,"Exvessel Price Database.csv")) %>%
    filter(Year > 2010) %>%
    group_by(scientific_name) %>%
    summarise(mean_exvessel_price = mean(exvessel, na.rm = T)) %>%
    mutate(price_units = "usd_per_ton") %>%
    mutate(price_per_kg = mean_exvessel_price / 1000)

  possible_prices <-na.omit(fish_prices$price_per_kg) %>% as.numeric()


  # load in some sample RAM data

  load(file.path(data_dir,"DBdata.RData")) # as of copy, RAM v4.40 with model fits

  # r0 in numbers

  r0 <- bioparams %>%
    filter(bioid == "R0-E00") %>%
    mutate(biovalue = as.numeric(biovalue)) %>%
    mutate(ptile_biovalue = percent_rank(biovalue)) %>%
    filter(ptile_biovalue < 0.75, ptile_biovalue > 0.25)

  possible_r0 <- r0$biovalue

  # estiamtes of catchability

  f <- er.data %>%
    mutate(year = rownames(.)) %>%
    select(year, everything()) %>%
    gather(stock, f, -year)

  effort <- effort.data %>%
    mutate(year = rownames(.)) %>%
    select(year, everything()) %>%
    gather(stock, effort, -year)


  e_and_f <- f %>%
    left_join(effort, by = c('year','stock')) %>%
    filter(is.na(f) == F & is.na(effort) == F) %>%
    mutate(q = f / effort) %>%
    mutate(year = as.numeric(year))

  possible_q <- e_and_f$q


  u <- divupref.data %>%
    mutate(year = rownames(.)) %>%
    select(year, everything()) %>%
    gather(stock, u_v_umsy, -year) %>%
    mutate(year = as.numeric(year)) %>%
    arrange(stock, year) %>%
    group_by(stock) %>%
    mutate(lead_u = lead(u_v_umsy,1)) %>%
    filter(is.na(u_v_umsy) == F & is.na(lead_u) == F) %>%
    mutate(delta_u = lead_u / u_v_umsy - 1) %>%
    group_by(stock) %>%
    summarise(max_delta_u = pmin(.4,max(delta_u))) %>%
    filter(max_delta_u > 0, is.finite(max_delta_u))

  possible_delta_u <- u$max_delta_u

  max_f_v_fmsy_increase <- median(possible_delta_u)

  # create simulated fisheries ----------------------------------------------


  if (sim_fisheries == T)
  {

    n_fisheries <- 10

    # "constant-effort",
    # "supplied-catch",

    fisheries_sandbox <- data_frame(
      sci_name = sample(
        c(
          "Atractoscion nobilis",
          "Scomber japonicus",
          "Lutjanus campechanus"
        ),
        n_fisheries,
        replace = T
      ),
      fleet_model = sample(
        c(
          "open-access"),
        n_fisheries,
        replace = T
      ),
      sigma_r = runif(n_fisheries, 0.01,.7),
      rec_ac = runif(n_fisheries,0,0.4),
      sigma_effort = runif(n_fisheries, 0,0),
      price_cv = runif(n_fisheries, 0,0.2),
      cost_cv = runif(n_fisheries, 0,0.2),
      q_cv = runif(n_fisheries, 0,0.2),
      price_slope = runif(n_fisheries, 0,0.005),
      cost_slope = runif(n_fisheries, -0.005,0.005),
      q_slope = runif(n_fisheries, 0,0.005),
      price_ac = runif(n_fisheries, 0.5,0.75),
      cost_ac = runif(n_fisheries, 0.5,0.75),
      q_ac = runif(n_fisheries, 0.5,0.75),
      steepness = runif(n_fisheries, 0.6,0.9),
      obs_error = runif(n_fisheries, 0,0),
      max_cp_ratio = runif(n_fisheries, 0.01,.95),
      price = sample(possible_prices, n_fisheries, replace = T),
      r0 = sample(possible_r0, n_fisheries, replace = T),
      q = sample(possible_q, n_fisheries, replace = T),
      max_perc_change_f = sample(possible_delta_u, n_fisheries, replace = T),
      profit_lags = sample(0, n_fisheries, replace = T),
      initial_f = sample(c(0.01,.1), n_fisheries, replace = T),
      beta = runif(n_fisheries, 2,2),
      percnt_loo_selected = runif(n_fisheries, 0.1, 0.6)
    )
    fleet_model_params <- data_frame(
      fleet_model = c(
        "constant-catch",
        "constant-effort",
        "supplied-catch",
        "open-access"
      ),
      fleet_params = list(
        list(target_catch = 10000),
        list(initial_effort = NA),
        list(catches = cdfw_catches$catch[cdfw_catches$sci_name == "semicossyphus pulcher"]),
        list(theta = NA, initial_effort = NA)
      )
    )

    fisheries_sandbox <- fisheries_sandbox %>%
      left_join(fleet_model_params, by = "fleet_model")
      # slice(1:4)

    if (file.exists("sim_progress.txt")){file.remove("sim_progress.txt")}

    a <- Sys.time()
    spf <- safely(prepare_fishery)

    doParallel::registerDoParallel(cores = n_cores)

    foreach::getDoParWorkers()

    prepped_fishery <- foreach::foreach(i = 1:nrow(fisheries_sandbox)) %dopar% {

      write(glue::glue("{i/nrow(fisheries_sandbox)*100}% done with sims"), "sim_progress.txt", append = T)

      out <- spf(
          sci_name = fisheries_sandbox$sci_name[i],
          fleet_model = fisheries_sandbox$fleet_model[i],
          fleet_params = fisheries_sandbox$fleet_params[[i]],
          sigma_r = fisheries_sandbox$sigma_r[[i]],
          sigma_effort = fisheries_sandbox$sigma_effort[i],
          price_cv = fisheries_sandbox$price_cv[i],
          cost_cv = fisheries_sandbox$cost_cv[i],
          q_cv = fisheries_sandbox$q_cv[i],
          price_slope = fisheries_sandbox$price_slope[i],
          cost_slope = fisheries_sandbox$cost_slope[i],
          q_slope = fisheries_sandbox$q_slope[i],
          price_ac = fisheries_sandbox$price_ac[i],
          cost_ac = fisheries_sandbox$cost_ac[i],
          q_ac = fisheries_sandbox$q_ac[i],
          steepness = fisheries_sandbox$steepness[i],
          percnt_loo_selected = fisheries_sandbox$percnt_loo_selected[i],
          obs_error = fisheries_sandbox$obs_error[i],
          initial_f = fisheries_sandbox$initial_f[i],
          r0 = fisheries_sandbox$r0[i],
          price = fisheries_sandbox$price[i],
          q = fisheries_sandbox$q[i],
          profit_lags = fisheries_sandbox$profit_lags[i],
          max_perc_change_f = fisheries_sandbox$max_perc_change_f[i],
          max_cp_ratio = fisheries_sandbox$max_cp_ratio[i],
          beta = fisheries_sandbox$beta[i],
          sim_years = 100,
          burn_years = 100,
          cv_len = 0.2,
          linf_buffer = 1.1
      )

    } # close dopar

    Sys.time() - a

    fisheries_sandbox$prepped_fishery <- prepped_fishery

    no_error <- map(fisheries_sandbox$prepped_fishery,'error') %>%
      map_lgl(is.null)

    fisheries_sandbox <- fisheries_sandbox %>%
      filter(no_error == T)

    fisheries_sandbox <- fisheries_sandbox %>%
      mutate(prepped_fishery = map(prepped_fishery,"result")) %>%
      mutate(summary_plot = map(prepped_fishery, plot_simmed_fishery))


    save(file = here::here("processed_data", "fisheries_sandbox.Rdata"),
         fisheries_sandbox)

  } else{

    if (in_clouds == F | local_data == T){
    load(file = here::here("processed_data", "fisheries_sandbox.Rdata"))
    } else{

      load(file = here::here("data","scrooge-data","fisheries_sandbox.Rdata"))

    }
  } # close sim fisheries stuff

# apply candidate assessment models ---------------------------------------

# test three simple cases to verify core model performance

if (run_tests == T) {
  # perfecto open access

  pfo <- fisheries_sandbox %>%
    filter(fleet_model == "open-access",
           sigma_r == min(sigma_r),
           sigma_effort == min(sigma_effort),
           price_cv == min(price_cv),
           cost_cv == min(cost_cv)) %>%
    slice(1) %>%
    mutate(prepped_fishery = map(prepped_fishery, subsample_data, window = 10, period = "end")) %>%
    mutate(scrooge_fit = pmap(
      list(
      data = map(prepped_fishery, "scrooge_data"),
      fish = map(prepped_fishery, "fish"),
      fleet = map(prepped_fishery, "fleet")),
      fit_scrooge,
        iter = 4000,
        warmup = 2000,
        adapt_delta = 0.8,
        economic_model = 1,
        scrooge_file = "scrooge",
        in_clouds = F,
      experiment = "pfo",
      max_f_v_fmsy_increase = 0.1
      )
    )


  pfo$summary_plot[[1]]


  pfo <- pfo %>%
    mutate(processed_scrooge = map2(
      scrooge_fit,
      map(prepped_fishery, "sampled_years"),
      process_scrooge
    )) %>%
    mutate(observed = map(prepped_fishery, "simed_fishery")) %>%
    mutate(scrooge_performance = pmap(
      list(observed = observed,
           predicted = processed_scrooge),
      judge_performance
    )) %>%
    mutate(
      scrooge_rec_performance = pmap(
        list(observed = observed,
             predicted = processed_scrooge),
        judge_performance,
        observed_variable = rec_dev,
        predicted_variable = "rec_dev_t"
      )
    )  %>%
    mutate(lcomps = map2(prepped_fishery, processed_scrooge, ~process_lcomps(.x, .y$n_tl))) %>%
    mutate(
      rmse = map_dbl(scrooge_performance, ~ .x$comparison_summary$rmse),
      bias =  map_dbl(scrooge_performance, ~ .x$comparison_summary$bias)
    ) %>%
    arrange(rmse)

  pfo$scrooge_performance[[1]]$comparison_plot

  pfo$scrooge_rec_performance[[1]]$comparison_plot


  # variable open access

  vfo <- fisheries_sandbox %>%
    dplyr::filter(
      fleet_model == "open-access",
      sigma_r == max(sigma_r),
      sigma_effort == min(sigma_effort),
      price_cv == max(price_cv),
      price_ac == max(price_ac),
      cost_ac == max(cost_ac),
      cost_cv == max(cost_cv),
      b_v_bmsy_oa == 0.5
    ) %>%
    slice(7) %>%
    mutate(prepped_fishery = map(prepped_fishery, subsample_data, window = 10, period = "end")) %>%
    mutate(scrooge_fit = pmap(
      list(
        data = map(prepped_fishery, "scrooge_data"),
        fish = map(prepped_fishery, "fish"),
        fleet = map(prepped_fishery, "fleet")),
      fit_scrooge,
      iter = 4000,
      warmup = 2000,
      adapt_delta = 0.9,
      economic_model = 1,
      scrooge_file = "scrooge",
      in_clouds = F,
      experiment = "vfo",
      max_f_v_fmsy_increase = 0.5

    )
    )


  vfo <- vfo %>%
    mutate(lime_fit = pmap(list(
      data = map(prepped_fishery, "scrooge_data"),
      fish = map(prepped_fishery, "fish"),
      fleet = map(prepped_fishery, "fleet")
    ), fit_lime))

  true_f <- vfo$prepped_fishery[[1]]$simed_fishery %>%
    group_by(year) %>%
    summarise(f = unique(f)) %>%
    filter(year >=64)

 a <- vfo$lime_fit[[1]]$Sdreport %>%
    summary() %>%
   as.data.frame() %>%
   mutate(variable = rownames(.)) %>%
   filter(variable == "lF_y") %>%
   mutate(mean = exp(Estimate),
          upper = exp(Estimate + 1.96 * `Std. Error`),
          lower = exp(Estimate -  1.96 * `Std. Error`)) %>%
   mutate(year = 1:nrow(.)) %>%
   ggplot() +
   geom_pointrange(aes(year, mean, ymin = lower, ymax = upper)) +
   geom_point(data = true_f, aes(1:11, f), color = "red")

  predicted_f <- data_frame(f = Report$F_y)

  true_f %>%
    ggplot() +
    geom_point(aes(year,f,color = "True")) +
    geom_line(aes(year, predicted_f,color = "Predicted"))

  vfo$summary_plot[[1]]

  vfo <- vfo %>%
    mutate(processed_scrooge = map2(
      scrooge_fit,
      map(prepped_fishery, "sampled_years"),
      process_scrooge
    )) %>%    mutate(observed = map(prepped_fishery, "simed_fishery")) %>%
    mutate(
      scrooge_rec_performance = pmap(
        list(observed = observed,
             predicted = processed_scrooge),
        judge_performance,
        observed_variable = rec_dev,
        predicted_variable = "rec_dev_t"
      )
    ) %>%
    mutate(scrooge_performance = pmap(
      list(observed = observed,
           predicted = processed_scrooge),
      judge_performance
    )) %>%
    mutate(lcomps = map2(prepped_fishery, processed_scrooge, ~process_lcomps(.x$length_comps, .y$n_tl))) %>%
    mutate(
      rmse = map_dbl(scrooge_performance, ~ .x$comparison_summary$rmse),
      bias =  map_dbl(scrooge_performance, ~ .x$comparison_summary$bias)
    ) %>%
    arrange(rmse)

    # vfo$scrooge_rec_performance[[1]]$comparison_plot

     vfo$scrooge_performance[[1]]$comparison_plot +
      labs(title = "econ")

    save(file = "scrooge_performance.Rdata", pfo, vfo)

  # constant and medium f


  cfo <- fisheries_sandbox %>%
    filter(
      fleet_model == "constant-effort",
      sigma_r == max(sigma_r),
      sigma_effort == max(sigma_effort),
      price_cv == 0,
      cost_cv == 0
    ) %>%
    slice(1) %>%
    mutate(summary_plot = map(prepped_fishery, plot_simmed_fishery)) %>%
    mutate(scrooge_fit = map(
      prepped_fishery,
      ~ fit_scrooge(
        data = .x$scrooge_data,
        iter = 10000,
        warmup = 4000
      )
    ))

  cfo$summary_plot[[1]]

  cfo$prepped_fishery[[1]]$length_comps %>%
    gather(age, numbers,-year) %>%
    mutate(age = as.numeric(age)) %>%
    group_by(year) %>%
    summarise(ncaught = sum(numbers)) %>%
    ggplot(aes(year, ncaught)) +
    geom_point()

  cfo <- cfo %>%
    mutate(processed_scrooge = map(scrooge_fit, process_scrooge)) %>%
    mutate(observed = map(prepped_fishery, "simed_fishery")) %>%
    mutate(
      scrooge_rec_performance = pmap(
        list(observed = observed,
             predicted = processed_scrooge),
        judge_performance,
        observed_variable = rec_dev,
        predicted_variable = "rec_dev_t"
      )
    ) %>%
    mutate(scrooge_performance = pmap(
      list(observed = observed,
           predicted = processed_scrooge),
      judge_performance
    )) %>%
    mutate(
      rmse = map_dbl(scrooge_performance, ~ .x$comparison_summary$rmse),
      bias =  map_dbl(scrooge_performance, ~ .x$comparison_summary$bias)
    ) %>%
    arrange(rmse)
  #
  cfo$scrooge_performance[[1]]$comparison_plot

  cfo$scrooge_rec_performance[[1]]$comparison_plot


} # close test runs

if (fit_models == T) {

res <- fisheries_sandbox

  experiments <- expand.grid(period = c("middle","end"), window = c(10),
                             economic_model = c(0,3),
                             fishery = 1:nrow(fisheries_sandbox), stringsAsFactors = F)

  fisheries_sandbox <- fisheries_sandbox %>%
    mutate(fishery = 1:nrow(.)) %>%
    left_join(experiments, by = "fishery") %>%
    mutate(prepped_fishery = pmap(list(
      prepped_fishery = prepped_fishery,
      window = window,
      period = period), subsample_data))

  fisheries_sandbox$experiment <- 1:nrow(fisheries_sandbox)

  sfs <- safely(fit_scrooge)

  doParallel::registerDoParallel(cores = n_cores)

  foreach::getDoParWorkers()

  set.seed(42)
  fits <- foreach::foreach(i = 1:nrow(fisheries_sandbox)) %dopar% {
    out <- sfs(
      data = fisheries_sandbox$prepped_fishery[[i]]$scrooge_data,
      fish = fisheries_sandbox$prepped_fishery[[i]]$fish,
      fleet = fisheries_sandbox$prepped_fishery[[i]]$fleet,
      experiment = fisheries_sandbox$experiment[i],
      economic_model = fisheries_sandbox$economic_model[i],
      scrooge_file = "scrooge",
      iter = 2000,
      warmup = 1000,
      adapt_delta = 0.8,
      max_treedepth = 12,
      max_perc_change_f = 0.2,
      in_clouds = in_clouds,
      cloud_dir = cloud_dir,
      chains = 2,
      cv_effort = 0.5,
      q_guess = mean(possible_q),
      r0 = 100,
      sd_sigma_r = 0.1,
      cores = 2,
      sigma_effort_guess = 0.1,
      burn_f = 0.12
    )

  } # close fitting loop


  fisheries_sandbox$scrooge_fit <- fits

  save(
    file = glue::glue("{run_dir}/fitted_fisheries_sandbox.Rdata"),
    fisheries_sandbox
  )


} else{

  load(file = glue::glue("{run_dir}/fitted_fisheries_sandbox.Rdata"))

} # close model fitting


  # fit lbspr

  # fisheries_sandbox <- fisheries_sandbox %>%
  #   mutate(lbspr_fit = pmap(list(
  #     data = map(prepped_fishery, "scrooge_data"),
  #     fish = map(prepped_fishery, "fish"),
  #     fleet = map(prepped_fishery, "fleet")
  #   ), fit_lbspr))

  # fit lime

  # fisheries_sandbox <- fisheries_sandbox %>%
  #   mutate(lime_fit = pmap(list(
  #     data = map(prepped_fishery, "scrooge_data"),
  #     fish = map(prepped_fishery, "fish"),
  #     fleet = map(prepped_fishery, "fleet")
  #   ), fit_lime))

  # fisheries_sandbox <- fisheries_sandbox %>%
  #   mutate(processed_lime = map2(
  #   lime_fit,
  #   map(prepped_fishery, "sampled_years"),
  #   safely(process_lime)
  # ))

# process_fits ------------------------------------------------------------


scrooge_worked <- map(fisheries_sandbox$scrooge_fit,'error') %>% map_lgl(is.null)


stan_worked <-  map_lgl(map(fisheries_sandbox$scrooge_fit,  "result"),~!(nrow(.x) %>% is.null()))


fisheries_sandbox <- fisheries_sandbox %>%
  filter(scrooge_worked) %>%  #,stan_worked, lime_worked, period != "middle")
   mutate(scrooge_fit = map(scrooge_fit, "result"))

scrooge_converged <- fisheries_sandbox$scrooge_fit %>%
  map_lgl(~class(rstan::get_logposterior(.x)) == "list")

fisheries_sandbox <- fisheries_sandbox %>%
  filter(scrooge_converged)

if(in_clouds == T){

  loadfoo <- function(experiment, cloud_dir){
    readRDS(glue::glue("{cloud_dir}/{experiment}"))
  }

  fisheries_sandbox <- fisheries_sandbox %>%
    mutate(scrooge_fit = map(scrooge_fit, safely(loadfoo), cloud_dir = cloud_dir))

}


processed_sandbox <-  fisheries_sandbox%>%
mutate(processed_scrooge = map2(
    scrooge_fit,
    map(prepped_fishery, "sampled_years"),
    process_scrooge
  )) %>%
  mutate(observed = map(prepped_fishery, "simed_fishery")) %>%
  mutate(scrooge_performance = pmap(
    list(observed = observed,
         predicted = processed_scrooge),
    judge_performance
  )) %>%
  mutate(
    scrooge_rec_performance = pmap(
      list(observed = observed,
           predicted = processed_scrooge),
      judge_performance,
      observed_variable = rec_dev,
      predicted_variable = "rec_dev_t"
    )
  )  %>%
  # mutate(lime_performance = pmap(
  #   list(observed = observed,
  #        predicted = processed_lime),
  #   judge_lime
  # )) %>%
  mutate(lcomps = map2(prepped_fishery, processed_scrooge, ~process_lcomps(.x, .y$n_tl))) %>%
  mutate(
    rmse = map_dbl(scrooge_performance, ~ .x$comparison_summary$rmse),
    bias =  map_dbl(scrooge_performance, ~ .x$comparison_summary$bias)
  ) %>%
  arrange(rmse)


# run diagnostics ---------------------------------------------------------


# make figures ------------------------------------------------------------


processed_sandbox %>%
  ggplot(aes(factor(economic_model), rmse)) +
  geom_boxplot(alpha = 0.5) +
  ggbeeswarm::geom_beeswarm()

processed_sandbox %>%
  ggplot(aes(factor(economic_model), bias)) +
  geom_boxplot(alpha = 0.5) +
  ggbeeswarm::geom_beeswarm()





