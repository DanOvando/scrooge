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
library(furrr)
library(tidyverse)
rstan::rstan_options(auto_write = TRUE)
extrafont::loadfonts()
functions <- list.files(here::here("functions"))

walk(functions, ~ here::here("functions", .x) %>% source()) # load local functions

in_clouds <-  F

run_name <- "d1.0"

local_data <- T

new_run <- F

scrooge_file <- "scrooge"

if (in_clouds == F) {
  run_dir <- here::here("results", run_name)

} else{
  if (new_run == T) {
    set.seed(sample(1:200, 1))
    run_name <- sample(fruit, 1)
    set.seed(42)
  }

  run_dir <- here::here("results", "scrooge-results", run_name)

}


if (dir.exists(run_dir) == F) {
  dir.create(run_dir, recursive = T)
}

run_description <- "making sure everything still works"


scrooge_theme <- theme_ipsum(base_size = 14, axis_title_size = 18)

theme_set(scrooge_theme)

# run options -------------------------------------------------------------

sim_fisheries <- FALSE

run_case_studies <- FALSE

run_clouds <- FALSE

fit_models <- FALSE

process_cloud_fits <- FALSE

n_cores <- 2

n_chains <- 1

if (in_clouds == T) {



  if (dir.exists("results/scrooge-results") == T) {

    system("umount results/scrooge-results")

    system("rm -r results/scrooge-results")

    system("mkdir results/scrooge-results")

  } else{

    system("mkdir results/scrooge-results")

  }

  system("gcsfuse scrooge-results results/scrooge-results")



  if (dir.exists("results/scrooge-data") == T) {

    system("umount data/scrooge-data")

    system("rm -r data/scrooge-data")

    system("mkdir data/scrooge-data")

  } else{

    system("mkdir data/scrooge-data")


  }

  # system("mkdir data/scrooge-data")

  system("gcsfuse scrooge-data data/scrooge-data")

  cloud_dir <- here::here("results", "scrooge-results", run_name)

  if (dir.exists(cloud_dir) == F) {
    dir.create(cloud_dir)

  }

}

write(run_description,
      file = glue::glue("{run_dir}/description.txt"))

# load data ---------------------------------------------------------------

if (in_clouds == F | local_data == T) {
  data_dir <- "data"
} else{
  data_dir <- "data/scrooge-data"
}

# load some sample cdfw catches

cdfw_data <-
  read_csv(file = file.path(data_dir, 'cfdw-catches.csv'))


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
  readr::read_csv(file.path(data_dir, "Exvessel Price Database.csv")) %>%
  filter(Year > 2010) %>%
  group_by(scientific_name) %>%
  summarise(mean_exvessel_price = mean(exvessel, na.rm = T)) %>%
  mutate(price_units = "usd_per_ton") %>%
  mutate(price_per_kg = mean_exvessel_price / 1000)

possible_prices <-
  na.omit(fish_prices$price_per_kg) %>% as.numeric()


# load in some sample RAM data

load(file.path(data_dir, "DBdata.RData")) # as of copy, RAM v4.40 with model fits

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
  gather(stock, f,-year)

effort <- effort.data %>%
  mutate(year = rownames(.)) %>%
  select(year, everything()) %>%
  gather(stock, effort,-year)


e_and_f <- f %>%
  left_join(effort, by = c('year', 'stock')) %>%
  filter(is.na(f) == F & is.na(effort) == F) %>%
  mutate(q = f / effort) %>%
  mutate(year = as.numeric(year))

possible_q <- e_and_f$q


u <- divupref.data %>%
  mutate(year = rownames(.)) %>%
  select(year, everything()) %>%
  gather(stock, u_v_umsy,-year) %>%
  mutate(year = as.numeric(year)) %>%
  arrange(stock, year) %>%
  group_by(stock) %>%
  mutate(lead_u = lead(u_v_umsy, 1)) %>%
  filter(is.na(u_v_umsy) == F & is.na(lead_u) == F) %>%
  mutate(delta_u = lead_u / u_v_umsy - 1) %>%
  group_by(stock) %>%
  summarise(max_delta_u = pmin(.4, max(delta_u))) %>%
  filter(max_delta_u > 0, is.finite(max_delta_u))

possible_delta_u <- u$max_delta_u

max_f_v_fmsy_increase <- median(possible_delta_u)

# create simulated fisheries ----------------------------------------------

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


if (sim_fisheries == T)
{
  n_fisheries <- 10

  candidate_species <- paste(FishLife::database$Z_ik$Genus,FishLife::database$Z_ik$Species) %>%
    unique() %>%
    sort()

    fisheries_sandbox <- data_frame(
    sci_name = sample(
      "Lepidopsetta polyxystra"
      ,
      n_fisheries,
      replace = T
    ),
    fleet_model = sample(c("open-access"),
                         n_fisheries,
                         replace = T),
    sigma_r = runif(n_fisheries, 0.01, .7),
    rec_ac = runif(n_fisheries, 0, 0.4),
    sigma_effort = runif(n_fisheries, 0, 0),
    price_cv = runif(n_fisheries, 0, 0.2),
    cost_cv = runif(n_fisheries, 0, 0.2),
    q_cv = runif(n_fisheries, 0, 0.2),
    price_slope = runif(n_fisheries, 0, 0.005),
    cost_slope = runif(n_fisheries,-0.005, 0.005),
    q_slope = runif(n_fisheries, 0, 0.005),
    price_ac = runif(n_fisheries, 0.5, 0.75),
    cost_ac = runif(n_fisheries, 0.5, 0.75),
    q_ac = runif(n_fisheries, 0.5, 0.75),
    steepness = runif(n_fisheries, 0.6, 0.9),
    obs_error = runif(n_fisheries, 0, 0),
    max_cp_ratio = runif(n_fisheries, 0.01, .95),
    price = sample(possible_prices, n_fisheries, replace = T),
    r0 = sample(possible_r0, n_fisheries, replace = T),
    q = sample(possible_q, n_fisheries, replace = T),
    max_perc_change_f = sample(possible_delta_u, n_fisheries, replace = T),
    profit_lags = sample(0, n_fisheries, replace = T),
    initial_f = sample(c(0.01, .1), n_fisheries, replace = T),
    beta = runif(n_fisheries, 2, 2),
    percnt_loo_selected = runif(n_fisheries, 0.1, 0.6)
  )

  fisheries_sandbox <- fisheries_sandbox %>%
    left_join(fleet_model_params, by = "fleet_model")
  # slice(1:4)

  if (file.exists("sim_progress.txt")) {
    file.remove("sim_progress.txt")
  }

  a <- Sys.time()
  spf <- safely(prepare_fishery)

  doParallel::registerDoParallel(cores = n_cores)

  foreach::getDoParWorkers()

  prepped_fishery <-
    foreach::foreach(i = 1:nrow(fisheries_sandbox)) %dopar% {
      write(
        glue::glue("{i/nrow(fisheries_sandbox)*100}% done with sims"),
        "sim_progress.txt",
        append = T
      )

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
        linf_buffer = 1.25
      )

    } # close dopar

  Sys.time() - a

  fisheries_sandbox$prepped_fishery <- prepped_fishery

  no_error <- map(fisheries_sandbox$prepped_fishery, 'error') %>%
    map_lgl(is.null)

  fisheries_sandbox <- fisheries_sandbox %>%
    filter(no_error == T)

  fisheries_sandbox <- fisheries_sandbox %>%
    mutate(prepped_fishery = map(prepped_fishery, "result")) %>%
    mutate(summary_plot = map(prepped_fishery, plot_simmed_fishery))

  fisheries_sandbox$fishery <- 1:nrow(fisheries_sandbox)


  if (in_clouds == FALSE) {
    saveRDS(fisheries_sandbox,
            file = here::here("processed_data", "fisheries_sandbox.RDS"),)
  } else {
    saveRDS(
      fisheries_sandbox,
      file = here::here("data", "scrooge-data", "fisheries_sandbox.RDS"),
    )

  }

} else{
  if (in_clouds == F | local_data == T) {
    fisheries_sandbox <-
      readRDS(file = here::here("processed_data", "fisheries_sandbox.RDS"))
  } else{
    fisheries_sandbox <-
      readRDS(file = here::here("data", "scrooge-data", "fisheries_sandbox.Rdata"))

  }
} # close sim fisheries stuff


# apply candidate assessment models ---------------------------------------

# test three cases  studies to demonstrate core model performance


if (run_case_studies == T) {
  # simple open access

  simple <- prepare_fishery(
    sci_name = "Lutjanus campechanus",
    fleet_model = "open-access",
    sigma_r = 0,
    rec_ac = 0,
    sigma_effort = 0,
    price_cv = 0,
    cost_cv = 0,
    q_cv = 0,
    price_slope = 0,
    cost_slope = 0,
    q_slope = 0,
    price_ac = 0,
    cost_ac = 0,
    q_ac = 0,
    steepness = 0.6,
    percnt_loo_selected = 0.1,
    obs_error = 0,
    initial_f = .01,
    r0 = 100,
    price = 2,
    q = 1e-3,
    profit_lags = 0,
    max_perc_change_f = 0.025,
    max_cp_ratio = 0.5,
    beta = 2,
    sim_years = 100,
    burn_years = 100,
    cv_len = 0.2,
    linf_buffer = 1.25,
    seed = 24
  )

  simple_plot <- plot_simmed_fishery(simple)

  simple_plot


  # realistic open access

  realistic <- prepare_fishery(
    sci_name = "Lutjanus campechanus",
    fleet_model = "open-access",
    sigma_r = 0.4,
    rec_ac = 0.5,
    sigma_effort = 0,
    price_cv = 0.4,
    cost_cv = 0.2,
    q_cv = 0.5,
    price_slope = 0,
    cost_slope =  0,
    q_slope = 0.0025,
    price_ac = 0.5,
    cost_ac = 0.5,
    q_ac = 0.75,
    steepness = 0.8,
    percnt_loo_selected = 0.3,
    obs_error = 0,
    initial_f = .025,
    r0 = 100,
    price = 4,
    q = 1e-3,
    profit_lags = 0,
    max_perc_change_f = 0.1,
    max_cp_ratio = .75,
    beta = 2,
    sim_years = 100,
    burn_years = 100,
    cv_len = 0.2,
    linf_buffer = 1.25,
    seed = 52
  )

  realistic_plot <- plot_simmed_fishery(realistic)

  realistic_plot

  realistic$simed_fishery %>%
    group_by(year) %>%
    summarise(rec_dev = unique(exp(rec_dev))) %>%
    ggplot(aes(year, rec_dev)) +
    geom_point()

  realistic$simed_fishery %>%
    group_by(year) %>%
    summarise(
      profits = sum(profits),
      effort = unique(effort),
      f = unique(f)
    ) %>%
    ungroup() %>%
    mutate(
      ppue = profits / effort,
      delta_effort = lead(effort) - effort,
      delta_f = lead(f) - f
    ) %>%
    ggplot(aes(year, ppue)) +
    geom_point()


  # decoupled

  decoupled <- prepare_fishery(
    sci_name = "Lutjanus campechanus",
    fleet_model = "supplied-catch",
    sigma_r = 0.5,
    rec_ac = 0.5,
    sigma_effort = 0,
    price_cv = 0.4,
    cost_cv = 0.2,
    q_cv = 0.5,
    price_slope = .0075,
    cost_slope =  -0.001,
    q_slope = 0.005,
    price_ac = 0.75,
    cost_ac = 0.75,
    q_ac = 0.75,
    steepness = 0.9,
    percnt_loo_selected = 0.3,
    obs_error = 0,
    initial_f = .025,
    r0 = 100,
    price = 4,
    q = 1e-3,
    profit_lags = 0,
    max_perc_change_f = 0.2,
    max_cp_ratio = 0.01,
    beta = 2,
    sim_years = 100,
    burn_years = 100,
    cv_len = 0.2,
    linf_buffer = 1.25,
    fleet_params = fleet_model_params$fleet_params[fleet_model_params$fleet_model == "supplied-catch"][[1]]
  )

  decoupled_plot <- plot_simmed_fishery(decoupled)

  case_studies <- data_frame(
    case_study = c(
                   "realistic",
                   "decoupled"),
    prepped_fishery  = list(
                            realistic,
                            decoupled)
  )

  experiments <- expand.grid(
    period = c("middle"),
    window = c(15),
    economic_model = c(0, 1),
    likelihood_model = c(0, 1, 2),
    prop_years_lcomp_data = c(0.2, 1),
    case_study = unique(case_studies$case_study),
    stringsAsFactors = F
  )

  case_studies <- experiments %>%
    left_join(case_studies, by = "case_study") %>%
    mutate(experiment = 1:nrow(.)) %>%
    mutate(prepped_fishery = pmap(
      list(
        prepped_fishery = prepped_fishery,
        window = window,
        period = period,
        experiment = experiment,
        prop_years_lcomp_data = prop_years_lcomp_data
      ),
      subsample_data
    )) %>%
    # filter(
    #   case_study == "realistic",
    #   period == "beginning",
    #   window == 15,
    #   prop_years_lcomp_data == 0.5
    # ) %>%
    filter(!(likelihood_model == 2 & economic_model == 3)) %>%
    filter(!(likelihood_model == 1 & economic_model == 2))

  scrooge_model <-
    rstan::stan_model(here::here("src", paste0(scrooge_file, ".stan")), model_name = scrooge_file)

  sfs <- safely(fit_scrooge)

  set.seed(42)

  # case_studies <- case_studies %>%
  #   filter(economic_model == 1, likelihood_model == 1) %>%
  #   slice(1:2)

  future::plan(future::multiprocess, workers = n_cores)

  case_study_fits <- future_pmap(
    list(
      data = map(case_studies$prepped_fishery, "scrooge_data"),
      fish = map(case_studies$prepped_fishery, "fish"),
      fleet = map(case_studies$prepped_fishery, "fleet"),
      experiment = case_studies$experiment,
      economic_model =  case_studies$economic_model,
      likelihood_model = case_studies$likelihood_model
    ),
    .f = sfs,
    chains = 1,
    cores = 1,
    refresh = 25,
    scrooge_model = scrooge_model,
    scrooge_file = "scrooge",
    iter = 2000,
    warmup = 1000,
    adapt_delta = 0.95,
    max_treedepth = 8,
    max_perc_change_f = 0.2,
    in_clouds = in_clouds,
    q_guess = mean(possible_q),
    r0 = 100,
    sd_sigma_r = 0.4,
    cv_obs = 2,
    .progress = TRUE
  )

  case_studies$scrooge_fit <- map(case_study_fits, "result")

  saveRDS(case_studies, file = glue::glue("{run_dir}/case_studies.RDS"))

  case_studies <- case_studies %>%
    mutate(performance = map2(prepped_fishery, scrooge_fit, safely(assess_fits)))

  case_studies$performance <- case_studies$performance %>% map("result")

  perf_summaries <- case_studies %>%
    mutate(others = map(performance, "others")) %>%
    select(-scrooge_fit,-prepped_fishery, -performance) %>%
    unnest() %>%
    group_by(
      year,
      variable,
      experiment,
      case_study,
      economic_model,
      likelihood_model,
      period,
      window,
      prop_years_lcomp_data
    ) %>%
    summarise(
      lower_90 = quantile(predicted, 0.05),
      upper_90 = quantile(predicted, 0.95),
      lower_50 = quantile(predicted, 0.25),
      upper_50 = quantile(predicted, 0.75),
      mean_predicted = mean(predicted),
      mean_predicted = median(mean_predicted),
      observed = mean(observed)
    )

  # perf_summaries %>%
  #   filter(variable == "f") %>%
  #   ggplot() +
  #   geom_ribbon(aes(year, ymin = lower_90, ymax = upper_90), fill = "lightgrey") +
  #   geom_ribbon(aes(year, ymin = lower_50, ymax = upper_50), fill = "darkgrey") +
  #   geom_line(aes(year, mean_predicted), color = "steelblue") +
  #   geom_point(
  #     aes(year, observed),
  #     fill = "tomato",
  #     size = 4,
  #     shape = 21
  #   ) +
  #   labs(y = "", x = "Year") +
  #   facet_wrap( ~ experiment, scales = "free") +
  #   theme_minimal()

  length_comps <- map(case_studies$performance, "length_comps")

#
#   length_comps[[1]] %>%
#     filter(source == "posterior_predictive") %>%
#     group_by(year, .chain, .iteration) %>%
#     mutate(predicted = predicted / sum(predicted)) %>%
#     group_by(year, .chain, length_bin) %>%
#     summarise(
#       lower_90 = quantile(predicted, 0.05),
#       upper_90 = quantile(predicted, 0.95),
#       mean = mean(predicted),
#       observed = unique(observed)
#     ) %>%
#     ggplot() +
#     geom_ribbon(aes(x = length_bin, ymin = lower_90, ymax = upper_90), fill = "lightgrey") +
#     geom_line(aes(length_bin, mean), color = "steelblue") +
#     geom_point(
#       aes(length_bin, observed),
#       size = .5,
#       alpha = 0.5,
#       color = "red"
#     ) +
#     facet_wrap( ~ year) +
#     theme_minimal()

  saveRDS(case_studies, file = glue::glue("{run_dir}/case_studies.RDS"))

  saveRDS(perf_summaries, file = glue::glue("{run_dir}/perf_summaries.RDS"))

} else {

  case_studies <- readRDS(glue::glue("{run_dir}/case_studies.RDS"))

  perf_summaries <- readRDS(glue::glue("{run_dir}/perf_summaries.RDS"))

}# close case studies runs

if (run_clouds == T){
if (fit_models == T) {
  experiments <- expand.grid(
    period = c("beginning", "middle", "end"),
    window = c(5, 10, 15),
    economic_model = c(0, 1, 2, 3),
    likelihood_model = c(0, 1, 2),
    prop_years_lcomp_data = c(0.1, .5, 1),
    fishery = unique(fisheries_sandbox$fishery),
    stringsAsFactors = F
  )

  fisheries_sandbox <- experiments %>%
    left_join(fisheries_sandbox, by = "fishery") %>%
    sample_n(2) %>%
    mutate(experiment = 1:nrow(.)) %>%
    mutate(prepped_fishery = pmap(
      list(
        prepped_fishery = prepped_fishery,
        window = window,
        period = period,
        experiment = experiment,
        prop_years_lcomp_data = prop_years_lcomp_data
      ),
      subsample_data
    ))


  scrooge_model <-
    rstan::stan_model(here::here("src", paste0(scrooge_file, ".stan")), model_name = scrooge_file)

  sfs <- safely(fit_scrooge)

  future::plan(future::multiprocess, workers = n_cores)

  fitted_sandbox <- future_pmap(
    list(
      data = map(fisheries_sandbox$prepped_fishery, "scrooge_data"),
      fish = map(fisheries_sandbox$prepped_fishery, "fish"),
      fleet = map(fisheries_sandbox$prepped_fishery, "fleet"),
      experiment = fisheries_sandbox$experiment,
      economic_model =  fisheries_sandbox$economic_model,
      likelihood_model = fisheries_sandbox$likelihood_model
    ),
    .f = sfs,
    chains = n_chains,
    cores = n_chains,
    refresh = 200,
    scrooge_model = scrooge_model,
    scrooge_file = "scrooge",
    iter = 2000,
    warmup = 1000,
    adapt_delta = 0.9,
    max_treedepth = 8,
    max_perc_change_f = 0.2,
    in_clouds = in_clouds,
    cloud_dir = cloud_dir,
    q_guess = mean(possible_q),
    r0 = 100,
    sd_sigma_r = 0.4,
    cv_obs = 2,
    .progress = TRUE
  )

  fisheries_sandbox$scrooge_fit <- fitted_sandbox

  save(file = glue::glue("{run_dir}/fitted_fisheries_sandbox.Rdata"),
       fisheries_sandbox)


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

scrooge_worked <-
  map(fisheries_sandbox$scrooge_fit, 'error') %>% map_lgl(is.null)

fisheries_sandbox <- fisheries_sandbox %>%
  filter(scrooge_worked) %>%  #,stan_worked, lime_worked, period != "middle")
  mutate(scrooge_fit = map(scrooge_fit, "result"))

if (in_clouds == T) {

  fisheries_sandbox <- fisheries_sandbox %>%
    mutate(performance = map2(scrooge_fit, prepped_fishery, summarise_performance, cloud_dir = cloud_dir))

} else {
  stop("must have in_clouds == T to run full diagnostics, never going to store them locally for now")
}

fisheries_sandbox <- fisheries_sandbox %>%
  mutate(mean_rmse = map_dbl(performance, ~sqrt(mean(.x$sq_er))))

fisheries_sandbox$performance[[1]] %>%
  group_by(year) %>%
  summarise(predicted = mean(predicted),
            observed = mean(observed)
            ) %>%
  ggplot() +
  geom_line(aes(year, predicted)) +
  geom_point(aes(year, observed), color = "red")

} # close cloud


# Process cloud fits ------------------------------------------------------
process_cloud_fits <-  FALSE
if (process_cloud_fits == T){


if (dir.exists(here::here("results","scrooge-results",run_name)) == F){

  system("gcsfuse scrooge-results results/scrooge-results")

}
# system("umount results/scrooge-results")
#
# system("rm -r results/scrooge-results")
#
# system("mkdir results/scrooge-results")
#
# system("gcsfuse scrooge-results results/scrooge-results")

processed_sandbox <- readRDS(file = glue::glue("results/scrooge-results/{run_name}/processed_fisheries_sandbox.RDS"))

reserve <- processed_sandbox$performance

processed_sandbox$performance <- map(processed_sandbox$performance, "result")

performance_worked <- map_dbl(map(processed_sandbox$performance, class), length) %>%
  map_lgl(~.x >1)

sandbox_performance <- processed_sandbox %>%
  select(-scrooge_fit, -prepped_fishery,-summary_plot,-fleet_params) %>%
  filter(performance_worked) %>%
  unnest()


sandbox_performance <- sandbox_performance%>%
  group_by(fishery,
           experiment,
           period,
           window,
           prop_years_lcomp_data,
           variable,
           economic_model,
           likelihood_model) %>%
  summarise(rmse = sqrt(mean(sq_er)),
            bias = median((predicted - observed) / observed)) %>%
  mutate(fit_model = glue::glue("em{economic_model}_lm{likelihood_model}"))


model_performance <- rstanarm::stan_glmer(rmse ~ fishery + (1|economic_model),
                                          data = sandbox_performance)



}


# the idea


# run diagnostics ---------------------------------------------------------


# make figures ------------------------------------------------------------


