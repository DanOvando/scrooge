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
library(ggridges)
library(scales)
library(recipes)
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

cloud_dir <- here::here("results", "scrooge-results", run_name)


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

run_lime <- FALSE

n_fisheries <- 10 # number of fisheries to simulate

n_cores <- 2

n_chains <- 1

max_realistic_f <- 3

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

  candidate_species <- paste(FishLife::database$Z_ik$Genus,FishLife::database$Z_ik$Species) %>%
    unique() %>%
    sort()

    fisheries_sandbox <- data_frame(
    sci_name = sample(
      candidate_species
      ,
      n_fisheries,
      replace = T
    ),
    fleet_model = sample(c("open-access",
                           "constant-effort",
                           "random-walk"),
                         n_fisheries,
                         replace = T),
    sigma_r = runif(n_fisheries, 0.01, .7),
    rec_ac = runif(n_fisheries, 0, 0.4),
    sigma_effort = runif(n_fisheries, 0, .1),
    price_cv = runif(n_fisheries, 0, 0.2),
    cost_cv = runif(n_fisheries, 0, 0.2),
    q_cv = runif(n_fisheries, 0, 0.1),
    price_slope = runif(n_fisheries, 0, 0.005),
    cost_slope = runif(n_fisheries,-0.005, 0.005),
    q_slope = runif(n_fisheries, 0, 0.005),
    price_ac = runif(n_fisheries, 0.25, 0.75),
    cost_ac = runif(n_fisheries, 0.25, 0.75),
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

  fisheries_sandbox$sigma_effort[fisheries_sandbox$fleet_model == "random-walk"] <- runif(sum(fisheries_sandbox$fleet_model == "random-walk"), 0.1, 0.4)


  fisheries_sandbox <- fisheries_sandbox %>%
    left_join(fleet_model_params, by = "fleet_model")


  spf <- safely(prepare_fishery)

  future::plan(future::multiprocess, workers = n_cores)

  fisheries_sandbox <- fisheries_sandbox %>%
    mutate(prepped_fishery =
             furrr::future_pmap(
               list(
                 sci_name = sci_name,
                 fleet_model = fleet_model,
                 fleet_params = fleet_params,
                 sigma_r = sigma_r,
                 sigma_effort = sigma_effort,
                 price_cv = price_cv,
                 cost_cv = cost_cv,
                 q_cv = q_cv,
                 price_slope = price_slope,
                 cost_slope = cost_slope,
                 q_slope = q_slope,
                 price_ac = price_ac,
                 cost_ac = cost_ac,
                 q_ac = q_ac,
                 steepness = steepness,
                 percnt_loo_selected = percnt_loo_selected,
                 obs_error = obs_error,
                 initial_f = initial_f,
                 r0 = r0,
                 price = price,
                 q = q,
                 profit_lags = profit_lags,
                 max_perc_change_f = max_perc_change_f,
                 max_cp_ratio = max_cp_ratio,
                 beta = beta),
               spf,
               sim_years = 100,
               burn_years = 100,
               cv_len = 0.2,
               linf_buffer = 1.25,
               .progress = T
             ))

  no_error <- map(fisheries_sandbox$prepped_fishery, 'error') %>%
    map_lgl(is.null)

  fisheries_sandbox <- fisheries_sandbox %>%
    filter(no_error == T)

  fisheries_sandbox <- fisheries_sandbox %>%
    mutate(prepped_fishery = map(prepped_fishery, "result")) %>%
    mutate(summary_plot = map(prepped_fishery, plot_simmed_fishery))

  fisheries_sandbox$fishery <- 1:nrow(fisheries_sandbox)

  bad_fisheries <- fisheries_sandbox %>%
    select(fishery,prepped_fishery) %>%
    mutate(max_f = map_dbl(prepped_fishery, ~.x$simed_fishery$f %>% max())) %>%
    filter(max_f > max_realistic_f)

  fisheries_sandbox <- fisheries_sandbox %>%
    filter(!(fishery %in% bad_fisheries$fishery)) %>%
    ungroup() %>%
    mutate(fishery = 1:nrow(.))

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
    sigma_effort = 0.1,
    price_cv = 0.4,
    cost_cv = 0.2,
    q_cv = 0.05,
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

  # realistic_plot

  # realistic$simed_fishery %>%
  #   group_by(year) %>%
  #   summarise(
  #     profits = sum(profits),
  #     effort = unique(effort),
  #     f = unique(f)
  #   ) %>%
  #   ungroup() %>%
  #   mutate(
  #     ppue = profits / effort,
  #     delta_effort = lead(effort) - effort,
  #     delta_f = lead(f) - f
  #   ) %>%
  #   ggplot(aes(year, ppue)) +
  #   geom_point()


  # decoupled

  decoupled <- prepare_fishery(
    sci_name = "Lutjanus campechanus",
    fleet_model = "random-walk",
    sigma_r = 0.5,
    rec_ac = 0.5,
    sigma_effort = 0.1,
    price_cv = 0.4,
    cost_cv = 0.2,
    q_cv = 0.05,
    price_slope = .0075,
    cost_slope =  -0.001,
    q_slope = 0.0025,
    price_ac = 0.75,
    cost_ac = 0.75,
    q_ac = 0.75,
    steepness = 0.9,
    percnt_loo_selected = 0.3,
    obs_error = 0,
    initial_f = 0.2,
    r0 = 100,
    price = 4,
    q = 1e-3,
    profit_lags = 0,
    max_perc_change_f = 0.1,
    max_cp_ratio = 0.01,
    beta = 2,
    sim_years = 100,
    burn_years = 100,
    cv_len = 0.2,
    linf_buffer = 1.25,
    seed = 32
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
    economic_model = c(0, 1, 2),
    likelihood_model = c(0, 1),
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
    filter(!(likelihood_model == 2 & economic_model == 3)) %>%
    filter(!(likelihood_model == 1 & economic_model == 2))

  scrooge_model <-
    rstan::stan_model(here::here("src", paste0(scrooge_file, ".stan")), model_name = scrooge_file)

  sfs <- safely(fit_scrooge)

  set.seed(42)

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
      median_predicted = median(mean_predicted),
      observed = mean(observed)
    )

  length_comps <- map(case_studies$performance, "length_comps")

  saveRDS(case_studies, file = glue::glue("{run_dir}/case_studies.RDS"))

  saveRDS(perf_summaries, file = glue::glue("{run_dir}/perf_summaries.RDS"))

} #else {
#
#   case_studies <- readRDS(glue::glue("{run_dir}/case_studies.RDS"))
#
#   perf_summaries <- readRDS(glue::glue("{run_dir}/perf_summaries.RDS"))
#
# }# close case studies runs

if (run_clouds == T){
if (fit_models == T) {
  experiments <- expand.grid(
    period = c("middle"),
    window = c(15),
    economic_model = c(0, 1, 2, 3),
    likelihood_model = c(0, 1, 2),
    prop_years_lcomp_data = c(0.25,1),
    fishery =   unique(fisheries_sandbox$fishery),
    stringsAsFactors = F
  ) %>%
    filter(!(economic_model == 2 & likelihood_model == 1),
           !(economic_model == 3 & likelihood_model == 2))

  fisheries_sandbox <- experiments %>%
    left_join(fisheries_sandbox, by = "fishery") %>%
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

  future::plan(future::multiprocess, workers = 5)

  processed_sandbox <- fisheries_sandbox %>%
    mutate(performance = future_map2(scrooge_fit, prepped_fishery, safely(summarise_performance), cloud_dir = cloud_dir, .progress = T))

  saveRDS(processed_sandbox, file = glue::glue("results/scrooge-results/{run_name}/processed_fisheries_sandbox.RDS"))

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
if (process_cloud_fits == T){


if (dir.exists(here::here("results","scrooge-results",run_name)) == F){

  system("gcsfuse scrooge-results results/scrooge-results")

}

processed_sandbox <- readRDS(file = glue::glue("results/scrooge-results/{run_name}/processed_fisheries_sandbox.RDS"))

# reserve <- processed_sandbox$performance

# processed_sandbox$performance <- map(processed_sandbox$performance, "result")

performance_worked <- map_lgl(map(processed_sandbox$performance, "error"), is.null)

processed_sandbox$performance <- map(processed_sandbox$performance,"result")

processed_sandbox <- processed_sandbox %>%
  filter(performance_worked)

# performance_worked2 <- map_dbl(map(processed_sandbox$performance, class), length) %>%
#   map_lgl(~.x >1)

sandbox_performance <- processed_sandbox %>%
  select(-scrooge_fit, -prepped_fishery,-summary_plot,-fleet_params) %>%
  unnest()

sandbox_performance <- sandbox_performance%>%
  group_by(fishery,
           experiment,
           fleet_model,
           period,
           window,
           prop_years_lcomp_data,
           variable,
           economic_model,
           likelihood_model) %>%
  summarise(rmse = sqrt(mean(sq_er)),
            mean_bias = mean(predicted - observed),
            median_bias = median(predicted-observed),
            percent_bias = median((predicted - observed) / observed),
            mare = median(abs(predicted - observed) / observed),
            rmedse = sqrt(median(sq_er)),
            w_rmse = sqrt(weighted.mean(sq_er, w = year))) %>%
  mutate(fit_model = glue::glue("em{economic_model}_lm{likelihood_model}"),
         log_rmse = log(rmse)) %>%
  left_join(processed_sandbox %>%
              select(experiment, sigma_r:percnt_loo_selected), by = "experiment")

sandbox_performance <- sandbox_performance %>%
  ungroup()


prep_performance <- sandbox_performance %>%
  select(w_rmse, fit_model,sigma_r:percnt_loo_selected) %>% {
  recipes::recipe(w_rmse ~ ., data = .)
  } %>%
  step_log(all_outcomes()) %>%
  step_nzv(all_predictors()) %>%
  step_center(all_numeric(),-all_outcomes()) %>%
  step_scale(all_numeric(),-all_outcomes()) %>%
  prep(data = sandbox_performance, retain = T) %>%
  juice()

model_formula <- paste0("w_rmse ~ ",paste(colnames(prep_performance %>% select(-w_rmse, -fit_model)), collapse = "+"),
"+ (1|fit_model)")

performance_model <- rstanarm::stan_glmer(model_formula,
                                          data = prep_performance, chains = 2, cores = 2)

tidy_performance <- tidybayes::spread_draws(performance_model, b[intercept,model])

rmse_effect_plot <- tidy_performance %>%
  group_by(model) %>%
  mutate(mrmse = mean(b)) %>%
  ungroup() %>%
  mutate(model = fct_reorder(model,b,.fun = mean, .desc = TRUE)) %>%
  ggplot(aes(b, fill = model)) +
  geom_vline(aes(xintercept = 0)) +
  geom_density_ridges(aes(x = b, y = model),alpha = 0.75, show.legend = FALSE)  +
  scale_fill_viridis_d(option = "E") +
  labs(x = "Effect on RMSE") +
  theme(axis.title.x = element_blank())

saveRDS(object = performance_model,file =  paste0(run_dir,"/performance_model.RDS"))

saveRDS(object = sandbox_performance,file =  paste0(run_dir,"/sandbox_performance.RDS"))


mod_select <- sandbox_performance %>%
  group_by(fishery,
           period,
           window,
           prop_years_lcomp_data) %>%
  filter(variable == "f") %>%
  filter(rmse == min(rmse)) %>%
  ungroup()


model_selected_plot <- mod_select %>%
  group_by(fit_model) %>%
  count() %>%
  ungroup() %>%
  mutate(p_selected = n / sum(n)) %>%
  mutate(fit_model = fct_reorder(fit_model, p_selected)) %>%
  ggplot(aes(fit_model, p_selected)) +
  geom_col(color = "black", fill = "lightgrey") +
  scale_y_continuous(labels = percent, name = "Percentage Best Model") +
  coord_flip() +
  theme(axis.title.y = element_blank()) +
  labs(title = "A")


econ_selected_plot <- mod_select %>%
  mutate(econ_model = !(fit_model == "em0_lm0")) %>%
  group_by(econ_model) %>%
  count() %>%
  ungroup() %>%
  mutate(p_selected = n / sum(n)) %>%
  # mutate(econ_model = fct_reorder(econ_model, p_selected)) %>%
  ggplot(aes(econ_model, p_selected)) +
  geom_col(color = "black", fill = "lightgrey") +
  scale_y_continuous(labels = percent, name = "Percentage Best Model") +
  labs(x = "Used Economics?", y = "Percentage Best Model",
       title = "B")

model_summary_plot <- model_selected_plot + econ_selected_plot


decision_tree = train(fit_model ~ .,
                      method = "rpart",
                      data =  prep_performance %>% select(-log_rmse),
                      tuneLength =40)

plot(decision_tree$finalModel)
text(decision_tree$finalModel, use.n=TRUE, all=TRUE, cex=0.5)
}


# compare to LIME

if (run_lime == TRUE){

  processed_sandbox <- readRDS(file = glue::glue("results/scrooge-results/{run_name}/processed_fisheries_sandbox.RDS"))

lime_fits <- pmap(list(
  fish = map(processed_sandbox$prepped_fishery,"fish"),
  fleet = map(processed_sandbox$prepped_fishery,"fleet"),
  data = map(processed_sandbox$prepped_fishery, "scrooge_data")
),
safely(fit_lime))

lime_worked <- map(lime_fits,"error") %>% map_lgl(is.null)

lime_fits <- map(lime_fits, "result")

saveRDS(object = lime_fits,file =  paste0(run_dir,"/lime_fits.RDS"))

processed_sandbox <- readRDS(file = glue::glue("results/scrooge-results/{run_name}/processed_fisheries_sandbox.RDS"))

lime_fits <- readRDS(file =  paste0(run_dir,"/lime_fits.RDS"))

lime_worked <- map_dbl(lime_fits, length)

lcomp_years <- map(processed_sandbox$prepped_fishery,c("scrooge_data","length_comps_years"))

processed_limes <- map2(lime_fits, lcomp_years, safely(process_lime))

saveRDS(object = processed_limes,file =  paste0(run_dir,"/processed_limes.RDS"))

}

if (run_lbspr == TRUE){

  # processed_sandbox <- readRDS(file = glue::glue("results/scrooge-results/{run_name}/processed_fisheries_sandbox.RDS"))


  experiments <- expand.grid(
    period = c("middle", "end"),
    window = c(15),
    economic_model = c(0),
    likelihood_model = c(0),
    prop_years_lcomp_data = c(1),
    fishery =   unique(fisheries_sandbox$fishery),
    stringsAsFactors = F
  ) %>%
    filter(
      !(economic_model == 2 & likelihood_model == 1),!(economic_model == 3 &
                                                         likelihood_model == 2)
    )

  sub_fisheries_sandbox <- experiments %>%
    left_join(fisheries_sandbox, by = "fishery") %>%
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


  lbspr_fits <- pmap(list(
    fish = map(sub_fisheries_sandbox$prepped_fishery, "fish"),
    fleet = map(sub_fisheries_sandbox$prepped_fishery, "fleet"),
    data = map(sub_fisheries_sandbox$prepped_fishery, "scrooge_data")
  ),
  safely(fit_lbspr))

lbspr_worked <- map(lbspr_fits,"error") %>% map_lgl(is_null)

lbspr_results <- map(lbspr_fits,"result") %>% keep(lbspr_worked)

fit <- lbspr_results[[1]]

truth <- sub_fisheries_sandbox$prepped_fishery[[1]]

process_lbspr <- function(fit,truth){

  ests <- fit@Ests %>%
    as_data_frame() %>%
    set_names(tolower)

  sampled_years <- truth$sampled_years

  true_f <- truth$simed_fishery %>%
    filter(year %in% sampled_years) %>%
    group_by(year) %>%
    summarise(true_f = unique(f))

  ests$true_fm <- true_f$true_f / truth$fish$m

  ests$true_sel50 <- truth$fleet$length_50_sel

  return(ests)

}

compare_lbspr <-  sub_fisheries_sandbox %>%
  filter(lbspr_worked) %>%
  mutate(lbspr_results = map2(lbspr_results, prepped_fishery, process_lbspr)) %>%
  select(experiment, rec_ac, sigma_r, economic_model, lbspr_results) %>%
  unnest() %>%
  group_by(experiment,rec_ac,sigma_r,economic_model) %>%
  summarise(rmse = sqrt(mean((true_fm - fm)^2))) %>%
  ungroup()

compare_lbspr %>%
  ggplot(aes(pmin(2,rmse))) +
  geom_histogram(aes(y = ..ncount..),fill = "lightgrey", color = "black") +
  labs(x = "rmse (units of F)", y = "Proportion")


compare_lbspr %>%
  ggplot(aes(sigma_r, pmin(2,rmse))) +
  geom_point()

compare_lbspr %>%
  ggplot(aes(rec_ac, pmin(2, rmse))) +
  geom_point()






}

# run diagnostics ---------------------------------------------------------


# make figures ------------------------------------------------------------

