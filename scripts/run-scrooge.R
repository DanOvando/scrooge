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
library(tidyverse)
library(extrafont)
rstan::rstan_options(auto_write = TRUE)
extrafont::loadfonts()
functions <- list.files(here::here("functions"))

walk(functions, ~ here::here("functions", .x) %>% source()) # load local functions

run_name <- "v1.0"

run_dir <- here::here("results", run_name)

if (dir.exists(run_dir) == F) {
  dir.create(run_dir, recursive = T)
}

run_description <- "Most basic version of results up and running"

write(run_description,
      file = here::here("results", run_name, "description.txt"))

scrooge_theme <- theme_ipsum(base_size = 14, axis_title_size = 18)

theme_set(scrooge_theme)

# run options -------------------------------------------------------------

sim_fisheries <- T

fit_models <- F

run_tests <- F

in_clouds <-  F

n_cores <- 1

# load data ---------------------------------------------------------------

if (in_clouds == F){

  data_dir <- "data"
} else{

  data_dir <- "data/google-bucket"
}

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

cdfw_data$common_name %>% unique()

# create simulated fisheries ----------------------------------------------


if (sim_fisheries == T)
{
  fisheries_sandbox <-
    purrr::cross_df(
      list(
        sci_name = c("Atractoscion nobilis"), # "Scomber japonicus"),
        fleet_model = c(
          # "constant-catch",
          "constant-effort",
          # "supplied-catch",
          "open-access"
        ),
        sigma_r = c(0.1),
        sigma_effort = c(0, 0.1),
        price_cv = c(0, 0.5),
        cost_cv = c(0),
        price_ac = 0,
        cost_ac = 0,
        economic_model = c(1,0)
      )
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
      list(initial_effort = 200),
      list(catches = cdfw_catches$catch[cdfw_catches$sci_name == "semicossyphus pulcher"]),
      list(theta = 0.1, theta_tuner = 0.05, initial_effort = 200)
    )
  )

  fisheries_sandbox <- fisheries_sandbox %>%
    left_join(fleet_model_params, by = "fleet_model") %>%
    # filter(fleet_model == "constant-effort") %>%
    # slice(1) %>%
    mutate(prepped_fishery = pmap(
      list(
        sci_name = sci_name,
        fleet_model = fleet_model,
        fleet_params = fleet_params,
        sigma_r = sigma_r,
        sigma_effort = sigma_effort,
        price_cv = price_cv,
        cost_cv = cost_cv,
        price_ac = price_ac,
        cost_ac = cost_ac
      ),
      prepare_fishery,
      sim_years = 20,
      burn_years = 50,
      price = 0.15
    )) %>%
    mutate(summary_plot = map(prepped_fishery, plot_simmed_fishery))


  save(file = here::here("results", run_name, "fisheries_sandbox.Rdata"),
       fisheries_sandbox)

} else{
  load(file = here::here("results", run_name, "fisheries_sandbox.Rdata"))
}

# apply candidate assessment models ---------------------------------------

# test three simple cases to verify core model performance

if (run_tests == T) {
  # perfecto open access

  pfo <- fisheries_sandbox %>%
    filter(fleet_model == "open-access",
           sigma_r == 0,
           sigma_effort == 0,
           price_cv == 0,
           cost_cv == 0) %>%
    slice(1) %>%
    mutate(summary_plot = map(prepped_fishery, plot_simmed_fishery)) %>%
    mutate(scrooge_fit = map(
      prepped_fishery,
      ~ fit_scrooge(
        data = .x$scrooge_data,
        iter = 2000,
        warmup = 1000,
        adapt_delta = 0.8,
        economic_model = 1,
        scrooge_file = "scrooge_v2.0"
      )
    ))

  fit <-
    rstan::stan(
      fit = pfo$scrooge_fit[[1]],
      data = pfo$prepped_fishery[[1]]$scrooge_data,
      chains = 1,
      refresh = 25,
      cores = 1,
      iter = 2000,
      warmup = 1000,
      control = list(
        adapt_delta = 0.8)    )
  pfo$scrooge_fit[[1]] %>% class()

  pfo$summary_plot[[1]]

  pfo$prepped_fishery[[1]]$length_comps %>%
    gather(age, numbers,-year) %>%
    mutate(age = as.numeric(age)) %>%
    group_by(year) %>%
    summarise(ncaught = sum(numbers)) %>%
    ggplot(aes(year, ncaught)) +
    geom_point()

  pfo <- pfo %>%
    mutate(processed_scrooge = map(scrooge_fit, process_scrooge)) %>%
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
    ) %>%
    mutate(
      rmse = map_dbl(scrooge_performance, ~ .x$comparison_summary$rmse),
      bias =  map_dbl(scrooge_performance, ~ .x$comparison_summary$bias)
    ) %>%
    arrange(rmse)

  pfo$scrooge_performance[[1]]$comparison_plot +
    lims(y = c(0,2))

  pfo$scrooge_rec_performance[[1]]$comparison_plot


  # variable open access
  #
  vfo <- fisheries_sandbox %>%
    filter(
      fleet_model == "open-access",
      sigma_r == max(sigma_r),
      sigma_effort == min(sigma_effort),
      price_cv == max(price_cv),
      price_ac == max(price_ac),
      cost_ac == max(cost_ac),
      cost_cv == max(cost_cv)
    ) %>%
    slice(1) %>%
    mutate(summary_plot = map(prepped_fishery, plot_simmed_fishery)) %>%
    mutate(scrooge_fit = map(
      prepped_fishery,
      ~ fit_scrooge(
        data = .x$scrooge_data,
        iter = 8000,
        warmup = 4000,
        scrooge_file = "scrooge_v3.0",
        adapt_delta = 0.8
      )
    ))

  vfo$summary_plot[[1]]

  vfo$prepped_fishery[[1]]$length_comps %>%
    gather(age, numbers,-year) %>%
    mutate(age = as.numeric(age)) %>%
    group_by(year) %>%
    summarise(ncaught = sum(numbers)) %>%
    ggplot(aes(year, ncaught)) +
    geom_point()

  vfo <- vfo %>%
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

    vfo$scrooge_rec_performance[[1]]$comparison_plot

    vfo$scrooge_performance[[1]]$comparison_plot +
      lims(y = c(0,1))


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



  sfs <- safely(fit_scrooge)

  doParallel::registerDoParallel(cores = n_cores)
  foreach::getDoParWorkers()

  fits <- foreach::foreach(i = 1:nrow(fisheries_sandbox)) %dopar% {
    sfs(
      data = fisheries_sandbox$prepped_fishery[[i]]$scrooge_data,
      economic_model = fisheries_sandbox$economic_model[[i]],
      scrooge_file = "scrooge_v2.0"
    )

  } # close fitting loop


  fisheries_sandbox$scrooge_fit <- fits

  # fisheries_sandbox <- fisheries_sandbox %>%
  #   slice(1) %>%
  #   mutate(summary_plot = map(prepped_fishery, plot_simmed_fishery)) %>%
  #   mutate(scrooge_fit = map2(
  #     prepped_fishery,
  #     economic_model,
  #     ~ sfs(data = .x$scrooge_data, economic_model = .y,scrooge_file = "scrooge_v2.0")
  #   ))


  save(
    file = here::here("results", run_name, "fitted_fisheries_sandbox.Rdata"),
    fisheries_sandbox
  )


} else{
  load(file = here::here("results", run_name, "fitted_fisheries_sandbox.Rdata"))

}

# process fits

fisheries_sandbox <- fisheries_sandbox %>%
  mutate(scrooge_error = map(fisheries_sandbox$scrooge_fit, "error")) %>%
  mutate(scrooge_fit = map(fisheries_sandbox$scrooge_fit, "result")) %>%
  filter(map_lgl(scrooge_error, is.null)) %>%
  mutate(processed_scrooge = map(scrooge_fit, process_scrooge)) %>%
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
  ) %>%
  mutate(
    rmse = map_dbl(scrooge_performance, ~ .x$comparison_summary$rmse),
    bias =  map_dbl(scrooge_performance, ~ .x$comparison_summary$bias)
  ) %>%
  arrange(rmse)


oa_check <- fisheries_sandbox %>%
  filter(fleet_model == "open-access")

huh <- fisheries_sandbox2 %>%
  mutate(rmse = map_dbl(scrooge_performance, ~ .x$comparison_summary$rmse)) %>%
  arrange(rmse)


# run diagnostics ---------------------------------------------------------


# make figures ------------------------------------------------------------
