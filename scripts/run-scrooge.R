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
library(tidyverse)
rstan::rstan_options(auto_write = TRUE)
extrafont::loadfonts()
functions <- list.files(here::here("functions"))

walk(functions, ~ here::here("functions", .x) %>% source()) # load local functions

run_name <- "v1.2"

in_clouds <-  T

local_data <- T

if (in_clouds == F){

  run_dir <- here::here("results", run_name)

} else{

  run_dir <- here::here("results","scrooge-results", run_name)

}


if (dir.exists(run_dir) == F) {
  dir.create(run_dir, recursive = T)
}

run_description <- "Most basic version of results up and running"


scrooge_theme <- theme_ipsum(base_size = 14, axis_title_size = 18)

theme_set(scrooge_theme)

# run options -------------------------------------------------------------

sim_fisheries <- F

fit_models <- T

run_tests <- F



n_cores <- 5

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
        sci_name = c("Atractoscion nobilis", "Scomber japonicus","Lutjanus campechanus"),
        fleet_model = c(
          # "constant-catch",
          "constant-effort",
          # "supplied-catch",
          "open-access"
        ),
        sigma_r = c(0.1,0.4),
        sigma_effort = c(0, 0.1),
        price_cv = c(0, 0.5),
        cost_cv = c(0,0.75),
        price_ac = 0.75,
        cost_ac = 0.75,
        # economic_model = c(1,0),
        steepness = c(0.6,0.9),
        obs_error = c(0,0.2)
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
      list(theta = 0.1, theta_tuner = 0.75, initial_effort = 10)
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
        cost_ac = cost_ac,
        steepness = steepness,
        obs_error = obs_error
      ),
      safely(prepare_fishery),
      sim_years = 50,
      burn_years = 25,
      price = 2,
      cost = 5
    ))


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
}

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
    mutate(scrooge_fit = map2(
      prepped_fishery,
      ~ fit_scrooge(
        data = .x$scrooge_data,
        iter = 2000,
        warmup = 1000,
        adapt_delta = 0.8,
        economic_model = 1,
        scrooge_file = "bioeconomic_scrooge",
        in_clouds = in_clouds
      )
    ))


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
    mutate(lcomps = map2(prepped_fishery, processed_scrooge, ~process_lcomps(.x$length_comps, .y$n_tl))) %>%
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
      cost_cv == max(cost_cv)
    ) %>%
    slice(1) %>%
    mutate(prepped_fishery = map(prepped_fishery, subsample_data, window = 10, period = "end")) %>%
    mutate(scrooge_fit = map(
      prepped_fishery,
      ~ fit_scrooge(
        data = .x$scrooge_data,
        iter = 2000,
        warmup = 1000,
        scrooge_file = "bioeconomic_scrooge",
        adapt_delta = 0.8,
        max_treedepth = 12
      )
    ))


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

    b = vfo$scrooge_performance[[1]]$comparison_plot +
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


  experiments <- expand.grid(period = c("beginning", "middle","end"), window = c(2,5,10),
                             experiment = 1:nrow(fisheries_sandbox), stringsAsFactors = F)

  fisheries_sandbox <- fisheries_sandbox %>%
    mutate(experiment = 1:nrow(.)) %>%
    left_join(experiments, by = "experiment") %>%
    slice(sample(1:nrow(fisheries_sandbox),4, replace = F)) %>%
    mutate(prepped_fishery = pmap(list(
      prepped_fishery = prepped_fishery,
      window = window,
      period = period), subsample_data))

  # fit lbspr

  # fisheries_sandbox <- fisheries_sandbox %>%
  #   mutate(lbspr_fit = pmap(list(
  #     data = map(prepped_fishery, "scrooge_data"),
  #     fish = map(prepped_fishery, "fish"),
  #     fleet = map(prepped_fishery, "fleet")
  #   ), fit_lbspr))

  sfs <- safely(fit_scrooge)

  doParallel::registerDoParallel(cores = n_cores)

  foreach::getDoParWorkers()

  fits <- foreach::foreach(i = 1:nrow(fisheries_sandbox)) %dopar% {
    out <- sfs(
      data = fisheries_sandbox$prepped_fishery[[i]]$scrooge_data,
      fish = fisheries_sandbox$prepped_fishery[[i]]$fish,
      fleet = fisheries_sandbox$prepped_fishery[[i]]$fleet,
      experiment = fisheries_sandbox$experiment[i],
      scrooge_file = "bioeconomic_scrooge",
      iter = 2000,
      warmup = 1000,
      adapt_delta = 0.8,
      max_treedepth = 12,
      in_clouds = in_clouds,
      cloud_dir = cloud_dir
    )

  } # close fitting loop

  fisheries_sandbox$scrooge_fit <- fits



  save(
    file = glue::glue("{run_dir}/fitted_fisheries_sandbox.Rdata"),
    fisheries_sandbox
  )


} else{

  load(file = glue::glue("{run_dir}/fitted_fisheries_sandbox.Rdata"))

}


# process_fits ------------------------------------------------------------


scrooge_worked <- map(fisheries_sandbox$scrooge_fit,'error') %>% map_lgl(is.null)

fisheries_sandbox <- fisheries_sandbox %>%
  filter(scrooge_worked) %>%
  mutate(scrooge_fit = map(scrooge_fit,"result"))


if(in_clouds == T){

  loadfoo <- function(experiment, cloud_dir){
    readRDS(glue::glue("{cloud_dir}/{experiment}"))
  }

  fisheries_sandbox <- fisheries_sandbox %>%
    mutate(scrooge_fit = map(scrooge_fit, loadfoo, cloud_dir = cloud_dir))

}

fisheries_sandbox <-  fisheries_sandbox%>%
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
  mutate(lcomps = map2(prepped_fishery, processed_scrooge, ~process_lcomps(.x$length_comps, .y$n_tl))) %>%
  mutate(
    rmse = map_dbl(scrooge_performance, ~ .x$comparison_summary$rmse),
    bias =  map_dbl(scrooge_performance, ~ .x$comparison_summary$bias)
  ) %>%
  arrange(rmse)


# run diagnostics ---------------------------------------------------------


# make figures ------------------------------------------------------------
