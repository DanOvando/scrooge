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
library(tidyverse)
rstan::rstan_options(auto_write = FALSE)

functions <- list.files(here::here("functions"))

walk(functions, ~here::here("functions",.x) %>% source()) # load local functions

run_name <- "v1.0"

run_dir <- here::here("results", run_name)

if (dir.exists(run_dir) == F){ dir.create(run_dir, recursive = T)}

run_description <- "Most basic version of results up and running"

write(run_description, file = here::here("results",run_name,"description.txt"))

scrooge_theme <- theme_ipsum(base_size = 14, axis_title_size = 18)

theme_set(scrooge_theme)

# run options -------------------------------------------------------------

sim_fisheries <- T

fit_models <- T


# load data ---------------------------------------------------------------

cdfw_data <- read_csv(file = file.path('data','cfdw-catches.csv'))


has_timeseries <- cdfw_data %>%
  group_by(sci_name) %>%
  summarise(has_catch = min(year) <=2000 & max(year) == 2015) %>%
  filter(has_catch == T)

fill_catches <- function(min_year, max_year, catches){

  full_frame <- data_frame(year = min_year:max_year) %>%
    left_join(catches, by = 'year') %>%
    mutate(catch = zoo::na.approx(catch_lbs, rule = 2)) %>%
    select(-catch_lbs)

}

cdfw_catches <- cdfw_data %>%
  filter(sci_name %in% unique(has_timeseries$sci_name)) %>%
  group_by(sci_name,year) %>%
  summarise(catch_lbs = sum(pounds_caught)) %>%
  group_by(sci_name) %>%
  nest(-sci_name, .key = 'catches') %>%
  filter(is.na(sci_name) == F) %>%
  mutate(catches = map(catches, fill_catches, min_year = 2000, max_year = 2015)) %>%
  unnest()

cdfw_data$common_name %>% unique()

# create simulated fisheries ----------------------------------------------


if (sim_fisheries == T)
  {

fisheries_sandbox <-
  purrr::cross_df(list(
    sci_name = c("Atractoscion nobilis"),# "Sebastes mystinus"),
    fleet_model = c(
      "constant-catch",
      "constant-effort",
      "supplied-catch",
      "open-access"
    ),
    sigma_r = c(0,1),
    sigma_effort = c(0,1),
    price_cv = c(0,0.2),
    cost_cv = c(0, 0.2),
    price_ac = 0.25,
    cost_ac = 0.25
  ))

  fleet_model_params <- data_frame(fleet_model = c(
    "constant-catch",
    "constant-effort",
    "supplied-catch",
    "open-access"
  ),
fleet_params = list(
  list(target_catch = 1000),
  list(initial_effort = 100),
  list(catches = cdfw_catches$catch[cdfw_catches$sci_name == "semicossyphus pulcher"]),
  list(theta = 0.1, theta_tuner = 0.1)
))

  fisheries_sandbox <- fisheries_sandbox %>%
    left_join(fleet_model_params, by = "fleet_model") %>%
    # filter(fleet_model == "constant-effort") %>%
    # slice(1) %>%
    mutate(prepped_fishery = pmap(list(
      sci_name = sci_name,
      fleet_model = fleet_model,
      fleet_params = fleet_params,
      sigma_r = sigma_r,
      sigma_effort = sigma_effort,
      price_cv = price_cv,
      cost_cv = cost_cv,
      price_ac = price_ac,
      cost_ac = cost_ac
    ), prepare_fishery))

  save(file = here::here("results",run_name,"fisheries_sandbox.Rdata"), fisheries_sandbox)

} else{
  load(file = here::here("results",run_name,"fisheries_sandbox.Rdata"))
}



# oa <- fisheries_sandbox %>%
#   filter(fleet_model == "open-access") %>%
#   slice(10) %>%
#   mutate(summary_plot = map(prepped_fishery, plot_simmed_fishery)) %>%
#   mutate(scrooge_fit = map(prepped_fishery, ~fit_scrooge(data = .x$scrooge_data)))
#
# oa <- oa %>%
#   mutate(processed_scrooge = map(scrooge_fit, process_scrooge)) %>%
#   mutate(observed = map(prepped_fishery, "simed_fishery")) %>%
#   mutate(scrooge_performance = pmap(list(
#     observed = observed,
#     predicted = processed_scrooge
#   ), judge_performance))
#


# apply candidate assessment models ---------------------------------------

# test three simple cases to verify core model performance

# perfecto open access

pfo <- fisheries_sandbox %>%
  filter(fleet_model == "open-access", sigma_r == 0, sigma_effort == 0, price_cv == 0, cost_cv == 0) %>%
  slice(1) %>%
  mutate(summary_plot = map(prepped_fishery, plot_simmed_fishery)) %>%
  mutate(scrooge_fit = map(prepped_fishery, ~fit_scrooge(data = .x$scrooge_data, iter = 10000, warmup = 4000)))

pfo$summary_plot[[1]]

pfo$prepped_fishery[[1]]$length_comps %>%
  gather(age, numbers, -year) %>%
  mutate(age = as.numeric(age)) %>%
  group_by(year) %>%
  summarise(ncaught = sum(numbers)) %>%
  ggplot(aes(year, ncaught)) +
  geom_point()

pfo <- pfo %>%
  mutate(processed_scrooge = map(scrooge_fit, process_scrooge)) %>%
  mutate(observed = map(prepped_fishery, "simed_fishery")) %>%
  mutate(scrooge_performance = pmap(list(
    observed = observed,
    predicted = processed_scrooge
  ), judge_performance)) %>%
  mutate(rmse = map_dbl(scrooge_performance, ~.x$comparison_summary$rmse),
         bias =  map_dbl(scrooge_performance, ~.x$comparison_summary$bias)) %>%
  arrange(rmse)

pfo$scrooge_performance[[1]]$comparison_plot

# variable open access
#
vfo <- fisheries_sandbox %>%
  filter(
    fleet_model == "open-access",
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

vfo$summary_plot[[1]]

vfo$prepped_fishery[[1]]$length_comps %>%
  gather(age, numbers, -year) %>%
  mutate(age = as.numeric(age)) %>%
  group_by(year) %>%
  summarise(ncaught = sum(numbers)) %>%
  ggplot(aes(year, ncaught)) +
  geom_point()

vfo <- vfo %>%
  mutate(processed_scrooge = map(scrooge_fit, process_scrooge)) %>%
  mutate(observed = map(prepped_fishery, "simed_fishery")) %>%
  mutate(scrooge_rec_performance = pmap(list(
    observed = observed,
    predicted = processed_scrooge
  ), judge_performance, observed_variable = rec_dev, predicted_variable = "rec_dev_t")) %>%
  mutate(scrooge_performance = pmap(list(
    observed = observed,
    predicted = processed_scrooge
  ), judge_performance)) %>%
  mutate(rmse = map_dbl(scrooge_performance, ~.x$comparison_summary$rmse),
         bias =  map_dbl(scrooge_performance, ~.x$comparison_summary$bias)) %>%
  arrange(rmse)

vfo$scrooge_performance[[1]]$comparison_plot


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
  gather(age, numbers, -year) %>%
  mutate(age = as.numeric(age)) %>%
  group_by(year) %>%
  summarise(ncaught = sum(numbers)) %>%
  ggplot(aes(year, ncaught)) +
  geom_point()

cfo <- cfo %>%
  mutate(processed_scrooge = map(scrooge_fit, process_scrooge)) %>%
  mutate(observed = map(prepped_fishery, "simed_fishery")) %>%
  mutate(scrooge_rec_performance = pmap(list(
    observed = observed,
    predicted = processed_scrooge
  ), judge_performance, observed_variable = rec_dev, predicted_variable = "rec_dev_t")) %>%
  mutate(scrooge_performance = pmap(list(
    observed = observed,
    predicted = processed_scrooge
  ), judge_performance)) %>%
  mutate(rmse = map_dbl(scrooge_performance, ~.x$comparison_summary$rmse),
         bias =  map_dbl(scrooge_performance, ~.x$comparison_summary$bias)) %>%
  arrange(rmse)

cfo$scrooge_performance[[1]]$comparison_plot

cfo$scrooge_rec_performance[[1]]$comparison_plot


if (fit_models == T){

  sfs <- safely(fit_scrooge)

  fisheries_sandbox <- fisheries_sandbox %>%
  mutate(summary_plot = map(prepped_fishery, plot_simmed_fishery)) %>%
  mutate(scrooge_fit = map(prepped_fishery, ~sfs(data = .x$scrooge_data)))


  save(file = here::here("results",run_name,"fitted_fisheries_sandbox.Rdata"), fisheries_sandbox)


} else{

  load(file = here::here("results",run_name,"fitted_fisheries_sandbox.Rdata"))

}

# process fits

fisheries_sandbox <- fisheries_sandbox %>%
  mutate(scrooge_error = map(fisheries_sandbox$scrooge_fit,"error")) %>%
  mutate(scrooge_fit = map(fisheries_sandbox$scrooge_fit,"result")) %>%
  filter(map_lgl(scrooge_error, is.null)) %>%
  mutate(processed_scrooge = map(scrooge_fit, process_scrooge)) %>%
  mutate(observed = map(prepped_fishery, "simed_fishery")) %>%
  mutate(scrooge_performance = pmap(list(
    observed = observed,
    predicted = processed_scrooge
  ), judge_performance)) %>%
  mutate(rmse = map_dbl(scrooge_performance, ~.x$comparison_summary$rmse),
         bias =  map_dbl(scrooge_performance, ~.x$comparison_summary$bias)) %>%
  arange(rmse)


huh <- fisheries_sandbox2 %>%
  mutate(rmse = map_dbl(scrooge_performance, ~.x$comparison_summary$rmse)) %>%
  arrange(rmse)


# run diagnostics ---------------------------------------------------------


# make figures ------------------------------------------------------------

