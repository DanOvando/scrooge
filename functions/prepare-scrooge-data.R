prepare_scrooge_data <- function(fish,
                                 fleet,
                                 bias,
                                 obs_error,
                                 linf_buffer = 1.25,
                                 sim,
                                 sample_type = "catch",
                                 percent_sampled = 1,
                                 cv_effort,
                                 economic_model,
                                 beta = 2
                                 ) {


true_min <- fish$min_age

true_max <- fish$max_age

time_step <- fish$time_step

fish <- map_if(fish, is.numeric, ~.x * exp(rnorm(1, rnorm(1,0, bias), obs_error)))

fish$min_age <- true_min

fish$max_age <- true_max

fish$time_step <- time_step

fish$length_at_age <-   fish$linf * (1 - exp(- fish$vbk * (seq(fish$min_age,fish$max_age, by = fish$time_step) -  fish$t0)))

fish$weight_at_age <-  fish$weight_a * fish$length_at_age ^ fish$weight_b

fish$maturity_at_age <-
  ((1 / (1 + exp(-log(
    19
  ) * ((seq(fish$min_age, fish$max_age, by = fish$time_step) - fish$age_50_mature) / pmax(0.01,fish$age_95_mature - fish$age_50_mature)
  )))))

fleet <- map_if(fleet, is.numeric, ~.x * exp(rnorm(1, rnorm(1,0, bias), obs_error)))

fleet$length_95_sel <- fleet$length_50_sel + fleet$delta

fleet$sel_at_age <-
  ((1 / (1 + exp(-log(
    19
  ) * ((fish$length_at_age - fleet$length_50_sel) / pmax(0.01,fleet$length_95_sel - fleet$length_50_sel)
  )))))


length_at_age_key <- generate_length_at_age_key(
  min_age = fish$min_age,
  max_age = fish$max_age,
  cv = fish$cv_len,
  linf = fish$linf,
  k = fish$vbk,
  t0 = fish$t0,
  time_step = fish$time_step,
  linf_buffer = linf_buffer
) %>%
  ungroup() %>%
  select(age, length_bin, p_bin) %>%
  spread(length_bin, p_bin) %>%
  select(-age)


length_comps <- sim %>%
  select(year, patch, age, numbers, numbers_caught) %>%
  nest(-year,-patch, .key = n_at_age) %>%
  mutate(catch_length = map(
    n_at_age,
    ~ spasm::sample_lengths(
      n_at_age = .x,
      cv = fish$cv_len,
      k = fish$vbk,
      linf = fish$linf,
      t0 = fish$t0,
      sample_type = sample_type,
      percent_sampled = percent_sampled,
      time_step = fish$time_step,
      linf_buffer = linf_buffer
    )
  )) %>%
  select(year, catch_length) %>%
  unnest() %>%
  mutate(numbers = round(numbers)) %>%
  spread(length_bin, numbers) %>%
  mutate(year = year - min(year) + 1)

price_and_cost_history <- sim %>%
  group_by(year) %>%
  summarise(price = unique(price),
            cost = unique(cost),
            q = unique(unique(f) / unique(effort))) %>%
  gather(variable, value,-year) %>%
  ungroup() %>%
  mutate(value = value * exp(rnorm(nrow(.), rnorm(nrow(.),0, bias), obs_error))) %>%
  group_by(variable) %>%
  mutate(max_scaled_value = value / max(value)) %>%
  ungroup()

price_t <- price_and_cost_history %>%
  filter(variable == "price")

cost_t <- price_and_cost_history %>%
  filter(variable == "cost")

q_t <- price_and_cost_history %>%
  filter(variable == "q")

data_t <- sim %>%
  group_by(year) %>%
  summarise(profits = sum(profits),
            effort = unique(effort)) %>%
  mutate(ppue = profits / effort) %>%
  mutate(effort = effort * exp(rnorm(nrow(.), rnorm(nrow(.),0, bias), obs_error)),
         ppue = ppue * exp(rnorm(nrow(.), rnorm(nrow(.),0, bias), obs_error))) %>%
  ungroup() %>%
  mutate(lead_effort = lead(effort, default = 1)) %>%
  mutate(perc_change_effort = lead_effort/effort)


scrooge_data <- list(
  economic_model = economic_model,
  estimate_recruits = 1,
  length_comps = length_comps,
  observed_effort = data_t$effort,
  perc_change_effort = data_t$perc_change_effort,
  ppue_t = data_t$ppue,
  length_comps_years  = length_comps$year,
  price_t = price_t$value,
  cost_t = cost_t$max_scaled_value,
  q_t = q_t$max_scaled_value,
  beta = beta,
  length_50_sel_guess = fleet$length_50_sel * exp(rnorm(1, rnorm(1,0, bias), obs_error)),
  delta_guess = 2,
  n_lcomps = nrow(length_comps),
  nt = length(length_comps$year),
  n_burn = 100,
  n_ages = fish$max_age + 1,
  n_lbins = ncol(length_at_age_key),
  ages = 1:(fish$max_age + 1),
  m = fish$m,
  h = fish$steepness,
  r0 = fish$r0,
  k = fish$vbk,
  loo = fish$linf,
  t0 = fish$t0,
  length_at_age_key = as.matrix(length_at_age_key),
  mean_length_at_age = fish$length_at_age,
  mean_weight_at_age = fish$weight_at_age,
  mean_maturity_at_age = fish$maturity_at_age,
  cv_effort = cv_effort
)

return(scrooge_data)

}