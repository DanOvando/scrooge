lime_inputs <- readRDS("~/PhD/scrooge/data/lime_inputs.RDS")

limecomps <- lime_inputs$LF %>% drop()

limedata <- data

limedata$length_comps <- limecomps

limedata$length_comps_years <- 1:nrow(limecomps)

limedata$n_lcomps <- nrow(limecomps)

limedata$nt <- nrow(limecomps)

limedata$n_ages <- length(lime_inputs$ages)

limedata$n_lbins <- ncol(limecomps)

limedata$ages <- lime_inputs$ages

limedata$m <- lime_inputs$M

limedata$h <- lime_inputs$h

limedata$r0 <- lime_inputs$R0

limedata$k <- lime_inputs$vbk

limedata$loo <- lime_inputs$linf

limedata$t0 <- lime_inputs$t0

limedata$mean_length_at_age <- lime_inputs$L_a

limedata$mean_weight_at_age <- lime_inputs$W_a

limedata$mean_maturity_at_age <- lime_inputs$Mat_a

length_at_age_key <- generate_length_at_age_key(
  min_age = min(lime_inputs$ages),
  max_age = max(lime_inputs$ages),
  cv = lime_inputs$CVlen,
  linf = lime_inputs$linf,
  k = lime_inputs$vbk,
  t0 = lime_inputs$t0,
  time_step = 1,
  linf_buffer = 1.5
) %>%
  ungroup() %>%
  select(age, length_bin, p_bin) %>%
  spread(length_bin, p_bin) %>%
  select(-age)

length_at_age_key <- length_at_age_key[,1:limedata$n_lbins]

length_at_age_key <-  length_at_age_key / rowSums(length_at_age_key)

limedata$length_at_age_key <- length_at_age_key

data$length_50_sel_guess <-  lime_inputs$SL50

limedata$economic_model <- 0

limedata$sigma_effort_guess <- lime_inputs$SigmaF

length_at_age_key %>%
  mutate(age = 1:nrow(.)) %>%
  gather(length, prob, -age) %>%
  mutate(length = as.numeric(length)) %>%
  ggplot(aes(length, prob, color = age, group = age)) +
  geom_line()

limedata$length_comps %>%
  as_data_frame() %>%
  mutate(year = 1:nrow(.)) %>%
  gather(lbin, number,-year) %>%
  mutate(lbin = as.numeric(lbin)) %>%
  ggplot(aes(lbin, number)) +
  geom_col() +
  facet_wrap(~year)

saveRDS(limedata, "limedata.RDS")

