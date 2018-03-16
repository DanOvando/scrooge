fit_scrooge <- function(data, scrooge_file = "scrooge.stan",
                        chains = 1, refresh = 25, cores = 1,
                        iter = 4000,
                        warmup = 2000,
                        adapt_delta = 0.9){

  fit <-
    rstan::stan(
      file = here::here("scripts", scrooge_file),
      data = data,
      chains = chains,
      refresh = refresh,
      cores = cores,
      iter = iter,
      warmup = warmup,
      control = list(
      adapt_delta = adapt_delta
    ))

}