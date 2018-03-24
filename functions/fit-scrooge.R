fit_scrooge <- function(data, scrooge_file = "scrooge.stan",
                        chains = 1, refresh = 25, cores = 1,
                        iter = 4000,
                        warmup = 2000,
                        adapt_delta = 0.9, economic_model = NA){

  if (is.na(economic_model) == F){
  data <- purrr::list_modify(data,economic_model = economic_model)
  }

  fit <-
    rstan::stan(
      file = here::here("src", scrooge_file),
      data = data,
      chains = chains,
      refresh = refresh,
      cores = cores,
      iter = iter,
      warmup = warmup,
      control = list(
      adapt_delta = adapt_delta)
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

  return(fit)

}