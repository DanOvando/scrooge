process_scrooge <- function(fit, to_tidy = c("f_t", "n_tl", "c_t","rec_dev_t")) {
  fits <- rstan::extract(fit)

  extract_thing <- function(thing, fits) {

     tidy_thing <- purrr::array_branch(purrr::pluck(fits, thing),1)

    tidy_thing <-
      data_frame(thing = map(tidy_thing, ~ as_data_frame(.x))) %>%
      mutate(iteration = 1:nrow(.))

    times <- nrow(tidy_thing$thing[[1]])

    tidy_thing <- tidy_thing %>%
      unnest() %>%
      mutate(year = rep(1:times, length(unique(.$iteration))))
  }

  if (any(str_detect(to_tidy, "all"))) {
    to_tidy <- names(fits)

  }

  tidy_fits <- map(to_tidy, extract_thing, fits = fits)

  names(tidy_fits) <- to_tidy

  return(tidy_fits)

}