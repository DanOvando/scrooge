est_msy <- function(fmsy, data, time, use = "fit") {

  nframe <- matrix(NA, nrow = time, ncol = data$n_ages)

  nframe[1, ] <- data$r0 * exp(-data$m * 0:(data$n_ages - 1))

  ssb0 <-
    sum(nframe[1, ] * data$mean_weight_at_age * data$mean_maturity_at_age)

  ssbtemp <- NA

  selguess <-
    as.numeric(data$mean_length_at_age > data$length_50_sel_guess)

  biomass <- rep(0, time)

  for (i in 2:time) {
    ssbtemp = sum(nframe[i - 1, ] * data$mean_weight_at_age * data$mean_maturity_at_age)

    recs <-
      ((0.8 * data$r0 * data$h * ssbtemp) / (0.2 * ssb0 * (1 - data$h) + (data$h - 0.2) * ssbtemp))

    nframe[i, 1] <- recs

    nframe[i, 2:data$n_ages] = nframe[i - 1, 1:(data$n_ages - 1)] * exp(-(data$m + fmsy * selguess[1:(data$n_ages - 1)]))

    nframe[i, data$n_ages] <-
      nframe[i, data$n_ages] + nframe[i - 1, data$n_ages] * exp(-(data$m + fmsy * selguess[data$n_ages]))

  }

  b_msy <- sum(nframe[time, ] * data$mean_weight_at_age)


  msy <-
    sum(((fmsy * selguess) / (data$m + fmsy * selguess)) * (nframe[time,] * (1 - exp(
      -(data$m + fmsy * selguess)
    )) * data$mean_weight_at_age))

  if (use == "fit") {
    out <- -msy
  } else{
    out <- b_msy
  }

  return(out)

}