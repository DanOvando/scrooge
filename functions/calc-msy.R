calc_msy <- function(fmsy,data, time){

  data_frame(n = rep(0, time), catch = rep(0,time))

  nframe <- matrix(NA, nrow = time, ncol = data$n_ages)

  nframe[1,] <- data$r0 * exp(-data$m * 0:(data$n_ages - 1))

  ssb0 <- sum(nframe[1,] * data$mean_weight_at_age * data$mean_maturity_at_age)

  ssbtemp <- NA

  selguess <- as.numeric(data$mean_length_at_age > data$length_50_sel_guess)

  for (i in 2:10){

    ssbtemp = sum(nframe[i - 1,] * data$mean_weight_at_age * data$mean_maturity_at_age)

    recs <- ((0.8 * data$r0 * data$h * ssbtemp) / (0.2 * ssb0 * (1 - data$h) + (data$h - 0.2) * ssbtemp))

    nframe[i,1] <- recs

    nframe[i,2:data$n_ages] = nframe[i - 1,1:(data$n_ages - 1)] * exp(-(data$m + fmsy * selguess[1:(data$n_ages - 1)]))

    nframe[i, data$n_ages] <- nframe[i, data$n_ages] + nframe[i - 1, data$n_ages] * exp(-(data$m + fmsy * selguess[data$n_ages]))

  }

  msy <- sum(((fmsy * selguess) / (data$m + fmsy * selguess)) *( nframe[time,] * (1 - exp(-(data$m + fmsy * selguess)))))

  return(-msy)

}