tune_costs <-
  function(cost,
           msy,
           e_msy,
           b_msy,
           data,
           time,
           p_response,
           price,
           b_v_bmsy_target = 0.5,
           use = "fit") {



    p_msy <-
     price * msy - cost * e_msy ^ data$beta

    nframe <- matrix(NA, nrow = time, ncol = data$n_ages)

    nframe[1,] <- data$r0 * exp(-data$m * 0:(data$n_ages - 1))

    nframe[1, data$n_ages]  <- nframe[1, data$n_ages] * exp(-(data$m)) / (1 - exp(-(data$m)))

    ssb0 <-
      sum(nframe[1,] * data$mean_weight_at_age * data$mean_maturity_at_age)

    ssbtemp <- NA

    selguess <-
      as.numeric(data$mean_length_at_age > data$length_50_sel_guess)

    mean_q <- mean(data$q_t$value)

    effort <- rep(0.001, time)

    f <- effort * mean_q

    profits <- rep(0, time)

    catches <- rep(0, time)

    biomass <- rep(0, time)

    biomass[1] <- sum(nframe[1,] * data$mean_weight_at_age)

    ssb <- rep(0, time)

    recs <- rep(0, time)


    catches[1] <-
      sum(((f[1] * selguess) / (data$m + f[1] * selguess)) * (nframe[1, ] * (1 - exp(
        -(data$m + f[1] * selguess)
      )) * data$mean_weight_at_age))


    profits[1] <-
      price * catches[1] - cost * effort[1] ^ data$beta

   # wtf <-
   #    price * catches[1] - 0:1000 * effort[1] ^ data$beta

    for (i in 2:time) {
      new_effort <-
        effort[i - 1] + e_msy * (p_response * (profits[i - 1] / p_msy))

      if (new_effort <= 0) {
        new_effort = -.01 / (new_effort - 1)


      }

      effort[i] <- new_effort

      f[i] <- new_effort * mean_q

      ssb[i] = sum(nframe[i - 1,] * data$mean_weight_at_age * data$mean_maturity_at_age)

      recs[i] <-
        ((0.8 * data$r0 * data$h *   ssb[i] ) / (0.2 * ssb0 * (1 - data$h) + (data$h - 0.2) *   ssb[i]))

      nframe[i, 1] <- recs[i]

      nframe[i, 2:data$n_ages] = nframe[i - 1, 1:(data$n_ages - 1)] * exp(-(data$m + f[i] * selguess[1:(data$n_ages - 1)]))

      nframe[i, data$n_ages] <-
        nframe[i, data$n_ages] + nframe[i - 1, data$n_ages] * exp(-(data$m + f[i] * selguess[data$n_ages]))

      biomass[i] = sum(nframe[i, ] * data$mean_weight_at_age)

      catches[i] <-
        sum(((f[i] * selguess) / (data$m + f[i] * selguess)) * (nframe[i, ] * (1 - exp(
          -(data$m + f[i] * selguess)
        )) * data$mean_weight_at_age))


      profits[i] <-
        price * catches[i] - cost * effort[i] ^ data$beta


    }

    if (use == "fit"){
    out <- ((biomass[i] / b_msy) - b_v_bmsy_target) ^ 2
    } else {
      out <- biomass
    }


    return(out)

  }