fit_lbspr<- function(data,fish, fleet){
  pars <- new("LB_pars")

  pars@Species <- fish$scientific_name

  pars@Linf <- fish$linf

  pars@L50 <- fish$length_50_mature

  pars@L95 <-  fish$length_95_mature

  pars@MK <- fish$m / fish$vbk

  pars@BinWidth <- 1

  pars@L_units <- "cm"

  pars@Steepness <- fish$steepness

  lcomps <- new("LB_lengths")


  bin_mids <- data_frame(lbin = colnames(data$length_comps) %>% as.numeric()) %>%
    mutate(ubin = dplyr::lead(lbin))

  bin_mids$ubin[nrow(bin_mids)] <- bin_mids$lbin[nrow(bin_mids)] + 1

  bin_mids <- map2_dbl(bin_mids$lbin, bin_mids$ubin, ~mean(c(.x,.y)))

  lcomps@LMids <- bin_mids

  lcomps@LData <- t(as.matrix(data$length_comps))

  lcomps@L_units <- pars@L_units

  lcomps@NYears <- nrow(data$length_comps)

  lcomps@Years <- 1:nrow(data$length_comps)

  lbspr_fit <- LBSPR::LBSPRfit(pars, lcomps)

  return(lbspr_fit)

}
