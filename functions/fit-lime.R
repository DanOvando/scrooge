fit_lime <- function(data,fish, fleet){

  scrooge_lh <- LIME::create_lh_list(vbk= fish$vbk,
                               linf= fish$linf,
                               t0= fish$t0,
                               lwa= fish$weight_a,
                               lwb= fish$weight_b,
                               S50= fleet$length_50_sel,
                               S95=fleet$length_95_sel,
                               selex_input="length",
                               selex_type=c("logistic"),
                               M50= fish$length_50_mature,
                               M95= fish$length_95_mature,
                               maturity_input="length",
                               M= fish$m,
                               binwidth=1,
                               CVlen= fish$cv_len,
                               SigmaR= fish$sigma_r + .001,
                               SigmaF= 0.1,
                               SigmaC=0.2,
                               SigmaI=0.2,
                               R0= fish$r0,
                               qcoef= fleet$q,
                               start_ages=0,
                               rho=0,
                               nseasons=1)

  temp_LF_matrix <- data$length_comps


  LF_matrix <- temp_LF_matrix %>%
    as.matrix()

  rownames(LF_matrix) <- 1:nrow(LF_matrix)


  scrooge_data_LF <-
    list("years" = 1:nrow(LF_matrix), "LF" = LF_matrix)


  res <- LIME::run_LIME(modpath=NULL,
                  lh=scrooge_lh,
                  input_data=scrooge_data_LF,
                  est_sigma="log_sigma_R",
                  data_avail="LC",
                  newtonsteps=3)

}