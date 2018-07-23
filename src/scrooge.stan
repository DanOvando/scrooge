/*scrooge_v4.0


*/


data{

//// run options ////

int<lower = 1> nt; // number of time steps to run

int<lower = 1> n_burn; // number of burn-in time steps

int n_ages; // number of ages

int n_lbins; // number of length bins;

int n_lcomps; // number of timesteps of length comps available

vector[n_ages] ages; // vector of ages

int<lower = 0> economic_model; // 0 = no bioeconomic prior, 1 = bioeconomic prior,

//// length data ////

int length_comps[n_lcomps,n_lbins];

int length_comps_years[n_lcomps]; // time steps in which length comps are available

//// economic data ////

row_vector[nt] perc_change_effort;

row_vector[nt] price_t;

row_vector[nt] relative_cost_t;

real max_cost_guess;

row_vector[nt] q_t;

vector[nt] ppue_t;

real beta;

real length_50_sel_guess;

real delta_guess;

real<lower = 0> p_response_guess;

real<lower = 0> cv_effort;

real<lower = 0> sigma_f; // penalty term for year to year variation in F, DO NOT ESTIMATE

//// biology ////

vector<lower=0>[n_lbins] bin_mids;

vector<lower=0>[n_ages] mean_length_at_age;

vector<lower=0>[n_ages] mean_weight_at_age;

vector<lower=0, upper=1>[n_ages] mean_maturity_at_age;

matrix[n_ages, n_lbins] length_at_age_key;

int age_sel; // estimate of the age at first selectivity

real m; //natural mortality

real h; //steepness

real r0; // virgin recruitment

real k; //vbk growth

real loo; // vbk linf

real t0; // t0

real sigma_r_guess;

real sd_sigma_r;


}

transformed data{

int n_total;

n_total = nt + n_burn;


}

parameters{

vector<lower = 0>[nt]  f_t; // f in time t

vector[nt + age_sel]  uc_rec_dev_t; //  recruitment deviates

real <lower = 0> sigma_r; // standard deviation of recruitment deviates

real<lower = 0> sigma_obs;

real  log_p_response;

real log_max_cost; // mean cost

real<lower = 0> p_length_50_sel; // length at 50% selectivity

} // close parameters block

transformed parameters{

  real ssb0;

  real ssb_temp;

  real sel_delta;

  real length_50_sel;

  real plus_group;

  real p_response;

  row_vector[n_ages] temp_a;

  row_vector[n_ages - 1] temp_a2;

  vector[nt]  effort_t; //  base effort multiplier

  vector[nt + age_sel]  rec_dev_t; //  base effort multiplier

  real temp_rec; // place holder for temporary recruitment deviates

  matrix[n_total, n_ages] n_ta; // numbers at time and age

  matrix[n_total, n_ages] ssb_ta; // spawning stock biomass at time and age

  matrix[n_total, n_ages] c_ta; // catch (biomass) at time and age

  matrix[n_total, n_ages] cn_ta; // catch (numbers) at time and age

  vector[n_total] c_t; // total catch at time step

  matrix[nt, n_lbins] p_lbin_sampled; // numbers at time and length bin

  row_vector[nt] cost_t;

  vector[nt] profit_t; // profits

  vector[nt] ppue_hat_t; // profit per unit effort

  vector[n_lbins] selectivity_at_bin; // mean selectivity at length bin midpoint

  vector[n_ages] mean_selectivity_at_age; // mean selectivity at age

  row_vector[n_ages] p_age_sampled;

  vector[nt - 1] delta_f;

  vector[nt - 1] ppue_hat;

  vector[nt - 1] perc_change_effort_hat;

  real temp_f;

  //////////////////////////
  //  set things up       //
  //////////////////////////

  p_response = exp(log_p_response);

  cost_t = exp(log_max_cost) * relative_cost_t;

  length_50_sel = loo * p_length_50_sel;

  sel_delta = 2;

  rec_dev_t = exp(uc_rec_dev_t);

  n_ta = rep_matrix(0,n_total, n_ages);

  ssb_ta = rep_matrix(0,n_total, n_ages);

  c_ta = rep_matrix(0,n_total, n_ages);

  cn_ta = rep_matrix(0,n_total, n_ages);

  p_lbin_sampled = rep_matrix(0,nt, n_lbins);

  selectivity_at_bin = 1.0 ./ (1 + exp(-log(19) * ((bin_mids - length_50_sel) / sel_delta))); // selectivity ogive at age

  mean_selectivity_at_age = length_at_age_key * selectivity_at_bin; // calculate mean selectivity at age given variance in length at age

  //////////////////////////
  // run population model //
  //////////////////////////

  effort_t[1] = f_t[1] / q_t[1];

  n_ta[1,1:n_ages] = r0 * exp(-m * (ages - 1))';

  plus_group = n_ta[1, n_ages - 1] * exp(-m) / (1 - exp(-m));  //plus group

  n_ta[1, n_ages] = plus_group;

  temp_a = n_ta[1, 1:n_ages] .* mean_maturity_at_age'.* mean_weight_at_age';

  ssb_ta[1, 1:n_ages] = temp_a; // calculate ssb at age

  ssb0 = sum(ssb_ta[1, 1:n_ages]);


  temp_rec = 1;

  temp_f = f_t[1];

  // order of events: spawn, grow and die, recruit

  for (t in 2:n_total){

  if (t > (n_burn - age_sel)){

   temp_rec = rec_dev_t[t - n_burn + age_sel];

  } // add in recruitment including some in the burn-in period if desired

  if (t > n_burn){

    temp_f = f_t[t - n_burn];

  }

    ssb_temp =  sum(ssb_ta[t - 1, 1:n_ages]);

    n_ta[t,1] = ((0.8 * r0 * h *ssb_temp) / (0.2 * ssb0 * (1 - h) + (h - 0.2) * ssb_temp)) * temp_rec; //calculate recruitment

    temp_a2 = n_ta[t - 1, 1:(n_ages -1)] .* (exp(-(m + temp_f * mean_selectivity_at_age[1:(n_ages - 1)])))';

    n_ta[t , 2:n_ages] = temp_a2; // grow and die

    n_ta[t, n_ages] = n_ta[t, n_ages] + n_ta[t - 1, n_ages] .* (exp(-(m + temp_f * mean_selectivity_at_age[n_ages])))'; // fill in plus group

  ssb_ta[t, 1:n_ages] = n_ta[t, 1:n_ages] .* mean_maturity_at_age'.* mean_weight_at_age'; // calculate ssb at age

   cn_ta[t - 1, 1:n_ages] = ((temp_f * mean_selectivity_at_age) ./ (m + temp_f * mean_selectivity_at_age))' .* n_ta[t - 1, 1:n_ages] .* (1 - exp(-(m + temp_f * mean_selectivity_at_age)))'; // .* mean_weight_at_age';

   c_ta[t - 1, 1:n_ages] = cn_ta[t - 1, 1:n_ages] .* mean_weight_at_age';

   c_t[t - 1] = sum(c_ta[t-1, 1:n_ages]);

  // calculate economic parameters

  if (t > n_burn + 1) {

  profit_t[t - 1 - n_burn] = price_t[t - 1 - n_burn] * c_t[t - 1] - cost_t[t - 1 - n_burn] * effort_t[t - 1 - n_burn] ^ beta;

  ppue_hat_t[t - 1 - n_burn] = profit_t[t - 1 - n_burn] / effort_t[t - 1 - n_burn];

  effort_t[t - n_burn] = f_t[t - n_burn] / q_t[t - n_burn];

  // sample lengths //

  p_age_sampled = cn_ta[t - 1, 1:n_ages] / sum(cn_ta[t - 1, 1:n_ages]);

  p_lbin_sampled[t - 1 - n_burn, 1:n_lbins] = (p_age_sampled * length_at_age_key) / sum(p_age_sampled * length_at_age_key);  // calculate the proportion of the sampled catch in each length bin

  }

  } // close time loop


// fill in final time step
    cn_ta[n_total, 1:n_ages] = ((f_t[nt] * mean_selectivity_at_age) ./ (m + f_t[nt] * mean_selectivity_at_age))' .* n_ta[n_total, 1:n_ages] .* (1 - exp(-(m + f_t[nt] * mean_selectivity_at_age)))'; //

   c_ta[n_total, 1:n_ages] = cn_ta[n_total, 1:n_ages] .* mean_weight_at_age';

   c_t[n_total] = sum(c_ta[n_total, 1:n_ages]);

  profit_t[nt] = price_t[nt] * c_t[n_total] - cost_t[nt] * effort_t[nt] ^ beta;

  ppue_hat_t[nt] = profit_t[nt] / effort_t[nt];

  p_age_sampled = cn_ta[n_total, 1:n_ages] / sum(cn_ta[n_total, 1:n_ages]);

  p_lbin_sampled[nt, 1:n_lbins] = (p_age_sampled * length_at_age_key) / sum(p_age_sampled * length_at_age_key);

    for (t in 1:(nt - 1)){

    delta_f[t] = effort_t[t + 1] - effort_t[t];

    perc_change_effort_hat[t] = effort_t[t + 1] / effort_t[t];


    }

  ppue_hat = delta_f / p_response;

}

model{

real new_f;

//// length comps likelihood ////

for (i in 1:(n_lcomps)){

  length_comps[i, 1:n_lbins] ~ multinomial(to_vector(p_lbin_sampled[length_comps_years[i], 1:n_lbins]));

} // close length likelihood

//// effort prior ////

if (economic_model == 0){

/*
apply simple penalty on year to year variation in F
*/

    for (t in 2:nt){

    f_t[t] ~ normal(f_t[t - 1] , sigma_f);

    } // close time loop


} // close effort 0


if (economic_model == 1) { // open access priors
/*
Set priors on the change in f based on open acces dynamics and
relative changes in price, cost, and q where available
*/

  for (t in 2:nt){

    new_f = (effort_t[t - 1] + (p_response * (ppue_hat_t[t - 1]))) * q_t[t];

    f_t[t] ~ normal(new_f,sigma_f);

    } // close time loop

} // close 1

if (economic_model == 2){

/*
Set priors on the change in F base on knowledge on the relative change in effort
*/

    for (t in 2:nt){

    new_f = q_t[t] * effort_t[t - 1] * perc_change_effort[t - 1];

    f_t[t] ~ normal(new_f , sigma_f);

    } // close time loop

} // close 2



if (economic_model == 3){

/*
Fit to ppue data, under the assumption that ppue is proportional to delta F
*/

  ppue_t[1:(nt - 1)] ~ normal(ppue_hat, sigma_obs);

    // for (t in 2:nt){
    //
    // f_t[t] ~ normal(f_t[t - 1] , sigma_f);
    //
    // } // close time loop

}


if (economic_model == 4){

  /*
Fit to ppue data, under the assumption that ppue is proportional to delta F,
and provide priors using economic knowledge as well
*/

  for (t in 2:nt){

    new_f = (effort_t[t - 1] + (p_response * (ppue_hat_t[t - 1]))) * q_t[t];

    f_t[t] ~ normal(new_f,sigma_f);

    } // close time loop

  ppue_t[1:(nt - 1)] ~ normal(ppue_hat, sigma_obs);

}


if (economic_model == 5){
  /*
fit to percentage change in effort
*/
  perc_change_effort[1:(nt - 1)] ~ normal(perc_change_effort_hat, sigma_f);

}

if (economic_model == 6){
  /*
fit to percentage change in effort with bioeconomic informed priors
*/

    for (t in 2:nt){

    new_f = (effort_t[t - 1] + (p_response * (ppue_hat_t[t - 1]))) * q_t[t];

    f_t[t] ~ normal(new_f,sigma_f);

    } // close time loop

  perc_change_effort[1:(nt - 1)] ~ normal(perc_change_effort_hat, sigma_obs);

}


log_max_cost ~ normal(log(max_cost_guess),1);

f_t ~ cauchy(0,2.5);

sigma_obs ~ cauchy(0, 2.5);

log_p_response ~ normal(log(.1),2);

//// recruitment prior ////

uc_rec_dev_t ~ normal(0, sigma_r);

sigma_r ~ normal(sigma_r_guess, sd_sigma_r);

//// selectivity likelihood ////

p_length_50_sel ~ normal(length_50_sel_guess/loo, .05);

// sel_delta ~ normal(delta_guess, 2);


} // close model block

generated quantities{

int n_tl[nt, n_lbins];

vector[nt - 1] pp_ppue;

vector[nt - 1] pp_perc_change_effort;

for (t in 1:nt){

    n_tl[t, 1:n_lbins] = multinomial_rng(to_vector(p_lbin_sampled[t, 1:n_lbins]),500); // generate length comp samples


    if (t < nt){

      pp_ppue[t] = normal_rng(ppue_hat[t], sigma_obs);

      pp_perc_change_effort[t] = normal_rng(perc_change_effort_hat[t], sigma_obs);
    }
}

} // close generated quantities block
