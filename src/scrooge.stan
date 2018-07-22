/*scrooge_v3.0

effort prior is open-access dynamics

non-centered recruitment and effort deviations

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

row_vector[nt] cost_t;

row_vector[nt] q_t;

vector[nt] ppue_t;

real beta;

real length_50_sel_guess;

real delta_guess;

real<lower = 0> p_response_guess;

real<lower = 0> cv_effort;

//// biology ////


vector<lower=0>[n_ages] mean_length_at_age;

vector<lower=0>[n_ages] mean_weight_at_age;

vector<lower=0, upper=1>[n_ages] mean_maturity_at_age;

matrix[n_ages, n_lbins] length_at_age_key;

real m; //natural mortality

real h; //steepness

real r0; // virgin recruitment

real k; //vbk growth

real loo; // vbk linf

real t0; // t0

real sigma_r_guess;

real sd_sigma_r;

real sigma_effort_guess;

}

transformed data{



}

parameters{

// real log_burn_effort; // burn in effort

// vector[nt]  log_effort_t; // effort in time t

real<lower = 0> burn_f;

vector<lower = 0>[nt]  f_t; // f in time t

vector[nt]  uc_rec_dev_t; //  recruitment deviates

real <lower = 0> sigma_r; // standard deviation of recruitment deviates

// real<lower = 0> sigma_effort;

real<lower = 0> sigma_ppue;

real  log_p_response;

real<lower = 0, upper = .9> p_length_50_sel; // length at 50% selectivity

// real<lower = 0, upper = 10> p_response;

} // close parameters block

transformed parameters{

  real ssb0;

  real ssb_temp;

  real sel_delta;

  real length_50_sel;

  real plus_group;

  row_vector[n_ages] temp_a;

  row_vector[n_ages - 1] temp_a2;

  vector[nt]  effort_t; //  base effort multiplier

  vector[nt]  rec_dev_t; //  base effort multiplier

  matrix[n_burn, n_ages] n_a_init; // numbers at time and age

  matrix[n_burn, n_ages] ssb_init; // numbers at time and age

  matrix[nt, n_ages] n_ta; // numbers at time and age

  matrix[nt, n_ages] ssb_ta; // spawning stock biomass at time and age

  matrix[nt, n_ages] c_ta; // catch (biomass) at time and age

  matrix[nt, n_ages] cn_ta; // catch (numbers) at time and age

  matrix[nt, n_lbins] p_lbin_sampled; // numbers at time and length bin

  vector[nt] c_t; // total catch at time step

  vector[nt] profit_t; // profits

  vector[nt] ppue_hat_t; // profit per unit effort

  vector[n_ages] mean_selectivity_at_age; // mean selectivity at age

  row_vector[n_ages] p_age_sampled;

  vector[nt - 1] delta_f;

  vector[nt - 1] ppue_hat;

  real sigma_effort;

  sigma_effort = sigma_effort_guess;

  // length_50_sel = length_50_sel_guess;

  length_50_sel = loo * p_length_50_sel;

  sel_delta = 2;

  rec_dev_t = exp(uc_rec_dev_t);

  // fill matrices with zeros //

  n_ta = rep_matrix(0,nt, n_ages);

  ssb_ta = rep_matrix(0,nt, n_ages);

  c_ta = rep_matrix(0,nt, n_ages);

  cn_ta = rep_matrix(0,nt, n_ages);

  p_lbin_sampled = rep_matrix(0,nt, n_lbins);

  mean_selectivity_at_age = 1.0 ./ (1 + exp(-log(19) * ((mean_length_at_age - length_50_sel) / sel_delta))); // selectivity ogive at age

// set up initial population //

  // f_t[1] = effort_t[1] .* q_t[1];

  effort_t[1] = f_t[1] / q_t[1];

  n_a_init[1,1:n_ages] = r0 * exp(-m * (ages - 1))';

  plus_group = n_a_init[1, n_ages - 1] * exp(-m) / (1 - exp(-m));  //plus group

  n_a_init[1, n_ages] = plus_group;

  temp_a = n_a_init[1, 1:n_ages] .* mean_maturity_at_age'.* mean_weight_at_age';

  ssb_init[1, 1:n_ages] = temp_a; // calculate ssb at age

  ssb0 = sum(ssb_init[1, 1:n_ages]);

 for (t in 2:n_burn){

    ssb_temp =  sum(ssb_init[t - 1, 1:n_ages]);

    n_a_init[t,1] = (0.8 * r0 * h *ssb_temp) / (0.2 * ssb0 * (1 - h) + (h - 0.2) * ssb_temp); //calculate recruitment

  temp_a2 = n_a_init[t - 1, 1:(n_ages -1)] .* (exp(-(m + burn_f* mean_selectivity_at_age[1:(n_ages - 1)])))';

    n_a_init[t , 2:n_ages] = temp_a2;  // grow and die

    n_a_init[t, n_ages] = n_a_init[t, n_ages] + n_a_init[t - 1, n_ages] .* (exp(-(m + burn_f * mean_selectivity_at_age[n_ages])))'; // assign to plus group

  ssb_init[t, 1:n_ages] = n_a_init[t, 1:n_ages] .* mean_maturity_at_age'.* mean_weight_at_age'; // calculate ssb at age

}

  // Start observed period

  n_ta[1, 1:n_ages] = n_a_init[n_burn, 1:n_ages];  // start at initial depletion

  n_ta[1,1] = n_ta[1,1] * rec_dev_t[1]; // initial recruitment deviate

  ssb_ta[1, 1:n_ages] = n_ta[1, 1:n_ages] .* mean_maturity_at_age' .* mean_weight_at_age'; // virgin ssb


  // order of events spawn, grow and die, recruit

  for (t in 2:nt){

    ssb_temp =  sum(ssb_ta[t - 1, 1:n_ages]);

    n_ta[t,1] = ((0.8 * r0 * h *ssb_temp) / (0.2 * ssb0 * (1 - h) + (h - 0.2) * ssb_temp)) * rec_dev_t[t]; //calculate recruitment

    temp_a2 = n_ta[t - 1, 1:(n_ages -1)] .* (exp(-(m + f_t[t - 1] * mean_selectivity_at_age[1:(n_ages - 1)])))';

    n_ta[t , 2:n_ages] = temp_a2; // grow and die

    n_ta[t, n_ages] = n_ta[t, n_ages] + n_ta[t - 1, n_ages] .* (exp(-(m + f_t[t - 1] * mean_selectivity_at_age[n_ages])))'; // assign to plus group

  ssb_ta[t, 1:n_ages] = n_ta[t, 1:n_ages] .* mean_maturity_at_age'.* mean_weight_at_age'; // calculate ssb at age

   cn_ta[t - 1, 1:n_ages] = ((f_t[t - 1] * mean_selectivity_at_age) ./ (m + f_t[t - 1] * mean_selectivity_at_age))' .* n_ta[t - 1, 1:n_ages] .* (1 - exp(-(m + f_t[t - 1] * mean_selectivity_at_age)))'; // .* mean_weight_at_age';

   c_ta[t-1, 1:n_ages] = cn_ta[t - 1, 1:n_ages] .* mean_weight_at_age';

   c_t[t - 1] = sum(c_ta[t-1, 1:n_ages]);

  // run economic model

  profit_t[t - 1] = price_t[t - 1] * c_t[t - 1] - cost_t[t - 1] *         effort_t[t - 1] ^ beta;

  ppue_hat_t[t - 1] = profit_t[t - 1] / effort_t[t - 1];

  effort_t[t] = f_t[t] / q_t[t];

  // sample lengths //

  p_age_sampled = cn_ta[t - 1, 1:n_ages] / sum(cn_ta[t - 1, 1:n_ages]);

  p_lbin_sampled[t - 1, 1:n_lbins] = (p_age_sampled * length_at_age_key) / sum(p_age_sampled * length_at_age_key);  // calculate the proportion of the sampled catch in each length bin


  } // close time loop


// fill in final time step
     cn_ta[nt, 1:n_ages] = ((f_t[nt] * mean_selectivity_at_age) ./ (m + f_t[nt] * mean_selectivity_at_age))' .* n_ta[nt, 1:n_ages] .* (1 - exp(-(m + f_t[nt] * mean_selectivity_at_age)))'; //

   c_ta[nt, 1:n_ages] = cn_ta[nt, 1:n_ages] .* mean_weight_at_age';

   c_t[nt] = sum(c_ta[nt, 1:n_ages]);

  profit_t[nt] = price_t[nt] * c_t[nt] - cost_t[nt] * effort_t[nt] ^ beta;

  ppue_hat_t[nt] = profit_t[nt] / effort_t[nt];

  p_age_sampled = cn_ta[nt, 1:n_ages] / sum(cn_ta[nt, 1:n_ages]);

  p_lbin_sampled[nt, 1:n_lbins] = (p_age_sampled * length_at_age_key) / sum(p_age_sampled * length_at_age_key);

    for (t in 1:(nt - 1)){

    delta_f[t] = effort_t[t + 1] - effort_t[t];

    }

  ppue_hat = delta_f / exp(log_p_response);


}

model{

real oa_prediction;

real effort_data_prediction;

real new_f;

real previous_max;


//// length comps likelihood ////

for (i in 1:(n_lcomps)){

  length_comps[i, 1:n_lbins] ~ multinomial(to_vector(p_lbin_sampled[length_comps_years[i], 1:n_lbins]));

} // close length likelihood

//// effort prior ////


if (economic_model == 1) { // open access priors

  for (t in 2:nt){

    new_f = (effort_t[t - 1] + (p_response_guess * (ppue_hat_t[t - 1]))) * q_t[t];

    f_t[t] ~ normal(new_f,1e-6);

    } // close time loop

} // close 1

if (economic_model == 2){

    for (t in 2:nt){

    new_f = q_t[t] * effort_t[t - 1] * perc_change_effort[t - 1];

    f_t[t] ~ normal(new_f , sigma_effort);

    } // close time loop

} // close 2

if (economic_model == 0){

  for (i in 2:nt){

  f_t[i] ~ normal(f_t[i - 1],sigma_effort);

}

} // close effort 0

if (economic_model == 3){

  // for (i in 2:nt){
  //
  //   f_t[i] ~ normal(f_t[i - 1],sigma_effort);
  //
  // }

  ppue_t[1:(nt - 1)] ~ normal(ppue_hat, sigma_ppue);

    f_t[nt] ~ normal(f_t[nt - 1], .001);


}

// if (economic_model == 3){
//
//   // for (i in 2:nt){
//   //
//   //   f_t[i] ~ normal(f_t[i - 1],sigma_effort);
//   //
//   // }
//
//   ppue_t[1:(nt - 1)] ~ normal(ppue_hat, sigma_ppue);
//
//   f_t[nt] ~ normal(f_t[nt - 1], .001);
//
// }


sigma_ppue ~ cauchy(0, 2.5);

// sigma_effort ~ normal(sigma_effort_guess,.1);

log_p_response ~ normal(log(.1),2);

// p_response ~ normal(p_response_guess,.001); // constrain p_response

//// recruitment prior ////

uc_rec_dev_t ~ normal(0, sigma_r);

sigma_r ~ normal(sigma_r_guess, sd_sigma_r);

//// selectivity likelihood ////

// p_length_50_sel ~ normal(length_50_sel_guess/loo, 5);

// sel_delta ~ normal(delta_guess, 2);


} // close model block

generated quantities{


int n_tl[nt, n_lbins];

for (t in 1:nt){

    n_tl[t, 1:n_lbins] = multinomial_rng(to_vector(p_lbin_sampled[t, 1:n_lbins]),500); // generate length comp samples
}

} // close generated quantities block
