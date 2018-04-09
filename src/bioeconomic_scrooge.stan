/*scrooge_v3.0

effort prior is open-access dynamics

non-centered recruitment and effort deviations

*/


data{

//// run parameters ////

int<lower = 1> nt; // number of time steps to run

int n_ages; // number of ages

int n_lbins; // number of length bins;

int n_lcomps; // number of timesteps of length comps availabl

vector[n_ages] ages; // vector of ages

int estimate_recruits; // 1 if you want to estimate recruitment deviates, 0 if not

int<lower = 0> economic_model; // 0 = no economic prior, 1 = economic shock prior,

//// length data ////

int length_comps[n_lcomps,n_lbins];

int length_comps_years[n_lcomps]; // time steps in which length comps are available

//// economic data ////

// vector<lower=0, upper=1>[n_ages] mean_selectivity_at_age;

matrix[nt,2] price_t;

matrix[nt,2] cost_t;

matrix[nt,2] q_t;

real beta;

real length_50_sel_guess;

real delta_guess;

// real max_change;

// real<lower  = 0> base_effort;
// vector<lower = 0>[nt] f_t;

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

}

transformed data{

  // int n_sampled;

  // n_sampled = 1;


}

parameters{

real <lower = 0.001, upper = 10> base_effort;

// vector<lower = -5, upper = 5>[nt - 1]  effort_shock_t; //  base effort multiplier

vector[nt - 1]  effort_shock_t; //  base effort multiplier

real<lower = -5, upper = 5> log_sigma_effort; // standard deviation of effort

vector[nt]  uc_rec_dev_t; //  recruitment deviates

real<lower = -5, upper = 2> log_sigma_r; // standard deviation of recruitment deviates

real<lower = 0, upper = 1.5> p_length_50_sel; // length at 50% selectivity

// real<lower = 0.001> sel_delta; // difference between length 50% selected and length 95% selected

// real<lower = 0, upper = 1> p_response;
// real<lower = 1, upper = 4*(m/.001)> base_effort;


} // close parameters block

transformed parameters{

  real ssb0;

  real ssb_temp;

  real sigma_r;

  real sigma_effort;

  real mean_price;

  real new_effort;

  real sel_delta;

  real length_50_sel;

  // real p_response;

  // vector[nt]  effort_t; //  base effort multiplier

  vector[nt]  rec_dev_t; //  base effort multiplier

  matrix[nt, n_ages] n_ta; // numbers at time and age

  matrix[nt, n_ages] b_ta; // biomass at time and age

  matrix[nt, n_ages] ssb_ta; // spawning stock biomass at time and age

  matrix[nt, n_ages] c_ta; // catch (biomass) at time and age

  matrix[nt, n_ages] cn_ta; // catch (numbers) at time and age

  matrix[nt, n_lbins] p_lbin_sampled; // numbers at time and length bin

  vector[nt] c_t; // total catch at time step

  vector[nt] f_t; // fishing mortality at time step

  vector[nt] total_effort_t; // effort in right units

  vector[nt] profit_t; // profits

  vector[nt] sq_profit_t; // profits under last years price and cost

  vector[nt] profit_shock_t; // profit shock

  vector[n_ages] mean_selectivity_at_age; // mean selectivity at age

  vector[nt] profit_per_unit_effort;

  // vector[nt] delta_f; // expected change in f from profit shock

  row_vector[n_ages] temp_n_a;

  row_vector[n_ages] p_age_sampled;

  // real approx_rmsy;
  //
  // row_vector[n_ages] approx_nmsy;
  //
  // real approx_msy;
  //
  // real approx_pmsy;
  // transform parameters

  length_50_sel = loo * p_length_50_sel;

  sel_delta = 1;

  sigma_r = exp(log_sigma_r);

  sigma_effort = exp(log_sigma_effort);

  // p_response = exp(log_p_response);
  // total_effort_t = base_effort;

  rec_dev_t = exp(sigma_r * uc_rec_dev_t - sigma_r^2/2);

  // fill matrices with zeros
  n_ta = rep_matrix(0,nt, n_ages);

  b_ta = rep_matrix(0,nt, n_ages);

  ssb_ta = rep_matrix(0,nt, n_ages);

  c_ta = rep_matrix(0,nt, n_ages);

  cn_ta = rep_matrix(0,nt, n_ages);

  p_lbin_sampled = rep_matrix(0,nt, n_lbins);

  // profit_shock = rep_matrix(0, nt, 2);
  // print("wtf")
  mean_selectivity_at_age = 1.0 ./ (1 + exp(-log(19) * ((mean_length_at_age - length_50_sel) / sel_delta))); // selectivity ogive at age

  // total_effort_t = effort_t * base_effort; // convert effort into correct scale


  // ssb0 = sum((r0 * exp(-m * (ages - 1)))  .* mean_maturity_at_age .* mean_weight_at_age); // virgin ssb

  // approx_rmsy = ((0.8 * r0 * h * 0.4 *ssb_temp) / (0.2 * ssb0 * (1 - h) + (h - 0.2) * ssb_temp * 0.4));
  //
  // approx_nmsy = approx_rmsy * exp(-2*m * (ages' - 1));
  //
  // approx_msy = sum(((m * mean_selectivity_at_age) ./ (m + m * mean_selectivity_at_age))' .* approx_nmsy .* (1 - exp(-(m + m * mean_selectivity_at_age)))' .* mean_weight_at_age'); // .* mean_weight_at_age';
  //
  // approx_pmsy = mean(price_t[1:nt,1]) * approx_msy - mean(cost_t[1:nt,1]) * (m / mean(q_t[1:nt,1])) ^ beta;


 ssb0 = sum((r0 * rec_dev_t[1] * exp(-m * (ages - 1)))  .* mean_maturity_at_age .* mean_weight_at_age); // virgin ssb
  //
    // ssb0 = sum((r0 * exp(-m * (ages - 1)))  .* mean_maturity_at_age .* mean_weight_at_age); // virgin ssb

  mean_price = mean(price_t[1:nt,1]);

  temp_n_a = (r0 *  rec_dev_t[1]) * exp(-m * (ages' - 1)); // transpose to specify row format

  // temp_n_a = (r0) * exp(-m * (ages' - 1)); // transpose to specify row format

  n_ta[1, 1:n_ages] = temp_n_a;  // add in virign recruitment to get population moving

  n_ta[1, n_ages] = (n_ta[1, n_ages - 1] * exp(-m)) / (1 - exp(-m));  //plus group

  b_ta[1, 1:n_ages] = n_ta[1, 1:n_ages] .* mean_weight_at_age'; // virgin biomass

  ssb_ta[1, 1:n_ages] = b_ta[1, 1:n_ages] .* mean_maturity_at_age'; // virgin ssb

  total_effort_t[1] = base_effort;

  f_t[1] = total_effort_t[1] .* q_t[1,1];

  // order of events spawn, grow and die, recruit

  for (t in 2:nt){

    ssb_temp =  sum(ssb_ta[t - 1, 1:n_ages]);

    n_ta[t,1] = ((0.8 * r0 * h *ssb_temp) / (0.2 * ssb0 * (1 - h) + (h - 0.2) * ssb_temp)) * rec_dev_t[t]; //calculate recruitment

    n_ta[t , 2:n_ages] = n_ta[t - 1, 1:(n_ages -1)] .* (exp(-(m + f_t[t - 1] * mean_selectivity_at_age[1:(n_ages - 1)])))'; // grow and die

    n_ta[t, n_ages] = n_ta[t - 1, n_ages] .* (exp(-(m + f_t[t - 1] * mean_selectivity_at_age[n_ages])))'; // assign to plus group

  b_ta[t, 1:n_ages] = n_ta[t, 1:n_ages] .* mean_weight_at_age'; // calculate biomass at age

  ssb_ta[t, 1:n_ages] = b_ta[t, 1:n_ages] .* mean_maturity_at_age'; // calculate ssb at age

   cn_ta[t - 1, 1:n_ages] = ((f_t[t - 1] * mean_selectivity_at_age) ./ (m + f_t[t - 1] * mean_selectivity_at_age))' .* n_ta[t - 1, 1:n_ages] .* (1 - exp(-(m + f_t[t - 1] * mean_selectivity_at_age)))'; // .* mean_weight_at_age';

   c_ta[t-1, 1:n_ages] = cn_ta[t - 1, 1:n_ages] .* mean_weight_at_age';

   c_t[t - 1] = sum(c_ta[t-1, 1:n_ages]);

  // print(c_t[t - 1])

  // run economic model

  profit_t[t - 1] = price_t[t - 1,1] * c_t[t - 1] - cost_t[t - 1,1] * total_effort_t[t - 1] ^ beta;

  sq_profit_t[t - 1] = price_t[t - 1,2] * c_t[t - 1] - cost_t[t - 1,2] * total_effort_t[t - 1] ^ beta;

  profit_shock_t[t - 1] = profit_t[t - 1] - sq_profit_t[t - 1];

// print(exp(p_response * (profit_t[t - 1] / total_effort_t[t - 1])))


// effort_expansion[t - 1]  = fabs(exp(p_response * (profit_t[t - 1] / total_effort_t[t - 1])) - 1);

profit_per_unit_effort[t - 1]  = profit_t[t - 1] / total_effort_t[t - 1];

new_effort = (total_effort_t[t - 1] + 0.01 * (profit_t[t - 1] / (ssb0 * mean_price))) + sigma_effort * effort_shock_t[t - 1];

if (new_effort <= 0){

  new_effort = -.01 / (new_effort - 1); // if effort is negative, have it asymptotically approach zero instead

}

 total_effort_t[t] = new_effort; //* exp(sigma_effort * effort_shock_t[t - 1]);

 // total_effort_t[t] = (total_effort_t[t - 1] * exp(.001 * (profit_t[t - 1] / (ssb0 * mean_price)))) * exp(sigma_effort * effort_shock_t[t - 1]);

 // print(profit_t[t - 1])


  f_t[t] = total_effort_t[t] * q_t[t,1];

  // sample lengths

  p_age_sampled = cn_ta[t - 1, 1:n_ages] / sum(cn_ta[t - 1, 1:n_ages]);

  p_lbin_sampled[t - 1, 1:n_lbins] = (p_age_sampled * length_at_age_key) / sum(p_age_sampled * length_at_age_key);  // calculate the proportion of the sampled catch in each length bin


  } // close time loop


// fill in final time step
     cn_ta[nt, 1:n_ages] = ((f_t[nt] * mean_selectivity_at_age) ./ (m + f_t[nt] * mean_selectivity_at_age))' .* n_ta[nt, 1:n_ages] .* (1 - exp(-(m + f_t[nt] * mean_selectivity_at_age)))'; //

   c_ta[nt, 1:n_ages] = cn_ta[nt, 1:n_ages] .* mean_weight_at_age';

   c_t[nt] = sum(c_ta[nt, 1:n_ages]);

  profit_t[nt] = price_t[nt,1] * c_t[nt] - cost_t[nt,1] * total_effort_t[nt] ^ beta;

  sq_profit_t[nt] = price_t[nt,2] * c_t[nt] - cost_t[nt,2] * total_effort_t[nt] ^ beta;

  profit_shock_t[nt] = profit_t[nt] - sq_profit_t[nt];

  p_age_sampled = cn_ta[nt, 1:n_ages] / sum(cn_ta[nt, 1:n_ages]);

  p_lbin_sampled[nt, 1:n_lbins] = p_age_sampled * length_at_age_key / sum(p_age_sampled * length_at_age_key);


  profit_per_unit_effort[nt]  = profit_t[nt] / total_effort_t[nt];

  // effort_expansion[nt]  = fabs(exp(p_response * (profit_t[nt] / total_effort_t[nt])) - 1);

}

model{

matrix[n_lcomps, n_lbins] pn_tl; // numbers at time and length bin

vector[nt] economic_prior;

int temp_n_tl[n_lbins];

//// length comps likelihood ////

pn_tl = rep_matrix(0,n_lcomps, n_lbins);

for (i in 1:(n_lcomps - 1)){

  length_comps[i, 1:n_lbins] ~ multinomial(to_vector(p_lbin_sampled[length_comps_years[i], 1:n_lbins]));

} // close length likelihood

//// effort prior ////

effort_shock_t ~ normal(0,1);

log_sigma_effort ~ normal(0,0.1);

// base_effort ~ normal(log(m / .001), 1);

// p_response ~ normal(0,.1); // constrain p_response

// p_response ~ normal(log(mean_change) / mean(profit_per_unit_effort), .1);

//// recruitment prior ////

uc_rec_dev_t ~ normal(0, 1);

// sigma_r ~ normal(0.7, 0.2);

log_sigma_r ~ normal(log(sigma_r_guess), 0.1);

//// selectivity likelihood ////

p_length_50_sel ~ normal(length_50_sel_guess/loo, .1 * loo);

// sel_delta ~ normal(delta_guess, 2);


} // close model block

generated quantities{


int n_tl[nt, n_lbins];

for (t in 1:nt){
    n_tl[t, 1:n_lbins] = multinomial_rng(to_vector(p_lbin_sampled[t, 1:n_lbins]),500); // generate length comp samples
}

} // close generated quantities block

