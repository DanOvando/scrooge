/*scrooge_v3.0

effort prior is open-access dynamics

non-centered recruitment and effort deviations

*/


data{

//// run parameters ////

int<lower = 1> nt; // number of time steps to run

int n_ages; // number of ages

int n_lbins; // number of length bins;

int n_lcomps; // number of timesteps of length comps available

vector[n_ages] ages; // vector of ages

int<lower = 0> economic_model; // 0 = no bioeconomic prior, 1 = bioeconomic prior,

int<lower = 0, upper = 1> use_effort_data; // 0 = don't use effort data, 1 = use effort data,



//// length data ////

int length_comps[n_lcomps,n_lbins];

int length_comps_years[n_lcomps]; // time steps in which length comps are available

//// economic data ////

row_vector[nt] relative_effort;

row_vector[nt] price_t;

row_vector[nt] cost_t;

row_vector[nt] q_t;

real beta;

real length_50_sel_guess;

real delta_guess;

real<lower = 0> p_response_guess;

real<lower = 0> p_msy;

real<lower = 0> e_msy;

real<lower = 0> f_init_guess;

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



}

parameters{

// real <lower = -5, upper = 3> log_base_effort;

vector[nt - economic_model]  effort_shock_t; //  base effort multiplier

real<lower = -5, upper = 5> log_sigma_effort; // standard deviation of effort

vector[nt]  uc_rec_dev_t; //  recruitment deviates

real<lower = -5, upper = 2> log_sigma_r; // standard deviation of recruitment deviates

real<lower = 0, upper = 1.5> p_length_50_sel; // length at 50% selectivity

real<lower = 1e-3, upper = 2> f_init; // f to get to initial depletion

real<lower = 0, upper = 1> p_response;

// real<lower = 0.001> sel_delta; // difference between length 50% selected and length 95% selected

} // close parameters block

transformed parameters{

  real ssb0;

  real ssb_temp;

  real sigma_r;

  real sigma_effort;

  real new_effort;

  real sel_delta;

  real length_50_sel;

  real base_effort;

  vector[nt]  rec_dev_t; //  base effort multiplier

  matrix[nt, n_ages] n_ta; // numbers at time and age

  matrix[nt, n_ages] ssb_ta; // spawning stock biomass at time and age

  matrix[nt, n_ages] c_ta; // catch (biomass) at time and age

  matrix[nt, n_ages] cn_ta; // catch (numbers) at time and age

  matrix[nt, n_lbins] p_lbin_sampled; // numbers at time and length bin

  vector[nt] c_t; // total catch at time step

  vector[nt] f_t; // fishing mortality at time step

  vector[nt] total_effort_t; // effort in right units

  vector[nt] profit_t; // profits

  vector[n_ages] mean_selectivity_at_age; // mean selectivity at age

  // vector[nt] profit_per_unit_effort;

  row_vector[n_ages] n_a_0;

  row_vector[n_ages] p_age_sampled;

  length_50_sel = loo * p_length_50_sel;

  sel_delta = 1;

  sigma_r = exp(log_sigma_r);

  sigma_effort = exp(log_sigma_effort);

  rec_dev_t = exp(sigma_r * uc_rec_dev_t - sigma_r^2/2);

  // base_effort = exp(log_base_effort);

  base_effort = f_init / q_t[1];

  // fill matrices with zeros //

  n_ta = rep_matrix(0,nt, n_ages);

  // b_ta = rep_matrix(0,nt, n_ages);

  ssb_ta = rep_matrix(0,nt, n_ages);

  c_ta = rep_matrix(0,nt, n_ages);

  cn_ta = rep_matrix(0,nt, n_ages);

  p_lbin_sampled = rep_matrix(0,nt, n_lbins);

  // profit_shock = rep_matrix(0, nt, 2);

  mean_selectivity_at_age = 1.0 ./ (1 + exp(-log(19) * ((mean_length_at_age - length_50_sel) / sel_delta))); // selectivity ogive at age

  ssb0 = sum((r0 * exp(-m * (ages - 1)))  .* mean_maturity_at_age .* mean_weight_at_age); // virgin ssb

  n_a_0[1] = r0;

  for (i in 2:n_ages){

    n_a_0[i] = n_a_0[i-1] * exp(-(m + f_init .* mean_selectivity_at_age[i])); // transpose to specify row format

  }


  n_ta[1, 1:n_ages] = n_a_0 * rec_dev_t[1];  // add in virign recruitment to get population moving

  n_ta[1, n_ages] = n_ta[1, n_ages - 1] * exp(-(m + f_init)) / (1 - exp(-(m + f_init)));  //plus group

  // b_ta[1, 1:n_ages] = n_ta[1, 1:n_ages] .* mean_weight_at_age'; // virgin biomass

  ssb_ta[1, 1:n_ages] = n_ta[1, 1:n_ages] .* mean_maturity_at_age' .* mean_weight_at_age'; // virgin ssb


  if (economic_model == 1) {
  total_effort_t[1] = base_effort;
  } else {

   total_effort_t[1] = base_effort + sigma_effort * effort_shock_t[1];

  }

  f_t[1] = total_effort_t[1] .* q_t[1];

  // order of events spawn, grow and die, recruit

  for (t in 2:nt){

    ssb_temp =  sum(ssb_ta[t - 1, 1:n_ages]);

    n_ta[t,1] = ((0.8 * r0 * h *ssb_temp) / (0.2 * ssb0 * (1 - h) + (h - 0.2) * ssb_temp)) * rec_dev_t[t]; //calculate recruitment

    n_ta[t , 2:n_ages] = n_ta[t - 1, 1:(n_ages -1)] .* (exp(-(m + f_t[t - 1] * mean_selectivity_at_age[1:(n_ages - 1)])))'; // grow and die

    n_ta[t, n_ages] = n_ta[t, n_ages] + n_ta[t - 1, n_ages] .* (exp(-(m + f_t[t - 1] * mean_selectivity_at_age[n_ages])))'; // assign to plus group

  // b_ta[t, 1:n_ages] = n_ta[t, 1:n_ages] .* mean_weight_at_age'; // calculate biomass at age

  ssb_ta[t, 1:n_ages] = n_ta[t, 1:n_ages] .* mean_maturity_at_age'.* mean_weight_at_age'; // calculate ssb at age

   cn_ta[t - 1, 1:n_ages] = ((f_t[t - 1] * mean_selectivity_at_age) ./ (m + f_t[t - 1] * mean_selectivity_at_age))' .* n_ta[t - 1, 1:n_ages] .* (1 - exp(-(m + f_t[t - 1] * mean_selectivity_at_age)))'; // .* mean_weight_at_age';

   c_ta[t-1, 1:n_ages] = cn_ta[t - 1, 1:n_ages] .* mean_weight_at_age';

   c_t[t - 1] = sum(c_ta[t-1, 1:n_ages]);

  // print(c_t[t - 1])

  // run economic model

  profit_t[t - 1] = price_t[t - 1] * c_t[t - 1] - cost_t[t - 1] * total_effort_t[t - 1] ^ beta;

if (economic_model == 1){

new_effort = total_effort_t[t - 1] + e_msy * (p_response * (profit_t[t - 1] / p_msy)) + sigma_effort * effort_shock_t[t - 1];

} else {

new_effort = base_effort + sigma_effort * effort_shock_t[t];

}


if (new_effort <= 0){

  new_effort = -.01 / (new_effort - 1); // if effort is negative, have it asymptotically approach zero instead

}

 total_effort_t[t] = new_effort; //* exp(sigma_effort * effort_shock_t[t - 1]);

  f_t[t] = total_effort_t[t] * q_t[t];

  // sample lengths //

  p_age_sampled = cn_ta[t - 1, 1:n_ages] / sum(cn_ta[t - 1, 1:n_ages]);

  p_lbin_sampled[t - 1, 1:n_lbins] = (p_age_sampled * length_at_age_key) / sum(p_age_sampled * length_at_age_key);  // calculate the proportion of the sampled catch in each length bin


  } // close time loop


// fill in final time step
     cn_ta[nt, 1:n_ages] = ((f_t[nt] * mean_selectivity_at_age) ./ (m + f_t[nt] * mean_selectivity_at_age))' .* n_ta[nt, 1:n_ages] .* (1 - exp(-(m + f_t[nt] * mean_selectivity_at_age)))'; //

   c_ta[nt, 1:n_ages] = cn_ta[nt, 1:n_ages] .* mean_weight_at_age';

   c_t[nt] = sum(c_ta[nt, 1:n_ages]);

  profit_t[nt] = price_t[nt] * c_t[nt] - cost_t[nt] * total_effort_t[nt] ^ beta;

  p_age_sampled = cn_ta[nt, 1:n_ages] / sum(cn_ta[nt, 1:n_ages]);

  p_lbin_sampled[nt, 1:n_lbins] = p_age_sampled * length_at_age_key / sum(p_age_sampled * length_at_age_key);

}

model{

vector[nt] relative_total_effort_t; // effort in right units

matrix[n_lcomps, n_lbins] pn_tl; // numbers at time and length bin

//// length comps likelihood ////

pn_tl = rep_matrix(0,n_lcomps, n_lbins);

for (i in 1:(n_lcomps - 1)){

  length_comps[i, 1:n_lbins] ~ multinomial(to_vector(p_lbin_sampled[length_comps_years[i], 1:n_lbins]));

} // close length likelihood

//// effort prior ////

if (use_effort_data == 1){

relative_total_effort_t = total_effort_t / max(total_effort_t); // relative effort

relative_effort ~ normal(relative_total_effort_t, 0.1);
}

f_init ~ normal(f_init_guess, 1);


if (economic_model == 1) {

effort_shock_t ~ normal(0,1);

} else{

  for (i in 2:nt){

  effort_shock_t[i] ~ normal(effort_shock_t[i - 1],1);

}


} // close effort loop

log_sigma_effort ~ normal(0,0.25);

p_response ~ normal(p_response_guess,.01); // constrain p_response

//// recruitment prior ////

uc_rec_dev_t ~ normal(0, 1);

log_sigma_r ~ normal(log(0.2), 0.01);

//// selectivity likelihood ////

p_length_50_sel ~ normal(length_50_sel_guess/loo, .01);

// sel_delta ~ normal(delta_guess, 2);


} // close model block

generated quantities{


int n_tl[nt, n_lbins];

for (t in 1:nt){
    n_tl[t, 1:n_lbins] = multinomial_rng(to_vector(p_lbin_sampled[t, 1:n_lbins]),500); // generate length comp samples
}

} // close generated quantities block

