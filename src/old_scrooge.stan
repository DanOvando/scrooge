data{

//// run parameters ////

int<lower = 1> nt; // number of time steps to run

int n_ages; // number of ages

int n_lbins; // number of length bins;

int n_lcomps; // number of timesteps of length comps availabl

vector[n_ages] ages; // vector of ages

int estimate_recruits; // 1 if you want to estimate recruitment deviates, 0 if not

//// length data ////

int length_comps[n_lcomps,n_lbins];

int length_comps_years[n_lcomps]; // time steps in which length comps are available

//// economic data ////

vector<lower=0, upper=1>[n_ages] mean_selectivity_at_age;

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

}

transformed data{

  // int n_sampled;

  // n_sampled = 1;


}

parameters{

vector<lower = 0, upper = 10>[nt]  f_t; //  fishing mortality

real<lower = 0> sigma_f; // standard deviation of fishing mortality

vector[nt]  rec_dev_t; //  recruitment deviates

real<lower = 0> sigma_r; // standard deviation of recruitment deviates

} // close parameters block

transformed parameters{

  real ssb0;

  real ssb_temp;

  real temp_rec_dev;

  matrix[nt, n_ages] n_ta; // numbers at time and age

  matrix[nt, n_ages] b_ta; // biomass at time and age

  matrix[nt, n_ages] ssb_ta; // spawning stock biomass at time and age

  matrix[nt, n_ages] c_ta; // catch (biomass) at time and age

  matrix[nt, n_ages] cn_ta; // catch (numbers) at time and age

  matrix[nt, n_lbins] p_lbin_sampled; // numbers at time and length bin

  vector[nt] c_t; // total catch at time step

  matrix[nt,2] profit_shock; // calculate profit

  row_vector[n_ages] temp_n_a;

  row_vector[n_ages] p_age_sampled;

  // fill matrices with zeros
  n_ta = rep_matrix(0,nt, n_ages);

  b_ta = rep_matrix(0,nt, n_ages);

  ssb_ta = rep_matrix(0,nt, n_ages);

  c_ta = rep_matrix(0,nt, n_ages);

  cn_ta = rep_matrix(0,nt, n_ages);

  p_lbin_sampled = rep_matrix(0,nt, n_lbins);

  profit_shock = rep_matrix(0, nt, 2);

  ssb0 = sum((r0 * exp(-m * (ages - 1)))  .* mean_maturity_at_age .* mean_weight_at_age); // virgin ssb

  temp_n_a = r0 * exp(-m * (ages' - 1)); // transpose to specify row format

  n_ta[1, 1:n_ages] = temp_n_a;  // add in virign recruitment to get population moving

  n_ta[1, n_ages] = (n_ta[1, n_ages - 1] * exp(-m)) / (1 - exp(-m));  //plus group

  b_ta[1, 1:n_ages] = n_ta[1, 1:n_ages] .* mean_weight_at_age'; // virgin biomass

  ssb_ta[1, 1:n_ages] = b_ta[1, 1:n_ages] .* mean_maturity_at_age'; // virgin ssb

  // order of events spawn, grow and die, recruit

  for (t in 2:nt){

  temp_rec_dev = 1;

  if (estimate_recruits == 1){

    temp_rec_dev = exp(rec_dev_t[t] - sigma_r^2/2);

  } // close estimate recruits

    ssb_temp =  sum(ssb_ta[t - 1, 1:n_ages]);

    n_ta[t,1] = ((0.8 * r0 * h *ssb_temp) / (0.2 * ssb0 * (1 - h) + (h - 0.2) * ssb_temp)) * temp_rec_dev; //calculate recruitment

    n_ta[t , 2:n_ages] = n_ta[t - 1, 1:(n_ages -1)] .* (exp(-(m + f_t[t - 1] * mean_selectivity_at_age[1:(n_ages - 1)])))'; // grow and die

    n_ta[t, n_ages] = n_ta[t - 1, n_ages] .* (exp(-(m + f_t[t - 1] * mean_selectivity_at_age[n_ages])))'; // assign to plus group

  b_ta[t, 1:n_ages] = n_ta[t, 1:n_ages] .* mean_weight_at_age'; // calculate biomass at age

  ssb_ta[t, 1:n_ages] = b_ta[t, 1:n_ages] .* mean_maturity_at_age'; // calculate ssb at age

   cn_ta[t - 1, 1:n_ages] = ((f_t[t - 1] * mean_selectivity_at_age) ./ (m + f_t[t - 1] * mean_selectivity_at_age))' .* n_ta[t - 1, 1:n_ages] .* (1 - exp(-(m + f_t[t - 1] * mean_selectivity_at_age)))'; // .* mean_weight_at_age';

   c_ta[t-1, 1:n_ages] = cn_ta[t - 1, 1:n_ages] .* mean_weight_at_age';

   c_t[t - 1] = sum(c_ta[t-1, 1:n_ages]);

  // run economic model


  // sample lengths

  p_age_sampled = cn_ta[t - 1, 1:n_ages] / sum(cn_ta[t - 1, 1:n_ages]);

  p_lbin_sampled[t - 1, 1:n_lbins] = (p_age_sampled * length_at_age_key) / sum(p_age_sampled * length_at_age_key);  // calculate the proportion of the sampled catch in each length bin


  } // close time loop


// fill in final time step
     cn_ta[nt, 1:n_ages] = ((f_t[nt] * mean_selectivity_at_age) ./ (m + f_t[nt] * mean_selectivity_at_age))' .* n_ta[nt, 1:n_ages] .* (1 - exp(-(m + f_t[nt] * mean_selectivity_at_age)))'; //

   c_ta[nt, 1:n_ages] = cn_ta[nt, 1:n_ages] .* mean_weight_at_age';

   c_t[nt] = sum(c_ta[nt, 1:n_ages]);

  p_age_sampled = cn_ta[nt, 1:n_ages] / sum(cn_ta[nt, 1:n_ages]);

  p_lbin_sampled[nt, 1:n_lbins] = p_age_sampled * length_at_age_key / sum(p_age_sampled * length_at_age_key);

}

model{

matrix[n_lcomps, n_lbins] pn_tl; // numbers at time and length bin

int temp_n_tl[n_lbins];

//// length comps likelihood ////

pn_tl = rep_matrix(0,n_lcomps, n_lbins);

for (i in 1:(n_lcomps - 1)){

  length_comps[i, 1:n_lbins] ~ multinomial(to_vector(p_lbin_sampled[length_comps_years[i], 1:n_lbins]));

} // close length likelihood

//// f likelihood ////

f_t[1] ~ normal(0,0.5);

for (t in 2:nt){

f_t[t] ~ normal(f_t[t - 1], sigma_f);

} // close f likelihood loop

sigma_f ~ normal(.1,3);

//// recruitment likelihood ////

rec_dev_t ~ normal(0, sigma_r);

sigma_r ~ normal(0.7, 0.2);


} // close model block

generated quantities{


int n_tl[nt, n_lbins];

for (t in 1:nt){
    n_tl[t, 1:n_lbins] = multinomial_rng(to_vector(p_lbin_sampled[t, 1:n_lbins]),100); // generate length comp samples
}

} // close generated quantities block

