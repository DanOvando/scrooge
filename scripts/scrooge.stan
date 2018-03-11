data{

//// run parameters ////

int<lower = 1> nt; // number of time steps to run

int n_ages; // number of ages

int n_lbins; // number of length bins;

vector[n_ages] ages;

//// length data ////

//// economic data ////

vector<lower=0, upper=1>[n_ages] mean_selectivity_at_age;


//// biology ////


vector<lower=0>[n_ages] mean_length_at_age;

vector<lower=0>[n_ages] mean_weight_at_age;

vector<lower=0, upper=1>[n_ages] mean_maturity_at_age;

real m; //natural mortality

real h; //steepness

real r0; // virgin recruitment

real k; //vbk growth

real loo; // vbk linf

real t0; // t0

}

transformed data{


}

parameters{

vector<lower=0, upper=3>[nt]  f_t; //  fishing mortality

real<lower = 0> sigma_f; // standard deviation of fishing mortality

}

transformed parameters{

  real ssb0;

  real ssb_temp;

  matrix[nt, n_ages] n_ta; // numbers at time and age

  matrix[nt, n_ages] b_ta; // biomass at time and age

  matrix[nt, n_ages] ssb_ta; // spawning stock biomass at time and age

  matrix[nt, n_ages] c_ta; // catch (biomass) at time and age

  matrix[nt, n_lbins] n_tl; // numbers at time and length bin

  row_vector[n_ages] temp_n_a;


  // fill matrices with zeros
  n_ta = rep_matrix(0,nt, n_ages);

  b_ta = rep_matrix(0,nt, n_ages);

  ssb_ta = rep_matrix(0,nt, n_ages);

  c_ta = rep_matrix(0,nt, n_ages);

  n_tl = rep_matrix(0,nt, n_lbins);

  ssb0 = sum((r0 * exp(-m * (ages - 1)))  .* mean_maturity_at_age .* mean_weight_at_age); // virgin ssb

  temp_n_a = r0 * exp(-m * (ages' - 1)); // transpose to specify row format

  n_ta[1, 1:n_ages] = temp_n_a;  // add in virign recruitment to get population moving

  b_ta[1, 1:n_ages] = n_ta[1, 1:n_ages] .* mean_weight_at_age'; // virgin biomass

  ssb_ta[1, 1:n_ages] = b_ta[1, 1:n_ages] .* mean_maturity_at_age'; // virgin ssb

  // order of events spawn, grow and die, recruit

  for (t in 2:nt){

    ssb_temp =  sum(ssb_ta[t - 1, 1:n_ages]);

    n_ta[t,1] = (0.8 * r0 * h *ssb_temp) / (0.2 * ssb0 * (1 - h) + (h - 0.2) * ssb_temp);

    // n_ta[t,1] = (4 * h * r0 * ssb_temp) / (ssb0 * (1 - h) + ssb_temp * (5*h - h - 1)); // add in recruits

    n_ta[t , 2:n_ages] = n_ta[t - 1, 1:(n_ages -1)] .* (exp(-(m + f_t[t - 1] * mean_selectivity_at_age[1:(n_ages - 1)])))'; // grow and die

    n_ta[t, n_ages] = n_ta[t - 1, n_ages] .* (exp(-(m + f_t[t - 1] * mean_selectivity_at_age[n_ages])))'; // assign to plus group

  b_ta[t, 1:n_ages] = n_ta[t, 1:n_ages] .* mean_weight_at_age'; // virgin biomass

  ssb_ta[t, 1:n_ages] = b_ta[t, 1:n_ages] .* mean_maturity_at_age'; // virgin ssb

   c_ta[t - 1, 1:n_ages] = ((f_t[t - 1] * mean_selectivity_at_age) ./ (m + f_t[t - 1] * mean_selectivity_at_age))' .* n_ta[t - 1, 1:n_ages] .* (1 - exp(-(m + f_t[t - 1] * mean_selectivity_at_age)))' .* mean_weight_at_age';


  } // close time loop

     c_ta[nt, 1:n_ages] = ((f_t[nt] * mean_selectivity_at_age) ./ (m + f_t[nt] * mean_selectivity_at_age))' .* n_ta[nt, 1:n_ages] .* (1 - exp(-(m + f_t[nt] * mean_selectivity_at_age)))' .* mean_weight_at_age';


}

model{

f_t ~ normal(0, 0.1);

}

generated quantities{



}
