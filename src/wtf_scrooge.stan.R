/*scrooge_v3.0

effort prior is open-access dynamics

non-centered recruitment and effort deviations

*/


  data{

    //// run parameters ////

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

      row_vector[nt] q_t;

      real beta;

      real length_50_sel_guess;

      real delta_guess;

      real <lower = 0> max_cost_guess;

      real<lower = 0> p_response_guess;

      real<lower = 0> cv_effort;

      real<lower = 0, upper  = 1> effort_data_weight;


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

  vector[nt]  log_effort_t; // effort in time t

  //vector<lower = 0.001, upper = 15>[nt]  effort_t; // effort in time t

  vector[nt]  uc_rec_dev_t; //  recruitment deviates

  // real <lower = 0> sigma_r; // standard deviation of recruitment deviates

  // real <lower = 0> sigma_effort; // standard deviation of effort process error

  // real<lower = 0, upper = .9> p_length_50_sel; // length at 50% selectivity

  // real<lower = 0, upper = 10> p_response;

  // real<lower = 0> max_cost;

} // close parameters block

transformed parameters{

  real ssb0;

  real ssb_temp;

  real sel_delta;

  real length_50_sel;

  vector[nt]  effort_t; //  base effort multiplier

  vector[nt]  cost_t; //  base effort multiplier

  vector[nt]  rec_dev_t; //  base effort multiplier

  matrix[n_burn, n_ages] n_a_init; // numbers at time and age

  matrix[n_burn, n_ages] ssb_init; // numbers at time and age

  matrix[nt, n_ages] n_ta; // numbers at time and age

  matrix[nt, n_ages] ssb_ta; // spawning stock biomass at time and age

  matrix[nt, n_ages] c_ta; // catch (biomass) at time and age

  matrix[nt, n_ages] cn_ta; // catch (numbers) at time and age

  matrix[nt, n_lbins] p_lbin_sampled; // numbers at time and length bin

  vector[nt] c_t; // total catch at time step

  vector[nt] f_t; // fishing mortality at time step

  vector[nt] profit_t; // profits

  vector[n_ages] mean_selectivity_at_age; // mean selectivity at age

  row_vector[n_ages] p_age_sampled;

  // length_50_sel = loo * p_length_50_sel;

  length_50_sel =  length_50_sel_guess;

  sel_delta = 2;

  cost_t = max_cost_guess * relative_cost_t';

  effort_t = exp(log_effort_t);

  // effort_t = exp(log_effort_t - sigma_effort^2/2);

  //sigma_r = exp(log_sigma_r);

  rec_dev_t = exp(sigma_r_guess * uc_rec_dev_t - sigma_r_guess^2/2);

  // fill matrices with zeros //

  n_ta = rep_matrix(0,nt, n_ages);

  ssb_ta = rep_matrix(0,nt, n_ages);

  c_ta = rep_matrix(0,nt, n_ages);

  cn_ta = rep_matrix(0,nt, n_ages);

  p_lbin_sampled = rep_matrix(0,nt, n_lbins);

  mean_selectivity_at_age = 1.0 ./ (1 + exp(-log(19) * ((mean_length_at_age - length_50_sel) / sel_delta))); // selectivity ogive at age

  // set up initial population //

  f_t[1] = effort_t[1] .* q_t[1];

  n_a_init[1,1:n_ages] = r0 * exp(-m * (ages - 1))';

  n_a_init[1, n_ages] = n_a_init[1, n_ages - 1] * exp(-m) / (1 - exp(-m));  //plus group

  ssb_init[1, 1:n_ages] = n_a_init[1, 1:n_ages] .* mean_maturity_at_age'.* mean_weight_at_age'; // calculate ssb at age

  ssb0 = sum(ssb_init[1, 1:n_ages]);

  for (t in 2:n_burn){

    ssb_temp =  sum(ssb_init[t - 1, 1:n_ages]);

    n_a_init[t,1] = (0.8 * r0 * h *ssb_temp) / (0.2 * ssb0 * (1 - h) + (h - 0.2) * ssb_temp); //calculate recruitment

    n_a_init[t , 2:n_ages] = n_a_init[t - 1, 1:(n_ages -1)] .* (exp(-(m + f_t[1] * mean_selectivity_at_age[1:(n_ages - 1)])))'; // grow and die

    n_a_init[t, n_ages] = n_a_init[t, n_ages] + n_a_init[t - 1, n_ages] .* (exp(-(m + f_t[1] * mean_selectivity_at_age[n_ages])))'; // assign to plus group

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

    n_ta[t , 2:n_ages] = n_ta[t - 1, 1:(n_ages -1)] .* (exp(-(m + f_t[t - 1] * mean_selectivity_at_age[1:(n_ages - 1)])))'; // grow and die

    n_ta[t, n_ages] = n_ta[t, n_ages] + n_ta[t - 1, n_ages] .* (exp(-(m + f_t[t - 1] * mean_selectivity_at_age[n_ages])))'; // assign to plus group

    ssb_ta[t, 1:n_ages] = n_ta[t, 1:n_ages] .* mean_maturity_at_age'.* mean_weight_at_age'; // calculate ssb at age

    cn_ta[t - 1, 1:n_ages] = ((f_t[t - 1] * mean_selectivity_at_age) ./ (m + f_t[t - 1] * mean_selectivity_at_age))' .* n_ta[t - 1, 1:n_ages] .* (1 - exp(-(m + f_t[t - 1] * mean_selectivity_at_age)))'; // .* mean_weight_at_age';

    c_ta[t-1, 1:n_ages] = cn_ta[t - 1, 1:n_ages] .* mean_weight_at_age';

    c_t[t - 1] = sum(c_ta[t-1, 1:n_ages]);

    // run economic model

    profit_t[t - 1] = price_t[t - 1] * c_t[t - 1] - cost_t[t - 1] * effort_t[t - 1] ^ beta;

    f_t[t] = effort_t[t] * q_t[t];

    // sample lengths //

      p_age_sampled = cn_ta[t - 1, 1:n_ages] / sum(cn_ta[t - 1, 1:n_ages]);

    p_lbin_sampled[t - 1, 1:n_lbins] = (p_age_sampled * length_at_age_key) / sum(p_age_sampled * length_at_age_key);  // calculate the proportion of the sampled catch in each length bin


  } // close time loop


  // fill in final time step
  cn_ta[nt, 1:n_ages] = ((f_t[nt] * mean_selectivity_at_age) ./ (m + f_t[nt] * mean_selectivity_at_age))' .* n_ta[nt, 1:n_ages] .* (1 - exp(-(m + f_t[nt] * mean_selectivity_at_age)))'; //

    c_ta[nt, 1:n_ages] = cn_ta[nt, 1:n_ages] .* mean_weight_at_age';

  c_t[nt] = sum(c_ta[nt, 1:n_ages]);

  profit_t[nt] = price_t[nt] * c_t[nt] - cost_t[nt] * effort_t[nt] ^ beta;

  p_age_sampled = cn_ta[nt, 1:n_ages] / sum(cn_ta[nt, 1:n_ages]);

  p_lbin_sampled[nt, 1:n_lbins] = (p_age_sampled * length_at_age_key) / sum(p_age_sampled * length_at_age_key);

}

model{

real oa_prediction;

real effort_data_prediction;

real new_effort;

real previous_max;

//// length comps likelihood ////

for (i in 1:(n_lcomps)){

length_comps[i, 1:n_lbins] ~ multinomial(to_vector(p_lbin_sampled[length_comps_years[i], 1:n_lbins]));

} // close length likelihood

//// recruitment prior ////


uc_rec_dev_t ~ normal(0, 1);

} // close model block

generated quantities{


int n_tl[nt, n_lbins];

for (t in 1:nt){

n_tl[t, 1:n_lbins] = multinomial_rng(to_vector(p_lbin_sampled[t, 1:n_lbins]),500); // generate length comp samples
}

} // close generated quantities block


