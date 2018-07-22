#include <tuple>
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{

  ///////////////////////////
  /////// data inputs ////////
  ///////////////////////////

  // options and misc

  DATA_INTEGER(nt); // number of time steps to run

  DATA_INTEGER(n_burn); // number of burn-in time steps

  DATA_INTEGER(n_ages); // number of ages

  DATA_INTEGER(n_lbins); // number of length bins;

  DATA_INTEGER(n_lcomps); // number of timesteps of length comps available

  DATA_VECTOR(ages); // vector of ages

  DATA_INTEGER(economic_model); // 0 = no bioeconomic prior, 1 = bioeconomic prior,

  // biological and length data

  DATA_IMATRIX(length_comps); // length compositiond data

  DATA_IVECTOR[length_comps_years]; // time steps in which length comps are available

  DATA_VECTOR(mean_length_at_age);

  DATA_VECTOR(mean_weight_at_age);

  DATA_VECTOR(mean_maturity_at_age);

  DATA_MATRIX(length_at_age_key);

  DATA_SCALAR(m); //natural mortality

  DATA_SCALAR(h); //steepness

  DATA_SCALAR(r0); // virgin recruitment

  DATA_SCALAR(k); //vbk growth

  DATA_SCALAR(loo); // vbk linf

  DATA_SCALAR(t0); // t0

  DATA_SCALAR(sigma_r_guess);

  DATA_SCALAR(sd_sigma_r);

  // economic data

  DATA_SCALAR(length_50_sel_guess);

  DATA_SCALAR(delta_guess);

  DATA_VECTOR(q_t);

  ///////////////////////////
  /////// parameters ////////
  ///////////////////////////

  PARAMETER(log_burn_effort); // log burn-in effort

  PARAMETER_VECTOR(log_effort_t); // log effort in time t

  PARAMETER_VECTOR(exp_rec_dev_t); // exp recruitment deviates in time t

  PARAMETER(log_sigma_r);

  PARAMETER(log_sigma_effort);

  ///////////////////////////
  ///////  model     ////////
  ///////////////////////////

  Type burn_f = q_t(1) * exp(log_burn_effort);

  vector<Type> effort_t = exp(log_effort_t);

  Type sigma_r = exp(log_sigma_r);

  vector<Type> rec_dev_t = exp(uc_rec_dev_t);

  vector<Type> mean_selectivity_at_age = 1 / (1 + exp(-log(Type(19)) * ((mean_length_at_age - length_50_sel_guess) / delta_guess))); // selectivity ogive at age

  // create some storage //

  matrix<Type> n_a_init(n_burn, n_ages); // numbers at time and age

  n_a_init = n_a_init.setZero();

  matrix<Type> ssb_init(n_burn, n_ages); // numbers at time and age

  ssb_init = ssb_init.setZero();

  matrix<Type> n_ta(nt, n_ages);

  n_ta = n_ta.setZero();

  matrix<Type> ssb_ta(nt, n_ages);

  ssb_ta = ssb_ta.setZero();

  matrix<Type> c_ta(nt, n_ages);

  c_ta = c_ta.setZero();

  matrix<Type> cn_ta(nt, n_ages);

  cn_ta = cn_ta.setZero();

  matrix<Type> p_lbin_sampled(nt, n_lbins);

  p_lbin_sampled = p_lbin_sampled.setZero();

  // burn in initial population //

  for (int a = 0; a<n_ages, a++){

  n_a_init(0,a) = r0 * exp(-m * a);

  ssb_init(0,a) = n_a_init(0,a) * mean_weight_at_age(a) * mean_maturity_at_age(a);

} //close virgin population loop

n_a_init(0, n_ages + 1) = n_a_init(0, n_ages) * exp(-m) / (1 - exp(-m));  //plus group

ssb_init(0,n_ages + 1) = n_a_init(0,n_ages + 1) * mean_weight_at_age(n_ages + 1) * mean_maturity_at_age(n_ages + 1);

Type ssb0 = ssb_init.row(0).sum(); // unfished ssb;

for (int t = 1; t < n_burn; t++){

  Type ssb_temp =  ssb_init.row(t - 1).sum();

  n_a_init(t,1) = (0.8 * r0 * h *ssb_temp) / (0.2 * ssb0 * (1 - h) + (h - 0.2) * ssb_temp); //calculate recruitment

  for (int a = 1; a<n_ages, a++){

   n_a_init(t,a) = n_a_init(t - 1, a - 1) * (exp(-(m + burn_f* mean_selectivity_at_age(a - 1))));

   ssb_init(t,a) = n_a_init(t,a) * mean_weight_at_age(a) * mean_maturity_at_age(a);

   }

   n_a_init(t, n_ages + 1) = n_a_init(t, n_ages + 1) + n_a_init(t - 1, n_ages + 1) * (exp(-(m + burn_f * mean_selectivity_at_age(n_ages + 1)))); // assign to plus group

   ssb_init(t,n_ages + 1) = n_a_init(t,n_ages + 1) * mean_weight_at_age(n_ages + 1) * mean_maturity_at_age(n_ages + 1);

} // close burn-in period

// start observed period //


for (int a = 0; a<n_ages, a++){

n_ta(0, a) = n_a_init(n_burn + 1, a);  // start at initial depletion

n_ta(0,0) = n_ta(0,0) * rec_dev_t(0); // initial recruitment deviate

ssb_ta(0, a) = n_ta(0, a) * mean_maturity_at_age(a) * mean_weight_at_age(a); // starting ssb

} // close filling in intial Population

// order of events spawn, grow and die, recruit

for (int t = 1; t < nt; t++) {


  ssb_temp =  ssb_ta.row(t - 1).sum();

   n_ta(t,0) = ((0.8 * r0 * h *ssb_temp) / (0.2 * ssb0 * (1 - h) + (h - 0.2) * ssb_temp)) * rec_dev_t(t); //calculate recruitment

   for (int a = 1; a<n_ages, a++){

     n_ta(t,a) = n_ta(t - 1, a - 1) * (exp(-(m + f_t(t - 1) * mean_selectivity_at_age(a - 1))));

     ssb_ta(t,a) = n_ta(t,a) * mean_weight_at_age(a) * mean_maturity_at_age(a);
   } // close grow and die

   n_ta(t, n_ages + 1) = n_ta(t, n_ages + 1) + n_ta(t - 1, n_ages + 1) * (exp(-(m + f_t(t - 1) * mean_selectivity_at_age(n_ages + 1)))); // assign to plus group

   ssb_ta(t,n_ages + 1) = n_ta(t,n_ages + 1) * mean_weight_at_age(n_ages + 1) * mean_maturity_at_age(n_ages + 1);

   // fill in economic data

   for (int a = 0; a<n_ages, a++){

     cn_ta(t - 1,a) = ((f_t[t - 1] * mean_selectivity_at_age(a)) / (m + f_t[t - 1] * mean_selectivity_at_age(a))) * n_ta(t - 1, a) .* (1 - exp(-(m + f_t(t - 1) * mean_selectivity_at_age(a))));

     c_ta(t-1, a) = cn_ta(t - 1, a) .* mean_weight_at_age(a);
   } // close catch age loop

  c_t(t - 1) = c_ta.row(t-1).sum();

 profit_t(t - 1) = price_t(t - 1) * c_t(t - 1) - cost_t(t - 1) * pow(effort_t(t - 1),beta);

 ppue_hat_t(t - 1) = profit_t(t - 1) / effort_t(t - 1);

 f_t(t) = effort_t(t) * q_t(t);

 // sample lengths //

Type total_caught = cn_ta.row(t - 1).sum();

vector<Type> p_age_sampled(n_ages);

p_age_sampled.setZero();

for (int a = 0; a<n_ages, a++){

p_age_sampled(a) = cn_ta(t - 1, a) / total_caught;

}


 p_lbin_sampled[t - 1, 1:n_lbins] = (p_age_sampled * length_at_age_key) / sum(p_age_sampled * length_at_age_key);  // calculate the proportion of the sampled catch in each length bin

} // close age loop

n_ta[t, n_ages] = n_ta[t, n_ages] + n_ta[t - 1, n_ages] .* (exp(-(m + f_t[t - 1] * mean_selectivity_at_age[n_ages])))'; // assign to plus group


}

} // close TMB
