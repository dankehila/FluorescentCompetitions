functions{
  //ODE system
  vector logistic(real t, //time
                         vector y, // state
                         //real r, //average growth rate
                         real K, //carrying capacity / maximum total density of both y
                         row_vector r // relative growth (dis)advantage of y1
                         //real r2 // relative growth (dis)advantage of y2
                        ) {
      vector[2] dydt;
      dydt[1] = (r[1])*y[1]*(1-(y[1]+y[2])/K);
      dydt[2] = (r[2])*y[2]*(1-(y[1]+y[2])/K);
    return dydt;
  }

}
data{
  int<lower=1> N; //rows in data
  int<lower=1> N_ID; //number of experiments
  array[N] int<lower=1,upper=N_ID> IDX; //experiment ID
  int<lower=1> N_y0; // how many y0 to expect (dilution factors, different days, different FPs, etc.)
  array[N_ID] int<lower=1> ID_y10; //inoculum information for competitor 1
  array[N_ID] int<lower=1> ID_y20; //inoculum information for competitor 2
  real t0; //initial timepoint (0)
  array[N] vector[2] yobs; //ob
  array[N] real ts; // timepoints
  array[N_ID] int<lower=1> start_time;
  array[N_ID] int<lower=1> end_time;
  
  //model for fitting selection coefficients
  int<lower=0> N_X; //number of cost components
  matrix[N_ID,N_X] X1;
   matrix[N_ID,N_X] X2;
  int<lower=0> N_pred;
  matrix[N_pred,N_X] X_pred; //design matrix for predictions
  matrix[N_pred,N_X] X_pred_no_epistasis;
  
  //prior information
  array[2] real prior_r;
  array[2] real prior_K;
  array[2] real prior_y0;
  real<lower=0> epsilon;
  int num_steps;
}

parameters{
  //parameters for fitting logistic ODE model
  real K; // carrying capacity / max density
  real<lower=0> sigma; //experimental noise hyperparameter
  //array[N_ID] vector[2]  r; //array of growth rates for each strain
  matrix<lower=0>[N_ID,2]  r;
  array[N_y0] real log_y0; //y0 is the same for both competitors (experimental constraint) but different for different dilutions
  
  //parameters for cost breakdown model
  vector[N_X] dc; //differential cost component
  real<lower=0> sigma_dc; //cost residual standard error
  real<lower=1> nu; //degrees of freedom for robust regression
  cholesky_factor_corr[2] L; //cholesky factor of correlation between competitors
} 
transformed parameters{
    //vector of multivariate means
   array[N_ID] vector[2] costs;
    for(i in 1:N_ID){
      costs[i,1] = X1[i]*dc;
      costs[i,2] = X2[i]*dc;
    }
    
     //identical sigmas for both selection coefficients (their position as reference/competitor is arbitrary)
  vector[2] L_Sigma = rep_vector(sigma_dc,2);

  // covariance matrix
  matrix[2,2] L_Cov = diag_pre_multiply(L_Sigma,L);
  
  //hierarchical maximum density: average density effect realized differently for each experiment
  array[N] vector[2] mu;
  real rbar = mean(r);//(mean(r1)+mean(r2))/2;
  array[N_ID] vector[2] s; 
    
  for(i in 1:N_ID){
    vector[2] y0_vec = exp([log_y0[ID_y10[i]],log_y0[ID_y20[i]]]');  //initialize y0 parameter vector based on index
    
    mu[start_time[i]:end_time[i],] = ode_rk45_tol(logistic,y0_vec,t0,ts[start_time[i]:end_time[i]],
                                                      epsilon,epsilon,num_steps,
                                                      exp(K),r[i,]); //numerical integration step
    
    s[i,] = [r[i,1]-rbar, r[i,2]-rbar]';//[r1[i]-rbar, r2[i]-rbar]';                                        
  }


 
}
model{
  //priors
  sigma ~ exponential(1); 
  K ~ normal(prior_K[1],prior_K[2]);
 // r ~ student_t(3,prior_r[1],prior_r[1]);
  for(i in 1:2){
    r[,i] ~ student_t(3,prior_r[1],prior_r[1]);
   // s ~  normal(0,0.1);  // expecting small differences from mean growth rate
  }
  log_y0 ~ normal(prior_y0[1],prior_y0[2]);

 // target+= reduce_sum(partial_sum,yobs,1,mu,sigma);

  //manual calculation of costs
  
  for(i in 1:2){
    target +=  lognormal_lpdf(yobs[,i]|log(mu[,i]),sigma); //likelihood for competitor y1
    //target +=  lognormal_lpdf(yobs[i,2]|log(mu[i,2]),sigma); //likelikhood for competitor y2
  }


  
  //priors for differential costs
  nu ~ gamma(2,0.1);
  dc ~ normal(0,0.1);
  L ~ lkj_corr_cholesky(2);
  sigma_dc ~ exponential(1);
  
  target+= multi_student_t_cholesky_lpdf(s|nu,costs,L_Cov);
    
 
 // target+=reduce_sum(partial_sum_mvt,s,1,costs,nu,L_Cov);
}
generated quantities {
  
  //convert cholesky matrix to correlation matrix
    corr_matrix[2] rho = multiply_lower_tri_self_transpose(L);
  
// log_lik for cross validation
vector[2*N] log_lik_ODE;
  {
  vector[N] log_lik_ODE1;
  vector[N] log_lik_ODE2;
// 
    for(i in 1:N) {
//     int start = ((IDX[i]-1)*culturing_time) + 1;
//     int end = IDX[i]*culturing_time;
//     vector[2] y0_vec = exp([log_y0[ID_y10[IDX[i]]],log_y0[ID_y20[IDX[i]]]]');  //initialize y0 parameter vector based on index
//     array[culturing_time] vector[2] mu = ode_rk45(logistic,y0_vec,t0,ts[start:end],r,K,s[IDX[i],1],s[IDX[i],2]); //numerical integration step
//     int position = i%culturing_time + 1;
    log_lik_ODE1[i] = lognormal_lpdf(yobs[i,1]|log(mu[i,1]),sigma); //likelihood for competitor y1
    log_lik_ODE2[i] = lognormal_lpdf(yobs[i,2]|log(mu[i,2]),sigma); //likelikhood for competitor y2
   }
   log_lik_ODE = append_row(log_lik_ODE1,log_lik_ODE2);
  }

  
    // array[N_ID] real log_lik_cost;
    // for(i in 1:N_ID){
    //   log_lik_cost[i] = multi_student_t_cholesky_lpdf(s[i] | nu, costs[i],L_Cov); 
    // }
    // 
  
      //posterior predictions
      
    //residuals: differences between s and costs
    array[N_ID] vector[2] resids;
    //y_rep: posterio prediction for diagnostics 
    array[N_ID] vector[2] y_rep;
    for(i in 1:N_ID){
      resids[i,] = s[i,] - costs[i,];
      y_rep[i,] = multi_student_t_cholesky_rng(nu,costs[i,],L_Cov);
    }
    
    
    //using X_pred matrix
    vector[N_pred] y_pred = X_pred*dc;
    array[N_pred] real l_pred = student_t_rng(nu,y_pred,sigma_dc);
    
    vector[N_pred] y_pred_no_epi = X_pred_no_epistasis*dc;
    array[N_pred] real l_pred_no_epi = student_t_rng(nu,X_pred_no_epistasis*dc,sigma_dc);
}

