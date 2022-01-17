
functions {  // runs model for one site.
  real[,] run_model(real EIR_y, real prob_um_plac, real prob_um_pmc, real w_scale_risk, 
  real w_shape_risk, 
    real max_beta_phi, real r, int runTime, int N_plac_init,
    int N_pmc_init, real[] prop_fup_plac, real[] prop_fup_pmc, 
    real[] prob_rp, int[] pmc_times,
    real[] p_protect_dp) {
    
  real res[runTime,4];
  real prob_sm_plac = 1-prob_um_plac;
  real prob_sm_pmc = 1-prob_um_pmc;
  real beta_phi[runTime];
  real current_pp_dp;
  real prob_symp_d[runTime];
  real al_eff;  // AL efficacy
  int dur_p;
  real current_prob_um;
  int dur_sev;
  real return_sev_plac;
  real return_um_plac;
  real return_sev_pmc;
  real return_um_pmc;

  ////// Placebo groups
  real Rd_plac[runTime+1];
  real S_plac[runTime+1];
  real Sm_plac[runTime+1];
  real Um_plac[runTime+1];
  real N_check_plac[runTime+1];

    /////// PMC group
  real Rd_pmc[runTime+1];
  real S_pmc[runTime+1];
  real Sm_pmc[runTime+1];
  real Um_pmc[runTime+1];
  real N_check_pmc[runTime+1];

    /// Fixed AL efficacy
  al_eff=0.98;
  dur_p=15;  // 2 days treatment seeking + 13 days proph.
  dur_sev=dur_p+3;  // dur_p + 3 days median hospital
  /// INITIAL CONDITIONS
    for(i in 1:runTime) {
      beta_phi[i]=max_beta_phi*exp(-r*EIR_y)*exp(-((i/w_scale_risk)^w_shape_risk));
      prob_symp_d[i]=1-exp(-beta_phi[i]*EIR_y/365);
    }
    
      
    // Initial conditions
    // Placebo initial conditions
    Rd_plac[1] = N_plac_init - prob_rp[1]*N_plac_init;
    S_plac[1]=prob_rp[1]*N_plac_init;
    Sm_plac[1]=0;
    Um_plac[1]=0;


    // PMC initial conditions
    Rd_pmc[1] = N_pmc_init - prob_rp[1]*N_pmc_init;
    S_pmc[1]=prob_rp[1]*N_pmc_init;
    Sm_pmc[1]=0;
    Um_pmc[1]=0;

    N_check_plac[1]=N_plac_init;
    N_check_pmc[1]=N_pmc_init;
  
    

    for(i in 2:(runTime+1)) {



    /// Run placebo arm
      Rd_plac[i] = prop_fup_plac[i-1]*(Rd_plac[i-1] - prob_rp[i-1]*Rd_plac[1]);
      if(i<=dur_p) {
        return_um_plac=0;
      } else {
        return_um_plac=prob_symp_d[i-dur_p]*al_eff*S_plac[i-dur_p]*prob_um_plac;
      }
      
      if(i<=dur_sev) {
        return_sev_plac=0;
      } else {
        return_sev_plac=prob_symp_d[i-dur_sev]*al_eff*S_plac[i-dur_sev]*prob_sm_plac;
      }

      Sm_plac[i] = prop_fup_plac[i-1]*(Sm_plac[i-1] + prob_symp_d[i-1]*al_eff*S_plac[i-1]*prob_sm_plac - return_sev_plac);
      Um_plac[i] = prop_fup_plac[i-1]*(Um_plac[i-1] + prob_symp_d[i-1]*al_eff*S_plac[i-1]*prob_um_plac - return_um_plac);
      S_plac[i] = prop_fup_plac[i-1]*(S_plac[i-1] + prob_rp[i-1]*Rd_plac[1] - prob_symp_d[i-1]*al_eff*S_plac[i-1] + return_sev_plac + return_um_plac);
        
            
      res[i-1,1] = prob_symp_d[i-1]*S_plac[i-1]*prob_sm_plac;
      res[i-1,3] = prob_symp_d[i-1]*S_plac[i-1]*prob_um_plac;
      
      N_check_plac[i] = Rd_plac[i]+S_plac[i]+Um_plac[i]+Sm_plac[i];
      
      //print(N_check_plac[i]);
                  
      /// Run PMC arm
     Rd_pmc[i] = prop_fup_pmc[i-1]*(Rd_pmc[i-1] - prob_rp[i-1]*Rd_pmc[1]);
          
    if(i<pmc_times[1]) {
      
      if(i<=dur_p) {
        return_um_pmc=0;
      } else {
        return_um_pmc=prob_symp_d[i-dur_p]*al_eff*S_pmc[i-dur_p]*prob_um_plac;
      }
      
      if(i<=dur_sev) {
        return_sev_pmc=0;
      } else {
        return_sev_pmc=prob_symp_d[i-dur_sev]*al_eff*S_pmc[i-dur_sev]*prob_sm_plac;
      }

      S_pmc[i] = prop_fup_pmc[i-1]*(S_pmc[i-1] + prob_rp[i-1]*Rd_pmc[1] - prob_symp_d[i-1]*al_eff*S_pmc[i-1] + return_sev_pmc + return_um_pmc);
      Um_pmc[i] = prop_fup_pmc[i-1]*(Um_pmc[i-1] + prob_symp_d[i-1]*al_eff*S_pmc[i-1]*prob_um_plac - return_um_pmc);
      Sm_pmc[i] = prop_fup_pmc[i-1]*(Sm_pmc[i-1] + prob_symp_d[i-1]*al_eff*S_pmc[i-1]*prob_sm_plac - return_sev_pmc);
          
      res[i-1,2] = prob_symp_d[i-1]*S_pmc[i-1]*prob_sm_plac;
      res[i-1,4] = prob_symp_d[i-1]*S_pmc[i-1]*prob_um_plac;
  
      } else {
          // p_protect_dp is adjusted for adherence to each dose in Kwambai et al.
          if(i<pmc_times[2]) current_pp_dp=p_protect_dp[i-pmc_times[1]+1];
          else if(i>=pmc_times[2] && i< pmc_times[3]) current_pp_dp= 0.985*p_protect_dp[i-pmc_times[2]+1];
          else current_pp_dp=0.973*p_protect_dp[i-pmc_times[3]+1];
          
          if(current_pp_dp>0.01) current_prob_um=prob_um_pmc;
          else current_prob_um=prob_um_plac;
        
        
          if(i<=dur_p) {
            return_um_pmc=0;
          } else {
            return_um_pmc=prob_symp_d[i-dur_p]*al_eff*S_pmc[i-dur_p]*(1-current_pp_dp)*current_prob_um;

          }
          
          if(i<=dur_sev) {
            return_sev_pmc=0;
          } else {
            return_sev_pmc=prob_symp_d[i-dur_sev]*al_eff*S_pmc[i-dur_sev]*(1-current_pp_dp)*(1-current_prob_um);
          }

          
          S_pmc[i] = prop_fup_pmc[i-1]*(S_pmc[i-1] + prob_rp[i-1]*Rd_pmc[1] - prob_symp_d[i-1]*al_eff*S_pmc[i-1]*(1-current_pp_dp) + return_um_pmc +return_sev_pmc);
          Um_pmc[i] = prop_fup_pmc[i-1]*(Um_pmc[i-1] + prob_symp_d[i-1]*al_eff*S_pmc[i-1]*(1-current_pp_dp)*current_prob_um - return_um_pmc);
          Sm_pmc[i] = prop_fup_pmc[i-1]*(Sm_pmc[i-1] + prob_symp_d[i-1]*al_eff*S_pmc[i-1]*(1-current_pp_dp)*(1-current_prob_um) - return_sev_pmc);

          res[i-1,2] = prob_symp_d[i-1]*S_pmc[i-1]*(1-current_pp_dp)*(1-current_prob_um);
          res[i-1,4] = prob_symp_d[i-1]*S_pmc[i-1]*(1-current_pp_dp)*current_prob_um;
        } // end of if else clauses.
        
        N_check_pmc[i] = Rd_pmc[i]+S_pmc[i]+Um_pmc[i] + Sm_pmc[i];
        //print(i, " ", N_check_pmc[i]);
      
      } // end of time loop
    
    return(res);

  } // end of run model function

}

data{
  int nSites;
  int N_plac_init[nSites];
  int N_pmc_init[nSites];
  int runTime;
  real prop_fup_plac[runTime,nSites];
  real prop_fup_pmc[runTime,nSites];
  // 2 trt groups x 2 outcomes (UM,SM) x nSites
  int daily_events[(runTime-14),4,nSites]; //start fitting at day 14
  
  real shape_EIR_y;
  real rate_EIR_y[nSites];
  real prob_rp[runTime];
  int pmc_times[3];
  int length_pp_dp;
  real p_protect_dp[length_pp_dp];
  int nPred;
  real EIRpred[nPred];
  real prop_fup_pred[runTime];

}


parameters{
  real<lower=0.00001, upper=100> EIR_y[nSites];
  real<lower=0.00001, upper=1> prob_um_plac;
  real<lower=0.00001, upper=1> prob_um_pmc;
  real<lower=0.00001, upper=2> max_beta_phi;
  real<lower=0.0001> r; 
  real<lower=0.0001> w_scale_risk;
  real<lower=0.0001> w_shape_risk;
  real res[runTime,4,nSites];

}


// transformed parameters{
//   
//   real res[runTime,4,nSites];  // dimensions = days x 4 outcomes: SM plac, SM PMC, UM plac, UM PMC. x nSites
//   //vector[8] res0;
//   
//   for(s in 1:nSites) {
//     res[,,s]=run_model(EIR_y[s], prob_um_plac, prob_um_pmc, w_scale_risk, w_shape_risk,  
//       max_beta_phi, r, runTime, N_plac_init[s],
//       N_pmc_init[s], prop_fup_plac[,s], prop_fup_pmc[,s], prob_rp, pmc_times,
//       p_protect_dp);
//   }
// }



model{
  // priors
  EIR_y ~ gamma(shape_EIR_y, rate_EIR_y);
  r ~ normal(0, 1);
  w_scale_risk ~ normal(0, 500);
  w_shape_risk ~ normal(0, 500);
  
  
  // likelihood
  for(s in 1:nSites) {
    for(j in 1:4) {
      daily_events[,j,s] ~ poisson(res[15:runTime,j,s]);
    }
  }
}

generated quantities {
     real res2[runTime,4,nPred];
     vector[8] events_tot;
  
       for(j in 1:8) events_tot[j]=0.0;
    
          //  add up total events of each type across days and sites for the two time periods
    //   // res order of results: SM plac, SM PMC, UM plac, UM PMC.
    //   // events_tot order of results: SM plac3, SM PMC3, UM plac3, UM PMC3, SM plac6, SM PMC6, UM plac6, UM PMC6.
       for(s in 1:nSites) {
           for(d in 15:104) {
               for(j in 1:4) {
                    events_tot[j] += res[d,j,s];
                 }
             }
           for(d in 105:runTime) {
               for(j in 1:4) {
                    events_tot[j+4] += res[d,j,s];
                 }
             }
         }
    
           for(s in 1:nPred) {
               res2[,,s]=run_model(EIRpred[s], prob_um_plac, prob_um_pmc, w_scale_risk, w_shape_risk,
                                            max_beta_phi, r, runTime, 100,
                                            100, prop_fup_pred, prop_fup_pred,
                                            prob_rp, pmc_times,
                                            p_protect_dp);
             }
     }
