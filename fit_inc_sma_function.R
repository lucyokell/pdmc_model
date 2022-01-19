##### Function for recalibrating to SMA incidence
fit_inc_sma<-function(inc_sma_genpop,
                      target_tot_inc_sma,
                      eir,prop_hosp0,prop_hosp_paton0,
                      r_eir=r_eir,
                      prob_um_plac=prob_um_plac,
                      prob_um_pmc=prob_um_pmc,
                      max_beta_phi=max_beta_phi,
                      w_scale_risk=w_scale_risk,
                      w_shape_risk=w_shape_risk) {
  curr_inc<-run_repeat_mod(inc_sma_genpop, 
                           eir,pmc=F,
                           ft_um_pd=0.5,ft_sev_pd = 0.5,prop_hosp=prop_hosp0,
                           prop_hosp_paton=prop_hosp_paton0,
                           output="summary2R",
                           r_eir=r_eir,
                           prob_um_plac=prob_um_plac,
                           prob_um_pmc=prob_um_pmc,
                           max_beta_phi=max_beta_phi,
                           w_scale_risk=w_scale_risk,
                           w_shape_risk=w_shape_risk,
                           adm="test")
  ind<-which(names(curr_inc)=="inc_sma_all_paton_eq")
  return((1000*target_tot_inc_sma - 1000*curr_inc[ind])^2)
}


optimise_fit_sma<-function(filestart,adm,inc_sma_genpop0,
                           prop_hosp_paton0,eir0,
                           prop_hosp0,r_eir,prob_um_plac,
                           prob_um_pmc,max_beta_phi,w_scale_risk,w_shape_risk) {
  x<-optimise(f=fit_inc_sma, interval=c(0,inc_sma_genpop0/prop_hosp_paton0),
           target_tot_inc_sma = inc_sma_genpop0,eir=eir0,
           prop_hosp0=prop_hosp0,prop_hosp_paton0=prop_hosp_paton0,
           r_eir=r_eir,
           prob_um_plac=prob_um_plac,
           prob_um_pmc=prob_um_pmc,
           max_beta_phi=max_beta_phi,
           w_scale_risk=w_scale_risk,
           w_shape_risk=w_shape_risk,
           tol=1e-4)
  y<-x$minimum
  saveRDS(y,
          file=paste0(filestart,adm,".rds"))
}
