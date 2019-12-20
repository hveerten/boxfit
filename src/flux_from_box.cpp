////////////////////////////////////////////////////////////////////////////////
//
// flux_from_box.cpp
//
// Created, Jun 23 2011, by HJvE
// Last Modified, July 29, 2016 by HJvE
//
////////////////////////////////////////////////////////////////////////////////

#include "flux_from_box.h"

////////////////////////////////////////////////////////////////////////////////

c_flux_box :: c_flux_box()
{
  res = 500; // some standard value for radiative transfer resolution
  
  // disable overruling of estimates
  ur_overrule = -1.;
  t0_overrule = -1.;
  t1_overrule = -1.;
  ur0_overrule = -1.;
  
  save_emission_profile = false;
  save_image = false;
}

////////////////////////////////////////////////////////////////////////////////

c_flux_box :: ~c_flux_box()
{}

////////////////////////////////////////////////////////////////////////////////

int c_flux_box :: set_flux_th(int i_box, double &F_th)
{
  int i_t; // dummy time loop index
  
  c_BM BM; // provide access to BM solution for purpose of computing time domain

  // variables used when determining distance between snapshots
  double Dt_sim; // under consideration
  double t_current, t_prev, t_next; // emission times used in loop
  double t0, t1; // range of times ('time domain')
  double Ftemp, Fold; // for bookkeeping purposes
  double taux; // auxiliary variable for computing emission time domain
  
  char log_filename[14];
  char image_filename[15];
  FILE *p_log_file; // log file
 
  //----------------------------------------------------------------------------
  // If a flux is requested that is not covered by the opening angles,
  // return 0 immediately

  if (i_box == -1)
  {
    F_th = 0.;
    return 0;
  }

  //----------------------------------------------------------------------------
  // set EDS boost

  // We require a single simulation grid boost for the multibox. We take the
  // boost from the first box in the multibox. All box boosts are equal
  #if BOOST_ == ENABLED_
    eds.boost = mbox.box[i_box].boost;
    eds.boostsqr = eds.boost * eds.boost;
    eds.beta_sim = sqrt(1.0 - 1.0 / (eds.boost * eds.boost));
  #endif // BOOST_

  //----------------------------------------------------------------------------
  // set emission time range and EDS size, either from argument of from estimate
  // Argument values < 0 mean that the estimate will be applied

  mbox.box[i_box].set_scale_factors();  
  
  if (t1_overrule < 0.)
  {
    if (mbox.box[i_box].jet != RECEDING_)
    {
      // Compute final time assuming shock Lorentz factor remains
      // fixed at sqrt(2) * lfac_initial
      t1 = p_Obs->t / (p_Obs->z + 1.) * lfac_initial * lfac_initial * 4.;
      #if BOOST_ == ENABLED_
        // Lorentz transform from lab to boosted frame, if needed
        t1 = t1 * eds.boost *
          (1. - eds.beta_sim + eds.beta_sim / 
          (4. * lfac_initial * lfac_initial));
      #endif
    }
    else
    {
      t1 = p_Obs->t / (p_Obs->z + 1.) * 2;
      #if BOOST_ == ENABLED_
        // Lorentz transform from lab to boosted frame, if needed
        t1 = t1 * eds.boost;
      #endif
    }
  }
  else
  {
    t1 = t1_overrule;
  }
  
  if (mbox.box[i_box].BOX_enabled)
  {
    // Since the blast wave decelerates, the estimate will eventually even
    // overshoot the final BOX time, so use that as an additional cap. Also
    // prevent overly large values provided by the user from being applied.
    t1 = fmin(t1, mbox.box[i_box].t[mbox.box[i_box].tres - 1] *
      mbox.box[i_box].Sr );
  }
  else // BM only is used
  {
    BM.E = mbox.E;
    BM.n_ext = mbox.n;
    BM.k = mbox.box[i_box].k;
    BM.set_global_from_lfac(sqrt(2.0) * lfac_final);
    
    #if BOOST_ == ENABLED_
  
      // in boosted frame, times refer to simulation times, not lab time
      taux = eds.boost * (BM.t - eds.beta_sim * BM.R_shock * invv_light);

    #else
  
      taux = BM.t; 
  
    #endif // BOOST_

    t1 = fmin(t1, taux);
  }

  if (t0_overrule < 0.)
  {
    if (!mbox.box[i_box].BM_enabled)
    {
      t0 = mbox.box[i_box].t[0] * mbox.box[i_box].Sr;
    }
    else
    {
      BM.E = mbox.E;
      BM.n_ext = mbox.n;
      BM.k = mbox.box[i_box].k;
      BM.set_global_from_lfac(sqrt(2.0) * lfac_initial);
    
      #if BOOST_ == ENABLED_
  
        // in boosted frame, times refer to simulation times, not lab time
        t0 = eds.boost * (BM.t - eds.beta_sim * BM.R_shock * invv_light);

      #else
  
        t0 = BM.t; 
  
      #endif // BOOST_
    }
  }
  else
  {
    t0 = t0_overrule;
  }

  if (ur_overrule < 0.)
  {
    // Compute EDS size assuming a sphere expanding with fixed Lorentz factor  
    if (mbox.box[i_box].jet == RECEDING_)
    {
      eds.ur_max = v_light * p_Obs->t / (p_Obs->z + 1.);
    }
    else
    {
      eds.ur_max = v_light * p_Obs->t / (p_Obs->z + 1.) * lfac_initial *
        lfac_initial * 4;
    }
  }
  else
  {
    eds.ur_max = ur_overrule;
  }

  // correct if final blast wave size is exceeded
  eds.ur_max = fmin(eds.ur_max, 
    mbox.box[i_box].R_max[0][mbox.box[i_box].tres - 1] *
    mbox.box[i_box].Sr * 1.2);

  
  if (ur0_overrule < 0.)
  {  
    eds.ur_min = eds.ur_max * 1e-6;
  }
  else
  {
    eds.ur_min = ur0_overrule;
  }
  
  eds.reset();
  
  // in case no emission received for this observer time
  if (t0 >= t1)
  {
    F = 0.;
    return 0;
  } 

  //----------------------------------------------------------------------------
  // point directly to individual BOX if theta_0 interpolation not necessary

  emab.p_fluid = &(mbox.box[i_box]);

  //----------------------------------------------------------------------------
  // solve radiative transfer problem

  // reset first and last contributing emission times
  tfirst = 1e300;
  tlast = -1e300;
  Fold = 0.;

  // prepare log file, if applicable
  if (save_emission_profile)
  {
    if (mbox.box[i_box].theta_0 < mbox.th0)
      sprintf(log_filename, "obs%04d-0.log", counter);
    else
      sprintf(log_filename, "obs%04d-1.log", counter);
    
    p_log_file = fopen(log_filename, "w");
    if (p_log_file == NULL)
    {
      printf("ERROR: failure to open log file obs%04d.log by ID %d\n", counter,
        myid);
      fflush(stdout);
      return 1;
    }
    
    fprintf(p_log_file, "##################################################\n");
    fprintf(p_log_file, "# Log file for processor %d. Data point %d\n", myid,
      counter);
    fprintf(p_log_file, "##################################################\n");
    fprintf(p_log_file, "# Observer settings:\n");
    fprintf(p_log_file, "#   r_obs = %e cm\n", p_Obs->dL);
    fprintf(p_log_file, "#   theta_obs = %e rad\n", p_Obs->theta);
    fprintf(p_log_file, "#   z = %e\n", p_Obs->z);
    fprintf(p_log_file, "#   nu_obs = %e Hz\n", p_Obs->nu);
    fprintf(p_log_file, "#   t_obs = %e s (%e days)\n",
      p_Obs->t, p_Obs->t * sec2day);
    fprintf(p_log_file, "# Radiation settings:\n");
    if (emab.electron_cooling)
      fprintf(p_log_file, "#   electron cooling is enabled\n");
    else
      fprintf(p_log_file, "#   electron cooling is disabled\n");
    if(emab.absorption)
      fprintf(p_log_file, "#   self-absorption is enabled\n");
    else
      fprintf(p_log_file, "#   self-absorption is disabled\n");
    if (emab.absorption)
    {
      #if LOCAL_SELF_ABSORPTION_  == DISABLED_
        fprintf(p_log_file, "#   self-absorption computed globally\n");
      #else
        fprintf(p_log_file, "#   self-absorption computed locally\n");
      #endif
    }
    fprintf(p_log_file, "#   synchrotron slope p = %e\n", emab.p_synch);
    fprintf(p_log_file, "#   epsilon_E = %e\n", emab.epsilon_E);
    fprintf(p_log_file, "#   epsilon_B = %e\n", emab.epsilon_B);
    fprintf(p_log_file, "#   ksi_N = %e\n", emab.ksi_N);

    fprintf(p_log_file, "# BOX and resolution settings:\n");
    fprintf(p_log_file, "#   Included emission time range: %e - %e (s)\n",
      t0, t1);
    fprintf(p_log_file, "#   Included radial extent emission image: %e (cm)\n",
      eds.ur_max);
    fprintf(p_log_file, "#   ur_rays = %d\n", eds.ur_rays);
    fprintf(p_log_file, "#   uphi_rays = %d\n", eds.uphi_rays);

    if (mbox.box[i_box].BM_enabled)
    {
      fprintf(p_log_file, "#   Blandford-McKee analytical solution included\n");
      fprintf(p_log_file, "#   Initial BM Lorentz factor: %e\n", lfac_initial );
      if (lfac_final > 0.)
        fprintf(p_log_file, "#   Initial BM Lorentz factor: %e\n", lfac_final);
    }
    
    if (mbox.box[i_box].BOX_enabled)
    {
      fprintf(p_log_file, "#   BOX solution included.\n");
      fprintf(p_log_file, "#   BOX filename: %s\n",
        mbox.box[i_box].BOX_filename);
    }

    if (mbox.box[i_box].jet == FORWARD_)
    {
      fprintf(p_log_file, "#   Computing forward jet only\n");
    }
    if (mbox.box[i_box].jet == BOTH_)
    {
      fprintf(p_log_file, "#   Computing both forward jet and receding jet\n");
    }
    if (mbox.box[i_box].jet == RECEDING_)
    {
      fprintf(p_log_file, "#   Computing receding jet only\n");
    }

    fprintf(p_log_file, "#-------------------------------------------------\n");
    fprintf(p_log_file, "# i, t_e (sim frame), F(t_e), D F, D t_e"
      ", extent ur\n");
  }

  // Global approach to synchrotron self-absorption
  #if LOCAL_SELF_ABSORPTION_ == DISABLED_
  
    t_prev = t0;
    t_current = t0;
    mbox.box[i_box].t_e = t_current;
    mbox.box[i_box].set_global();
    eds.set_coordinates(t_current);
    eds.prepare_update();
  
    for (i_t = 0; i_t < res; i_t++)
    {
      // determine total stepsize
      t_next = fmin(t0 * pow(t1 / t0, (double) (i_t + 1) / (res - 1)), t1);
      Dt_sim = 0.5 * (t_next - t_prev);
    
      // take first substep
      eds.finalize_update(Dt_sim / 6.);
    
      // move to midpoint and prepare next update
      mbox.box[i_box].t_e = t_current + 0.5 * Dt_sim;
      mbox.box[i_box].set_global();
      eds.set_coordinates(mbox.box[i_box].t_e);
      eds.prepare_update();
    
      // take second sub step and low resolution full step
      eds.finalize_update(Dt_sim * 2. / 3.);
      eds.finalize_lores_update(Dt_sim);
    
      // move to endpoint and prepare next update
      mbox.box[i_box].t_e = t_current + Dt_sim;
      mbox.box[i_box].set_global();
      eds.set_coordinates(mbox.box[i_box].t_e);
      eds.prepare_update();
    
      // take last sub step
      eds.finalize_update(Dt_sim / 6.);

      // bookkeeping of first and last contributing emission times
      Ftemp = eds.get_total_flux();
      if (Ftemp != Fold)
      {
        tfirst = fmin(tfirst, t_current);
        tlast = t_current + Dt_sim;
      }
      
      // update log file if appropriate
      if (save_emission_profile)
      {
        eds.set_R();
        fprintf(p_log_file, "%d, %e, %e, %e, %e, %e\n", i_t, 
          mbox.box[i_box].t_e, Ftemp, Ftemp - Fold, Dt_sim, eds.R_100);
      }
      
      Fold = Ftemp;

      // shift times
      t_prev = t_current;
      t_current = t_next;
    }
  
  #endif // LOCAL_SELF_ABSORPTION_ == DISABLED_

  // local approach to synchrotron self-absorption. Warning, very vulnerable to
  // numerical deviation is fluid flow dynamics and only recommended for use
  // with analytical solutions
  
  #if LOCAL_SELF_ABSORPTION_ == ENABLED_
  
    // 4th order Runge-Kutta integrator, acting on F
    t_prev = t0;
    t_current = t0;
    
    mbox.box[i_box].t_e = t_current;
    mbox.box[i_box].set_global();
    eds.set_coordinates(t_current);

    for (i_t = 0; i_t < res; i_t++)
    {
      t_next = fmin(t0 * pow(t1 / t0, (double) (i_t + 1) / (res - 1)), t1);
      Dt_sim = 0.5 * (t_next - t_prev);

      // step 1
      eds.prepare_update(1, Dt_sim);
      
      // step 2
      mbox.box[i_box].t_e = t_current + 0.5 * Dt_sim;
      mbox.box[i_box].set_global();
      eds.set_coordinates(mbox.box[i_box].t_e);
      eds.prepare_update(2, Dt_sim);
      
      // step 3
      eds.prepare_update(3, Dt_sim);
      
      // step 4
      mbox.box[i_box].t_e = t_current + Dt_sim;
      mbox.box[i_box].set_global();
      eds.set_coordinates(mbox.box[i_box].t_e);
      eds.prepare_update(4, Dt_sim);
      
      // combine for full step
      eds.finalize_RK_update();

      // bookkeeping of first and last contributing emission times
      Ftemp = eds.get_total_flux();
      if (Ftemp != Fold)
      {
        tfirst = fmin(tfirst, t_current);
        tlast = t_current;
      }

      // update log file if appropriate
      if (save_emission_profile)
      {
        eds.set_R();
        fprintf(p_log_file, "%d, %e, %e, %e, %e, %e\n", i_t, 
          mbox.box[i_box].t_e, Ftemp, Ftemp - Fold, Dt_sim, eds.R_100);
      }

      Fold = Ftemp;
    
      // shift times
      t_prev = t_current;
      t_current = t_next;
      
    }
  
  #endif // LOCAL_SELF_ABSORPTION_ == ENABLED_
  
  //----------------------------------------------------------------------------
  // Finally, set the flux from the emerging rays

  F_th = eds.get_total_flux();
  eds.set_R();
  
  // close log file, if applicable
  if (save_emission_profile)
  {
    fprintf(p_log_file, "#-------------------------------------------------\n");
    fprintf(p_log_file, "# Earliest contributing emission time estimate: "
      "%e (s)\n", tfirst);
    fprintf(p_log_file, "# Latest contributing emission time estimate: "
      "%e (s)\n", tlast);
    fprintf(p_log_file, "# Radial extent of emission image: %e (cm)\n",
      eds.R_100);
    fprintf(p_log_file, "# extent found / extent set = %e\n",
      eds.R_100 / eds.ur_max);
    fclose(p_log_file);
  }
  
  // now save the eds profile too.
  if (save_image)
  {
    if (mbox.box[i_box].theta_0 < mbox.th0)
      sprintf(image_filename, "image%04d-0.h5", counter);
    else
      sprintf(image_filename, "image%04d-1.h5", counter);
    if (save_image) eds.save_image(image_filename);
  }

  // return success by default
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int c_flux_box :: set_flux()
{
  double F_th0, F_th1; // fluxes of tabulated boxes with surrounding theta_0's
  int status; // zero if function call successfull
  
  // reset flux
  F = 0.; F_th0 = 0.; F_th1 = 0.;
  
  // set surrounding box flux values
  if (mbox.fractheta0 < 1. - 1e-3)
  {
    status = set_flux_th(mbox.i_th0, F_th0);
    if (status != 0) return status;
  }
  
  if (mbox.fractheta0 > 1e-3)
  {
    status = set_flux_th(mbox.i_th0 + 1, F_th1);
    if (status != 0) return status;
  }

  if (mbox.fractheta0 <= 1e-3)
  {
    F = F_th0;
    return 0;
  }
  
  if (mbox.fractheta0 >= 1. - 1e-3)
  {
    F = F_th1;
    return 0;
  }
  
  // combine fluxes
  if (F_th0 <= 0. or F_th1 <= 0.)
  {
    F = 0.;
    return 0;
  }

  // weigh the surrounding theta_0 fluxes in log space
  F = exp(mbox.fractheta0 * log(F_th1) + (1. - mbox.fractheta0) * log(F_th0));

  // return success by default
  return 0;
}
