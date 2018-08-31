////////////////////////////////////////////////////////////////////////////////
//
// fluid_BM.cpp
//
// Created June 29, 2011, by HJvE
// Last Modified Nov 24, 2014 by HJvE
//
// Contains the analytical grid
// Based on Blandford & McKee 1976, Phys. Fluid 19, 1130
//
////////////////////////////////////////////////////////////////////////////////

#include "fluid_special.h"

////////////////////////////////////////////////////////////////////////////////

c_fluid_special :: c_fluid_special()
{
  #if BOOST_ == ENABLED_
    boost = 1.0;
    boostsqr = 1.0;
    beta_sim = 0.0;
  #endif // BOOST_
}

////////////////////////////////////////////////////////////////////////////////

void c_fluid_special :: init_times(double toverz)
{
  // set up timespan and global settings

  if (toverz < 0. or lfac_initial > 0.) 
  {
    // determine timespan for user-provided lfac's
  
    // final emission time for analytical coverage, but first emission time
    // for snapshot coverage
    BM.set_global_from_lfac(sqrt(2.0) * lfac_final);
 
    #if BOOST_ == ENABLED_
  
      // in boosted frame, times refer to simulation times, not lab time
      t_e_final = boost * (BM.t - beta_sim * BM.R_shock * invv_light);

    #else
  
      t_e_final = BM.t; 
  
    #endif // BOOST_

    //--------------------------------------------------------------------------

    // initial emission time for analytical coverage
    BM.set_global_from_lfac(sqrt(2.0) * lfac_initial);

    #if BOOST_ == ENABLED_
    
      // in boosted frame, times refer to simulation times, not lab time
      t_e_initial = boost * (BM.t - beta_sim * BM.R_shock * invv_light);

    #else

      t_e_initial = BM.t;

    #endif // BOOST_
  }
  else // FOR NOW, same thing even when times are provided. TO DO ADD PROPER TWEAK HERE!
  {

    // final emission time for analytical coverage, but first emission time
    // for snapshot coverage
    BM.set_global_from_lfac(sqrt(2.0) * lfac_final);
 
    #if BOOST_ == ENABLED_
  
      // in boosted frame, times refer to simulation times, not lab time
      t_e_final = boost * (BM.t - beta_sim * BM.R_shock * invv_light);

    #else
  
      t_e_final = BM.t; 
  
    #endif // BOOST_

    //--------------------------------------------------------------------------

    // initial emission time for analytical coverage
    BM.set_global_from_lfac(sqrt(2.0) * lfac_initial);

    #if BOOST_ == ENABLED_
    
      // in boosted frame, times refer to simulation times, not lab time
      t_e_initial = boost * (BM.t - beta_sim * BM.R_shock * invv_light);

    #else

      t_e_initial = BM.t;

    #endif // BOOST_

  }
}

////////////////////////////////////////////////////////////////////////////////

void c_fluid_special :: set_global()
{
  BM.t = t_e;
  BM.set_global_from_time(BM.t);
}

////////////////////////////////////////////////////////////////////////////////

void c_fluid_special :: set_local(s_coordinates cor)
{
  //----------------------------------------------------------------------------
  // if the frame is boosted, lab frame time differs across grid, and we can no
  // longer assume a global state to have been initialized in a way that
  // properly reflects local conditions
  #if BOOST_ == ENABLED_
  
    BM.t = cor.t;
    BM.set_global_from_time(cor.t);

  #endif

  // if outside of opening angle, use a value in front of shock
  if (cor.theta > theta_max and theta_max > 0.0) BM.set_local(2.0 * BM.R_shock);
    else BM.set_local(cor.r);

  local[lfac_] = BM.lfac;
  local[rho_] = BM.rho;
  local[eint_] = BM.e_therm;  

  local[v_x_] = BM.beta;
  local[v_y_] = 0.0;
  local[v_] = local[v_x_];
  local[N_] = local[rho_] * invm_p;
  local[D_] = local[lfac_] * local[rho_];
}

////////////////////////////////////////////////////////////////////////////////

void c_fluid_special :: load_settings(char *parfilename, int argc, char* argv[])
{
  BM.k = double_from_parfile(parfilename, "special_k");
  if (parse("-special_k=", argc, argv)) // command line overrules
    parse_double("-special_k=", BM.k, argc, argv);
  
  if (BM.k > 0.0)
  {
    BM.R_scale = double_from_parfile(parfilename, "special_R_ref");
    if (parse("-special_R_ref=", argc, argv)) // command line overrules
      parse_double("-special_R_ref=", BM.R_scale, argc, argv);
  }

  BM.E = double_from_parfile(parfilename, "special_E_iso");
  if (parse("-special_E_iso=", argc, argv)) // command line overrules
    parse_double("-special_E_iso=", BM.E, argc, argv);

  BM.n_ext = double_from_parfile(parfilename, "special_n_ref");
  if (parse("-special_n_ref=", argc, argv)) // command line overrules
    parse_double("-special_n_ref=", BM.n_ext, argc, argv);

  lfac_final = double_from_parfile(parfilename, "special_lfac_final");
  if (parse("-special_lfac_final=", argc, argv)) // command line overrules
    parse_double("-special_lfac_final=", lfac_final, argc, argv);

  lfac_initial = double_from_parfile(parfilename, "special_lfac_initial");
  if (parse("-special_lfac_initial=", argc, argv)) // command line overrules
    parse_double("-special_lfac_initial=", lfac_initial, argc, argv);

  theta_max = double_from_parfile(parfilename, "special_theta_max");
  if (parse("-special_theta_max=", argc, argv)) // command line overrules
    parse_double("-special_theta_max=", theta_max, argc, argv);

  // call set times to be able to report effect of lfac_final and lfac_initial
  init_times(-1.);

  //----------------------------------------------------------------------------
  // error-check the settings
  if (myid == 0)
  {
    if (lfac_final < 1.)
    {
      printf("ERROR reading user settings BM solution. lfac_final < 1\n");
      fflush(stdout);
      exit(1);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void c_fluid_special :: print_log(FILE *p_logfile)
{
  fprintf(p_logfile, "# Special fluid is Blandford-McKee solution\n");
  fprintf(p_logfile, "#   Simulation frame times: %e - %e (s)\n", t_e_initial,
    t_e_final);
  fprintf(p_logfile, "#       k: %e\n", BM.k);
  if (BM.k > 0.0)
    fprintf(p_logfile, "#   R_ref: %e (cm)\n", BM.R_scale);
  fprintf(p_logfile, "#   E_iso: %e (erg)\n", BM.E);
  fprintf(p_logfile, "#   n_ext: %e (cm^-3)\n", BM.n_ext);
  fprintf(p_logfile, "#   theta_max: %e (rad)\n", theta_max);
  fprintf(p_logfile, "#   Initial Lorentz factor: %e\n", lfac_initial);
  fprintf(p_logfile, "#   final Lorentz factor: %e\n", lfac_final);

  #if BOOST_ == ENABLED_
  
    fprintf(p_logfile, "#   Boosted frame Lorentz factor: %e\n", boost);
    BM.set_global_from_lfac(sqrt(2.0) * lfac_final);
    fprintf(p_logfile, "#   Final tip radius (lab): %e\n",
      BM.R_shock);

  #endif // BOOST_
}

////////////////////////////////////////////////////////////////////////////////

double c_fluid_special :: get_var(int vartype, s_coordinates cor)
{
  // if outside of opening angle, use a value in front of shock
  if (cor.theta > theta_max) BM.set_local(2.0 * BM.R_shock);
    else BM.set_local(cor.r);

  switch( vartype )
  {
    case eint_: return BM.e_therm;
    case rho_: return BM.rho;
    case tau_: return BM.tau;
    case lfac_: return BM.lfac;
    case v_: return BM.beta;
    case v_x_: return BM.beta;
    case v_y_: return 0.0;
    default: 
      printf("c_fluid_BM :: get_var ERROR unknown var.\n"); 
      fflush(stdout);
      abort();
  }

  return -1.0;
}
