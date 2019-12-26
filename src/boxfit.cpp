////////////////////////////////////////////////////////////////////////////////
//
// boxfit.cpp
//
// Created, Jun 28 2011, by HJvE
//
//
// Reference: van Eerten, van der Horst & MacFadyen 2012, ApJ 749, 44
//
////////////////////////////////////////////////////////////////////////////////

#include "boxfit.h"

////////////////////////////////////////////////////////////////////////////////
// global variables

int numprocs; // stores the number of processors used for this run
int myid; // identity of single processor, 0 for non-parallel run
c_boxfit *p_boxfit; // pointer to boxfit class for external fit function wrapper
int noderank; // identity of single core within node

#if OPEN_MPI_ == ENABLED_
  int nodesize;
  MPI_Comm nodecom;
#endif

////////////////////////////////////////////////////////////////////////////////

double funk_wrapper(double fv[]) // external wrapper for fit function
{
  return p_boxfit->funk(fv);
}

////////////////////////////////////////////////////////////////////////////////

double selftest_funk_wrapper(double fv[])
{
  return p_boxfit->bkn_power(fv);
}

////////////////////////////////////////////////////////////////////////////////

c_dataset :: c_dataset()
{
  initialized = false;
}

////////////////////////////////////////////////////////////////////////////////     

c_dataset :: ~c_dataset()
{
  if (initialized)
  {
    free_memory();
  }
}

////////////////////////////////////////////////////////////////////////////////     

void c_dataset :: free_memory()
{
  if (initialized)
  {
    delete_array_1d(F);
    delete_array_1d(Falt);    
    delete_array_1d(Fbox);    
    delete_array_1d(t);    
    delete_array_1d(nu);    
    delete_array_1d(sigma);    
    
    delete_array_1d(logF);
    delete_array_2d(dlogFdx);
    
    initialized = false;
  }
}

////////////////////////////////////////////////////////////////////////////////     

int c_dataset :: print(FILE *p_file)
{
  int chk;
  int i;

  // print number of datapoints
  chk = fprintf(p_file, "# datapoints = %d\n", n);
  if (chk < 0)
  {
    printf("c_dataset :: save ERROR: could not save no of datapoints\n");
    fflush(stdout);
    return 1;
  }

  // print data
  fprintf(p_file, "# t (s), nu (Hz), data F (mJy), data dF (mJy), box F "
    "(mJy), alt. F (mJy)\n");
  for (i = 0; i < n; i++)
  {
    chk = fprintf(p_file, "%1.4le, %1.4le, %1.4le, %1.4le, %1.4le, %1.4le\n", 
      t[i], nu[i], F[i], sigma[i], Fbox[i], Falt[i]);
    if (chk < 0)
    {
      printf("c_dataset :: save ERROR. could not save datapoint %d\n", i);
      fflush(stdout);
      return 1;
    }
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////////////     

int c_dataset :: load(const char* filename)
{
  FILE * p_file;
  int chk;
  int i;

  p_file = fopen(filename, "r");
  
  if (p_file == NULL)
  {
    if (myid == 0) printf("c_dataset :: load ERROR. Could not open %s\n",
      filename);
    fflush(stdout);
    return 1;
  }

  // read header
  chk = fscanf(p_file, "# datapoints = %d\n", &n);
  if (chk != 1)
  {
    printf("c_dataset :: load ERROR: could not read no of datapoints\n");
    fflush(stdout);
    fclose(p_file);
    return 1;
  }

  // initialize memory
  initialize();

  // read data
  for (i = 0; i < n; i++)
  {
    chk = fscanf(p_file, "%le, %le, %le, %le\n", &t[i], &nu[i], &F[i], 
      &sigma[i]);
    Falt[i] = F[i]; // set default values for alternative fluxes.
    if (chk != 4)
    {
      printf("c_dataset :: load ERROR. could not read datapoint %d\n", i);
      fflush(stdout);
      fclose(p_file);
      return 1;
    }
    t[i] = t[i];
  }

  // close file 
  fclose(p_file);
  
  return 0;
}

////////////////////////////////////////////////////////////////////////////////     

void c_dataset :: initialize()
{
  if (!initialized)
  {
    F = array_1d<double>(n);
    Falt = array_1d<double>(n);
    Fbox = array_1d<double>(n);
    sigma = array_1d<double>(n);
    t = array_1d<double>(n);
    nu = array_1d<double>(n);
    
    // Monte Carlo error determination related
    logF = array_1d<double>(n);
    dlogFdx = array_2d<double>(n, no_fit_vars_);
    
    initialized = true;
  }
  else
  {
    printf("c_dataset ERROR: data already initialized\n"); 
    fflush(stdout); 
    abort();
  }
}

////////////////////////////////////////////////////////////////////////////////     

c_boxfit :: c_boxfit()
{
  no_boxes = 0;

  fitvar_allocated = false;
  varmax_allocated = false;
  varmin_allocated = false;
  subvarmin_allocated = false;
  subvarmax_allocated = false;
  varfrozen_allocated = false;
  fitsubentry_allocated = false;
  fitvar_initial_allocated = false;
  initial_simplex_allocated = false;
  simplex_max_allocated = false;
  
  data_filename_allocated = false;
  box_filename_base_allocated = false;
  simplex_filename_allocated = false;
  parfilename_allocated = false;
  
  box_filename_allocated = false;

  temp_string = new char[1]; // dummy allocation
}

////////////////////////////////////////////////////////////////////////////////     

c_boxfit :: ~c_boxfit()
{
  release_memory(); // call release_memory once more, to be sure no memory leaks
    // remain on unexpected shutdown.

  delete[] temp_string; // always allocated (dummy allocation in constructor)
}

////////////////////////////////////////////////////////////////////////////////     

int c_boxfit :: main(int argc, char* argv[])
{
  double x;
  FILE *p_file;
  int i; // dummy loop index
  bool parfilegiven = false; // assume no parameterfile given at on command line
  int success; // 0 if OK
  
  //----------------------------------------------------------------------------
  // print opening message on screen
  
  if (myid == 0) 
  {
    print_opening_message(stdout);
    printf("# function call: ");
    for (i = 0; i < argc; i++)
    {
      printf(" %s", argv[i]);
    }
    printf("\n");
  }

  //----------------------------------------------------------------------------
  // Set up user-provided parameter file IO

  // check if parameter file given
  for (i = 1; i < argc; i++)
    if ((argv[i])[0] != '-') parfilegiven = true;
    
  if (!parfilegiven) // if not, assume one
  {
    parfilename = new char[19];
    sprintf(parfilename, "boxfitsettings.txt");
  }
  else // if given, store it
  {
    parfilename = new char[strlen(argv[1])+1];
    sprintf(parfilename, "%s", argv[1]);
  }
  parfilename_allocated = true;

  if (!file_exists(parfilename))
  {
    if (myid == 0)
    {
      printf("ERROR: user parameter file %s not found.\n", parfilename);
      fflush(stdout);
    }
    return 1;
  }

  #if OPEN_MPI_ == ENABLED_
    MPI_Barrier(MPI_COMM_WORLD); // wait until all cores have caught up
  #endif

  //----------------------------------------------------------------------------

  success = initialize(argc, argv);
  if (success > 0) return success;
  
  //----------------------------------------------------------------------------
  // select action to perform
   
  switch(what_to_do)
  {
    case -1: // self test
      success = selftest();
      if (success != 0) return success;
    break;
    
    case 0: // calculate chi^2 once
      calculate(); // calculate the box flux at the datapoints
      x = set_chisqr(); // calculate chi^2 based on most recent box calculation
      if (myid == 0)
      { 
        // print output on screen
        print_settings(stdout);
        printf("Chi^2 = %e\n", x);
        printf("Reduced Chi^2 = %e\n", x / (double) (dataset.n - thawed));
        fflush(stdout);
        
        // write output to file
        p_file = fopen("single.txt", "w");
        print_opening_message(p_file);
        print_settings(p_file);
        fprintf(p_file, "# Chi^2 = %e\n", x);
        fprintf(p_file, "# Reduced Chi^2 = %e\n", x / (double) 
          (dataset.n - thawed));
        success = dataset.print(p_file);
        fclose(p_file);
        if (!success) return 1;
      }
    break;
    
    case 1: // calculate light curve

      fluxes = array_1d<double>(no_points);
      if (myid == 0) print_settings(stdout);
      lightcurve();

      if (myid == 0)
      {
        // write output to screen
        print_fluxes(stdout);
        
        // write output to file
        p_file = fopen("lightcurve.txt", "w");
        print_opening_message(p_file);
        print_settings(p_file);
        print_fluxes(p_file);
        fclose(p_file);
      }
      
      delete_array_1d(fluxes);
    break;
    
    case 2: // calculate spectrum
      fluxes = array_1d<double>(no_points);
    
      spectrum();

      if (myid == 0)
      {
        // write output to screen
        print_settings(stdout);
        print_fluxes(stdout);
        
        // write output to file
        p_file = fopen("spectrum.txt", "w");
        print_opening_message(p_file);
        print_settings(p_file);
        print_fluxes(p_file);
        fclose(p_file);
      }
      
      delete_array_1d(fluxes);
    break;
    
    case 3: // perform iterative fit to data
      initialize_simplex();
     
      if (myid == 0)
      {
        print_settings(stdout);
        fprintf(stdout, "i, temp, theta_0, E, n, theta_obs, p, epsilon_B, "
          "epsilon_E, ksi_N, chi^2, red. chi^2\n");
      }
      success = fit_dataset();
    break;
    
    case 4: // Monte Carlo error determination based on full box calculations
      initialize_simplex();

      if (myid == 0) print_settings(stdout);
      success = montecarlo();
    break;
    
    case 5: // determine partial derivatives and exit
      initialize_simplex();

      if (myid == 0) print_settings(stdout);
      if (setup_derivatives() != 0) return 1; 
    break;
    
    case 6: // Monte Carlo error determination based on partial derivatives
      initialize_simplex();

      if (myid == 0) print_settings(stdout);
      load_derivatives();
      fromderivatives = true;
      success = montecarlo();
    break;
  }
  
  //----------------------------------------------------------------------------
  
  // release memory, close files etc.
  release_memory();

  if (myid == 0) { printf("# Program completed.\n"); fflush(stdout); }

  return success;
}

////////////////////////////////////////////////////////////////////////////////

double c_boxfit :: bkn_power(double fv[])
{
  int i;
  
  // update fitvars
  for (i = 0; i < thawed; i++)
    fitvar[fitsubentry[i]] = fv[i];
  
  // loop over dataset
  for (i = 0; i < dataset.n; i++)
  {
    if (dataset.t[i] * sec2day < fitvar[fit_t_break_])
      dataset.Fbox[i] = fitvar[fit_S_] * pow(dataset.t[i] * sec2day / 
        fitvar[fit_t_break_], -fitvar[fit_alpha_1_]);
    else
      dataset.Fbox[i] = fitvar[fit_S_] * pow(dataset.t[i] * sec2day / 
        fitvar[fit_t_break_], -fitvar[fit_alpha_2_]);
  }
  
  return set_chisqr();
}

////////////////////////////////////////////////////////////////////////////////

double c_boxfit :: funk(double fv[])
{
  int i;
  
  // update fit variables
  for (i = 0; i < thawed; i++)
    fitvar[fitsubentry[i]] = fv[i];
  
  // recalculate light curve
  calculate();
  
  // return chi^2
  return set_chisqr();
}

////////////////////////////////////////////////////////////////////////////////     

int c_boxfit :: selftest()
{
  int i;
  int success = 0; // track error flags
  
  // fflush(stdout); return;
  
  // going to do a broken power law fit to the data
  thawed = 4.0;
  for (i = 0; i < 4; i++)
    fitsubentry[i] = i;
  
  fitvar[fit_S_] = 1e-2;
  varmin[fit_S_] = 1e-5;
  varmax[fit_S_] = 1e3;
  
  fitvar[fit_alpha_1_] = 1.0;
  varmin[fit_alpha_1_] = 1e-2;
  varmax[fit_alpha_1_] = 10.0;
  
  fitvar[fit_t_break_] = 2.0;
  varmin[fit_t_break_] = 1e-3;
  varmax[fit_t_break_] = 200;
    
  fitvar[fit_alpha_2_] = 4.0;
  varmin[fit_alpha_2_] = 1e-2;
  varmax[fit_alpha_2_] = 50;
  
  // set up the dataset
  for (i = 0; i < dataset.n; i++)
  {
    if (dataset.t[i] * sec2day < 3.0)
      dataset.F[i] = fitvar[fit_S_] * pow(dataset.t[i] * sec2day / 
        3.0, -fitvar[fit_alpha_1_]);
    else
      dataset.F[i] = fitvar[fit_S_] * pow(dataset.t[i] * sec2day / 
        3.0, -fitvar[fit_alpha_2_]);
  } 
  
  // freeze all, unfreeze broken power law fit parameters
  for (i = 0; i < no_fit_vars_; i++) varfrozen[i] = true;
  varfrozen[fit_S_] = false;
  varfrozen[fit_alpha_1_] = false;
  varfrozen[fit_t_break_] = false;
  varfrozen[fit_alpha_2_] = false;

  simplex.clear_memory(); 
  simplex.initialize(4);

  for (i = 0; i < thawed; i++)
  {
    subvarmin[i] = varmin[fitsubentry[i]];
    subvarmax[i] = varmax[fitsubentry[i]];
    simplex.x_min[i] = varmin[fitsubentry[i]];
    simplex.x_max[i] = varmax[fitsubentry[i]];
  }
  
  initialize_simplex();
  simplex.p_f = selftest_funk_wrapper;
  success = fit_dataset();
  if (success != 0) return success;

  // write output to file
  FILE *p_file;
  p_file = fopen("lightcurve.txt", "w");
  print_opening_message(p_file);
  print_settings(p_file);
  success = dataset.print(p_file);
  fclose(p_file);
    
  return success;
}

////////////////////////////////////////////////////////////////////////////////

void c_boxfit :: print_opening_message(FILE *p_where)
{
  fprintf(p_where, "#######################################################\n");
  fprintf(p_where, "#              BOX-FIT. Release 2 \n");
  fprintf(p_where, "# Last modified: July 29, 2016, by HJvE \n");
  fprintf(p_where, "#\n");
  fprintf(p_where, "# reference: \"Gamma-ray burst afterglow broadband \n");
  fprintf(p_where, "#   fitting based directly on hydrodynamics simulations\""
    "\n");
  fprintf(p_where, "#   H.J. van Eerten, A.J. van der Horst, A.I. MacFadyen\n");
  fprintf(p_where, "#   ApJ (2012) Issue 749, Page 44\n");
  fprintf(p_where, "#   ArXiv: 1110.5089\n");
  fprintf(p_where, "#\n");
  fprintf(p_where, "# Development of the Boxfit code was supported in part \n");
  fprintf(p_where, "# by NASA through grant NNX10AF62G issued through the \n");
  fprintf(p_where, "# Astrophysics Theory Program and by the NSF through \n");
  fprintf(p_where, "# grant AST-1009863.\n");
  fprintf(p_where, "#######################################################\n");
}

////////////////////////////////////////////////////////////////////////////////

void c_boxfit :: print_fit_result(FILE *p_where)
{
  fprintf(p_where, "%05d, %1.3e, %1.3e, %1.3e, %1.3e, %1.3e, %1.3e, %1.3e, "
    "%1.3e, %1.3e, %1.3e, %1.3e\n", lccounter, simplex.temp, 
    fitvar[fit_theta_0_], pow(10.0, fitvar[fit_E_]), 
    pow(10.0, fitvar[fit_n_]), fitvar[fit_theta_obs_], fitvar[fit_p_], 
    pow(10.0, fitvar[fit_epsilon_B_]), pow(10.0, fitvar[fit_epsilon_E_]), 
    pow(10.0, fitvar[fit_ksi_N_]), chi_sqr, 
    chi_sqr / (double) (dataset.n - thawed));
  fflush(p_where);
}

////////////////////////////////////////////////////////////////////////////////

void c_boxfit :: print_fluxes(FILE *p_where)
{
  int i;
  double t, nu;
  
  fprintf(p_where, "# i, t (s), nu (Hz), F (mJy)\n");
  
  switch(what_to_do)
  {
    case 1: // light curve
      for (i = 0; i < no_points; i++)
      {
        t = t_0 * pow(t_1 / t_0, (double) i / (double) (no_points - 1));
        fprintf(p_where, "%d, %e, %e, %e\n", i, t, nu_0, fluxes[i]);
      }
    break;
    case 2: // spectrum
      for (i = 0; i < no_points; i++)
      {
        nu = nu_0 * pow(nu_1 / nu_0, (double) i / (double) (no_points - 1));
        fprintf(p_where, "%d, %e, %e, %e\n", i, t_0, nu, fluxes[i]);
      }
  }
}

////////////////////////////////////////////////////////////////////////////////

void c_boxfit :: print_settings(FILE *p_where)
{
  int i; // dummy loop index
  
  if (myid == 0)
  {
    #if BOOST_ == ENABLED_
    
      fprintf(p_where, "# BOXFIT Compiled WITH boosted frames support\n");
    
    #else

      fprintf(p_where, "# BOXFIT Compiled WITHOUT boosted frames support\n");

    #endif
    fprintf(p_where, "# I/O settings:\n");
    fprintf(p_where, "#   using parameter file %s\n", parfilename);
    fprintf(p_where, "#   base filename BOX files: %s\n", box_filename_base);

    if (what_to_do != -1 and what_to_do != 1 and what_to_do != 2)
      fprintf(p_where, "#   dataset file name: %s\n", data_filename);
    
    if (what_to_do == 1 or what_to_do == 2)
    {
      if (flux.save_emission_profile)
        fprintf(p_where, "#   Storing intermediate d F / d t_e\n");
      else
        fprintf(p_where, "#   Not storing intermediate d F / d t_e\n");
    }
    
    fprintf(p_where, "#####################################################\n");
    fprintf(p_where, "# radiation switches\n");

    if (usessa)
    {
      fprintf(p_where, "#   synchrotron self-absorption is enabled\n");
      #if LOCAL_SELF_ABSORPTION_ == DISABLED_
        fprintf(p_where, "#   synchrotron self-absorption is computed "
          "globally for each ray.\n");
      #else
        fprintf(p_where, "#   synchrotron self-absorption is computed "
          "locally for each ray.\n");
      #endif
    }
    else fprintf(p_where, "#   synchrotron self-absorption is disabled\n");

    if (usecooling) fprintf(p_where, "#   electron cooling is enabled\n");
      else fprintf(p_where, "#   electron cooling is disabled\n");
    if (usebox) fprintf(p_where, "#   box data is included\n");
      else fprintf(p_where, "#   box data is excluded\n");
    if (useBM) 
      fprintf(p_where, "#   Blandford-McKee analytical data is included\n");
      else 
        fprintf(p_where, "#   Blandford-McKee analytical data is excluded\n");

    // jet settings
    if (flux.mbox.box[0].jet == FORWARD_)
      fprintf(p_where, "#   Calculating forward jet only.\n");
    else if (flux.mbox.box[0].jet == RECEDING_)
      fprintf(p_where, "#   Calculating receding jet only.\n");
    else
      fprintf(p_where, "#   Calculating both forward and receding jets.\n");

    fprintf(p_where, "#####################################################\n");
    fprintf(p_where, "# Model parameter settings:\n");
    
    fprintf(p_where, "#   theta_0 = %e (rad)", fitvar_initial[fit_theta_0_]);
    if (what_to_do == 3 or what_to_do == 4 or what_to_do == 5 or 
      what_to_do == 6)
    {
      fprintf(p_where, ", between %e and %e.", varmin[fit_theta_0_], 
        varmax[fit_theta_0_]);
      if (varfrozen[fit_theta_0_]) fprintf(p_where, " Frozen.\n");
        else fprintf(p_where, " Thawed.\n");
    }
    else
      fprintf(p_where, "\n");
    
    fprintf(p_where, "#   E_iso = %e (erg)", pow(10.0, 
      fitvar_initial[fit_E_]));
    if (what_to_do == 3 or what_to_do == 4 or what_to_do == 5 or 
      what_to_do == 6)
    {
      fprintf(p_where, ", between %e and %e.", pow(10.0, varmin[fit_E_]), 
        pow(10.0, varmax[fit_E_]));
      if (varfrozen[fit_E_]) fprintf(p_where, " Frozen.\n");
        else fprintf(p_where, " Thawed.\n");
    }
    else
      fprintf(p_where, "\n");
    
    fprintf(p_where, "#   n_0 = %e (cm^-3)", pow(10.0, fitvar_initial[fit_n_]));
    if (what_to_do == 3 or what_to_do == 4 or what_to_do == 5 or 
      what_to_do == 6)
    {
      fprintf(p_where, ", between %e and %e.", pow(10.0, varmin[fit_n_]), 
        pow(10.0, varmax[fit_n_]));
      if (varfrozen[fit_n_]) fprintf(p_where, " Frozen.\n");
        else fprintf(p_where, " Thawed.\n");
    }
    else
      fprintf(p_where, "\n");

    fprintf(p_where, "#   theta_obs = %e (rad)", 
      fitvar_initial[fit_theta_obs_]);
    if (what_to_do == 3 or what_to_do == 4 or what_to_do == 5 or 
      what_to_do == 6)
    {
      fprintf(p_where, ", between %e and %e.", varmin[fit_theta_obs_], 
        varmax[fit_theta_obs_]);
      if (varfrozen[fit_theta_obs_]) fprintf(p_where, " Frozen.\n");
        else fprintf(p_where, " Thawed.\n");
    }
    else
      fprintf(p_where, "\n");
    
    fprintf(p_where, "#   p = %e", fitvar_initial[fit_p_]);
    if (what_to_do == 3 or what_to_do == 4 or what_to_do == 5 or 
      what_to_do == 6)
    {
      fprintf(p_where, ", between %e and %e.", varmin[fit_p_], varmax[fit_p_]);
      if (varfrozen[fit_p_]) fprintf(p_where, " Frozen.\n");
        else fprintf(p_where, " Thawed.\n");
    }
    else
      fprintf(p_where, "\n");
    
    fprintf(p_where, "#   epsilon_B = %e", pow(10., 
      fitvar_initial[fit_epsilon_B_]));
    if (what_to_do == 3 or what_to_do == 4 or what_to_do == 5 or 
      what_to_do == 6)
    {
      fprintf(p_where, ", between %e and %e.", pow(10., varmin[fit_epsilon_B_]), 
        pow(10., varmax[fit_epsilon_B_]));
      if (varfrozen[fit_epsilon_B_]) fprintf(p_where, " Frozen.\n");
        else fprintf(p_where, " Thawed.\n");
    }
    else
      fprintf(p_where, "\n");
    
    fprintf(p_where, "#   epsilon_E = %e", pow(10., 
      fitvar_initial[fit_epsilon_E_]));
    if (what_to_do == 3 or what_to_do == 4 or what_to_do == 5 or 
      what_to_do == 6)
    {
      fprintf(p_where, ", between %e and %e.", pow(10., varmin[fit_epsilon_E_]), 
        pow(10., varmax[fit_epsilon_E_]));
      if (varfrozen[fit_epsilon_E_]) fprintf(p_where, " Frozen.\n");
        else fprintf(p_where, " Thawed.\n");
    }
    else
      fprintf(p_where, "\n");
    
    fprintf(p_where, "#   ksi_N = %e", pow(10., fitvar_initial[fit_ksi_N_]));
    if (what_to_do == 3 or what_to_do == 4 or what_to_do == 5 or 
      what_to_do == 6)
    {
      fprintf(p_where, ", between %e and %e.", pow(10.0, varmin[fit_ksi_N_]), 
        pow(10.0, varmax[fit_ksi_N_]));
      if (varfrozen[fit_ksi_N_]) fprintf(p_where, " Frozen.\n");
        else fprintf(p_where, " Thawed.\n");
    }
    else
      fprintf(p_where, "\n");

    // only print for light curves or spectra
    if (what_to_do == 1 or what_to_do == 2)
    {
      fprintf(p_where, "###################################################\n");
      fprintf(p_where, "# Time and frequency settings\n");

      if (what_to_do == 2)
        fprintf(p_where, "#   nu_0 = %e (Hz), nu_1 = %e (Hz)\n", nu_0, nu_1);
      else
        fprintf(p_where, "#   nu_0 = %e (Hz)\n", nu_0);

      if (what_to_do == 1)
        fprintf(p_where, "#   t_0 = %e (s), t_1 = %e (s)\n", 
          t_0, t_1);
      else
        fprintf(p_where, "#   t_0 = %e (s)\n", t_0);
    }

    fprintf(p_where, "#####################################################\n");
    fprintf(p_where, "# Observer distance settings\n");
    fprintf(p_where, "#   luminosity distance: %e (cm)\n", Obs.dL);
    fprintf(p_where, "#   redshift: %e\n", Obs.z);

    fprintf(p_where, "#####################################################\n");
    fprintf(p_where, "# BOX data:\n");
    fprintf(p_where, "#   number of boxes: %d\n", no_boxes);
    for (i = 0; i < no_boxes; i++)
      fprintf(p_where, "#   box %d : theta_0 = %e\n", i, 
        flux.mbox.box[i].theta_0);

    #if BOOST_ == ENABLED_

      fprintf(p_where, "#   Box boost Lorentz factor: %e\n",
        flux.mbox.box[0].boost);

      if (flux.mbox.box[0].counterjet)
        fprintf(p_where, "#   Counter-jet data is available in the BOX files"
          ".\n");
      else
        fprintf(p_where, "#   Counter-jet data is not available in the BOX"
          " files.\n");

    #endif
    
    if (flux.mbox.box[0].k == 0)
      fprintf(p_where, "#   BOX data for k = 0 (ISM environment).\n");
    else
      fprintf(p_where, "#   BOX data for k = 2 (stellar-wind environment).\n");
    
    fprintf(p_where, "#####################################################\n");
    fprintf(p_where, "# Resolution settings:\n");
    fprintf(p_where, "#   BM / BOX emission time resolution 'res' = %d\n", 
      flux.res);
    fprintf(p_where, "#   EDS resolution: phi rays = %d, r rays = %d\n", 
      flux.eds.uphi_rays, flux.eds.ur_rays);

    if (flux.ur_overrule > 0.)
      fprintf(p_where, "#   dynamic EDS ur_max overruled and set to %e\n",
        flux.ur_overrule);
    else
      fprintf(p_where, "#   EDS ur_max computed dynamically\n");
    
    if (flux.t0_overrule > 0.)
      fprintf(p_where, "#   dynamic emission time t_0 overruled and set to %e\n"
        ,flux.t0_overrule);
    else
      fprintf(p_where, "#   Emission time t_0 computed dynamically\n");

    if (flux.t1_overrule > 0.)
      fprintf(p_where, "#   dynamic emission time t_1 overruled and set to %e\n"
        ,flux.t1_overrule);
    else
      fprintf(p_where, "#   Emission time t_1 computed dynamically\n");

    if (useBM)
    {
      fprintf(p_where, 
        "#   Starting Blandford-McKee fluid Lorentz factor: %e\n",
        flux.lfac_initial);
      if (flux.lfac_final > 1)
        fprintf(p_where, 
          "#   Stopping Blandford-McKee fluid Lorentz factor: %e\n",
          flux.lfac_final);
    }

    if (flux.ur_overrule > 0.)
    {
      fprintf(p_where, "#   EDS outer radius manually set to %e cm\n",
        flux.ur_overrule);
    }
    else
    {
      fprintf(p_where, "#   EDS outer radius estimated analytically\n");
    }
    
    if (flux.t0_overrule > 0.)
    {
      fprintf(p_where, "#   Emission frame time t0 manually set to %e\n",
        flux.t0_overrule);
    }
    else
    {
      fprintf(p_where, "#   Emission frame time t0 estimated analytically.\n");
    }
    
    if (flux.t1_overrule > 0.)
    {
      fprintf(p_where, "#   Emission frame time t1 manually set to %e\n",
        flux.t1_overrule);
    }
    else
    {
      fprintf(p_where, "#   Emission frame time t1 estimated analytically.\n");
    }

    if (what_to_do != -1 and what_to_do != 1 and what_to_do != 2)
    {
      fprintf(p_where, "#   Starting temperature annealing: %e\n", temp);
      fprintf(p_where, "#   Number of iterations at single temperature: %d\n",
        iter);
      fprintf(p_where, "#   Fraction drop in temperature per cycle: %e\n",
        temp_factor);
      fprintf(p_where, "#   Lowest nonzero annealing temperature: %e\n",
        temp_lowest);
      fprintf(p_where, "#   Number of Monte Carlo light curves calculated: %d\n"
        , MC_runs);
    }
      
    fprintf(p_where, "#####################################################\n");
    fflush(stdout);
  }
}

////////////////////////////////////////////////////////////////////////////////

void c_boxfit :: lightcurve()
{
  int i;
  
  //----------------------------------------------------------------------------
  // setup variables
  
  Obs.nu = nu_0; // light curve entirely observed at single frequency
  prepare();

  //----------------------------------------------------------------------------
  // calculate the actual light curve for serial computation, using openMP if
  // available

  #if OPEN_MPI_ == DISABLED_

    for (i = 0; i < no_points; i++)
    {
      Obs.t = t_0 * pow(t_1 / t_0, (double) i / (double) (no_points - 1));
    
      flux.counter = i; 
    
      if (flux.set_flux() != 0)
        printf("# WARNING, failure to set flux for data point %d\n", i);
      fluxes[i] = flux.F * cgs2mJy;
 
      printf("# %d, %e, %e, %e\n", i, Obs.t, Obs.nu, fluxes[i]); fflush(stdout);
    }

  #endif // OPEN_MPI_ == DISABLED_
  
  //----------------------------------------------------------------------------
  // Second possibility, parallel computation
  
  #if OPEN_MPI_ == ENABLED_

    int idle_id;
    MPI_Status err;
    int stopsignal = -1;
    double Ftemp[no_points]; // temp buffer to collect results
    
    // clean out Ftemp and Fresult
    for (i = 0; i < no_points; i++)
      { Ftemp[i] = 0.0; fluxes[i] = 0.0; }
  
    // myid 0 divides the different datapoints over all the processors
    if (myid == 0) 
    {
      for (i = 0; i < no_points; i++) 
      {
        // wait for a proc signalling that it is idle
        MPI_Recv(&idle_id, 1, MPI_INT, MPI_ANY_SOURCE, 10, MPI_COMM_WORLD,&err);
      
        // send the number of the datapoint to this machine
        MPI_Ssend(&i, 1, MPI_INT, idle_id, 20, MPI_COMM_WORLD);
      }
    
      // tell cores that are idle to stop signaling and expecting signals
      // by sending them a stop signal -1 instead of a datapoint number
      for (i = 1; i < numprocs; i++) 
      {
        // wait for a proc signalling that it is idle
        MPI_Recv(&idle_id, 1, MPI_INT, MPI_ANY_SOURCE, 10, MPI_COMM_WORLD,&err);
      
        // send special signal (-1) telling it to stop 
        MPI_Ssend(&stopsignal, 1, MPI_INT, idle_id, 20, MPI_COMM_WORLD);      
      }
    }
    else // any nonzero proc
    {
      for ( ; ; ) // infinite loop -until broken
      {
        // tell proc 0 that your ID is available for computations
        MPI_Ssend(&myid, 1, MPI_INT, 0, 10, MPI_COMM_WORLD);
      
        // receive the number of the datapoint to start working on
        MPI_Recv(&i, 1, MPI_INT, 0, 20, MPI_COMM_WORLD, &err);
      
        if (i >= 0)
        {
          Obs.t = t_0 * pow(t_1 / t_0, (double) i / (double) (no_points - 1));
          
          flux.counter = i;

          if (flux.set_flux() != 0)
            printf("# WARNING, failure to set flux for data point %d\n", i);

          Ftemp[i] = flux.F * cgs2mJy;
        }
        else
        {
          break;
        }
      }
    }
    
    // We now have an array Ftemp at each processor, filled with zeroes
    // except, for each processor, the ones that particular processor has
    // calculated. We collect the sums of all individual entries at myid = 0, 
    // with for each datapoint only a single processor contributing a nonzero 
    // term.
    MPI_Barrier(MPI_COMM_WORLD); // wait until all cores have caught up
  
    MPI_Reduce(&Ftemp[0], &fluxes[0], no_points, MPI_DOUBLE, MPI_SUM, 
      0, MPI_COMM_WORLD);

  #endif // OPEN_MPI_ == ENABLED_
}

////////////////////////////////////////////////////////////////////////////////

void c_boxfit :: spectrum()
{
  int i;
  
  //----------------------------------------------------------------------------
  // setup variables
  
  Obs.t = t_0; // use same time for entire spectrum
  prepare();

  //----------------------------------------------------------------------------
  // calculate the actual light curve for serial computation

  #if OPEN_MPI_ == DISABLED_

    for (i = 0; i < no_points; i++)
    {
      Obs.nu = nu_0 * pow(nu_1 / nu_0, (double) i / (double) (no_points - 1));
    
      flux.set_flux();
      fluxes[i] = flux.F * cgs2mJy;

      printf("# %d, %e, %e, %e\n", i, Obs.t, Obs.nu, fluxes[i]); fflush(stdout);
    }

  #endif // OPEN_MPI_ == DISABLED_
  
  //----------------------------------------------------------------------------
  // Second possibility, parallel computation
  
  #if OPEN_MPI_ == ENABLED_

    int idle_id;
    MPI_Status err;
    int stopsignal = -1;
    double Ftemp[no_points]; // temp buffer to collect results
    
    // clean out Ftemp and Fresult
    for (i = 0; i < no_points; i++)
      { Ftemp[i] = 0.0; fluxes[i] = 0.0; }
  
    // myid 0 divides the different datapoints over all the processors
    if (myid == 0) 
    {
      for (i = 0; i < no_points; i++) 
      {
        // wait for a proc signalling that it is idle
        MPI_Recv(&idle_id, 1, MPI_INT, MPI_ANY_SOURCE, 10, MPI_COMM_WORLD,&err);
      
        // send the number of the datapoint to this machine
        MPI_Ssend(&i, 1, MPI_INT, idle_id, 20, MPI_COMM_WORLD);
      }
    
      // tell cores that are idle to stop signaling and expecting signals
      // by sending them a stop signal -1 instead of a datapoint number
      for (i = 1; i < numprocs; i++) 
      {
        // wait for a proc signalling that it is idle
        MPI_Recv(&idle_id, 1, MPI_INT, MPI_ANY_SOURCE, 10, MPI_COMM_WORLD,&err);
      
        // send special signal (-1) telling it to stop 
        MPI_Ssend(&stopsignal, 1, MPI_INT, idle_id, 20, MPI_COMM_WORLD);      
      }
    }
    else // any nonzero proc
    {
      for ( ; ; ) // infinite loop -until broken
      {
        // tell proc 0 that your ID is available for computations
        MPI_Ssend(&myid, 1, MPI_INT, 0, 10, MPI_COMM_WORLD);
      
        // receive the number of the datapoint to start working on
        MPI_Recv(&i, 1, MPI_INT, 0, 20, MPI_COMM_WORLD, &err);
      
        if (i >= 0)
        {
          Obs.nu = nu_0 * 
            pow(nu_1 / nu_0, (double) i / (double) (no_points - 1));
    
          flux.set_flux();
          Ftemp[i] = flux.F * cgs2mJy;
        }
        else
        {
          break;
        }
      }
    }
    
    // We now have an array Ftemp at each processor, filled with zeroes
    // except, for each processor, the ones that particular processor has
    // calculated. We collect the sums of all individual entries at myid = 0, 
    // with for each datapoint only a single processor contributing a nonzero 
    // term.
    MPI_Barrier(MPI_COMM_WORLD); // wait until all cores have caught up
  
    MPI_Reduce(&Ftemp[0], &fluxes[0], no_points, MPI_DOUBLE, MPI_SUM, 
      0, MPI_COMM_WORLD);

  #endif // OPEN_MPI_ == ENABLED_
}

////////////////////////////////////////////////////////////////////////////////

void c_boxfit :: calculate_from_derivatives()
{
  int i, j;
  
  for (i = 0; i < dataset.n; i++)
  {
    dataset.Fbox[i] = dataset.logF[i]; // set zeroth order
    for (j = 0; j < thawed; j++)
      dataset.Fbox[i] += dataset.dlogFdx[i][fitsubentry[j]] *
        (fitvar[fitsubentry[j]] - fitvar_initial[fitsubentry[j]]);

    dataset.Fbox[i] = pow(10.0, dataset.Fbox[i]);
  }
}

////////////////////////////////////////////////////////////////////////////////

void c_boxfit :: calculate()
{
  if (fromderivatives) 
  {
    calculate_from_derivatives();
    return;
  }
  
  int i; // dummy index loop over data points
  
  prepare();
  
  //----------------------------------------------------------------------------
  
  #if OPEN_MPI_ == DISABLED_

    for (i = 0; i < dataset.n; i++)
    {
      // set observer frequeny and time
      Obs.t = dataset.t[i];
      Obs.nu = dataset.nu[i];
      
      flux.set_flux();
      dataset.Fbox[i] = flux.F * cgs2mJy;
    }

  #endif // OPEN_MPI_ == DISABLED_
  
  //----------------------------------------------------------------------------
  // Second possibility, parallel computation
  
  #if OPEN_MPI_ == ENABLED_

    int idle_id;
    MPI_Status err;
    int stopsignal = -1;
    double Ftemp[dataset.n]; // temp buffer to collect results
  
    // clear out dataset.Fbox and Ftemp on all processors
    for (i=0; i < dataset.n; i++)
    {
      dataset.Fbox[i] = 0.0;
      Ftemp[i] = 0.0;
    }

    // myid 0 divides the different datapoints over all the processors
    if (myid == 0) 
    {
      for (i = 0; i < dataset.n; i++) 
      {
        // wait for a proc signalling that it is idle
        MPI_Recv(&idle_id, 1, MPI_INT, MPI_ANY_SOURCE, 10, MPI_COMM_WORLD,&err);
      
        // send the number of the datapoint to this machine
        MPI_Ssend(&i, 1, MPI_INT, idle_id, 20, MPI_COMM_WORLD);
      }
    
      // tell cores that are idle to stop signaling and expecting signals
      // by sending them a stop signal -1 instead of a datapoint number
      for (i = 1; i < numprocs; i++) 
      {
        // wait for a proc signalling that it is idle
        MPI_Recv(&idle_id, 1, MPI_INT, MPI_ANY_SOURCE, 10, MPI_COMM_WORLD,&err);
      
        // send special signal (-1) telling it to stop 
        MPI_Ssend(&stopsignal, 1, MPI_INT, idle_id, 20, MPI_COMM_WORLD);      
      }
    }
    else // any nonzero proc
    {
      for ( ; ; ) // infinite loop -until broken
      {
        // tell proc 0 that your ID is available for computations
        MPI_Ssend(&myid, 1, MPI_INT, 0, 10, MPI_COMM_WORLD);
      
        // receive the number of the datapoint to start working on
        MPI_Recv(&i, 1, MPI_INT, 0, 20, MPI_COMM_WORLD, &err);
      
        if (i >= 0)
        {
          // set observer frequeny and time
          Obs.t = dataset.t[i];
          Obs.nu = dataset.nu[i];
          
          flux.set_flux();
          Ftemp[i] = flux.F * cgs2mJy;
        }
        else
        {
          break;
        }
      }
    }
    
    // We now have an array Ftemp at each processor, filled with zeroes
    // except, for each processor, the ones that particular processor has
    // calculated. We collect the sums of all individual entries at myid = 0, 
    // with for each datapoint only a single processor contributing a nonzero 
    // term.
    MPI_Barrier(MPI_COMM_WORLD); // wait until all cores have caught up
  
    MPI_Reduce(&Ftemp[0], &(dataset.Fbox)[0], dataset.n, MPI_DOUBLE, MPI_SUM, 
      0, MPI_COMM_WORLD);

  #endif // OPEN_MPI_ == ENABLED_
  
  // increase number of light curve calculations counter
  lccounter++;
}

////////////////////////////////////////////////////////////////////////////////

double c_boxfit :: set_chisqr()
{
  // calculate chi^2 based on the most recent box light curve. The full results
  // for these light curve are known only to core 0, so chi^2 has to be
  // calculated there and broadcasted to the others.
  double x = 0.0;
  double xi; // for individual measurement
  int i;
  
  if (myid == 0)
  {
    for (i = 0; i < dataset.n; i++)
    {
      if (!usealt) // use real data, not artificially perturbed data
      {
        xi = pow((dataset.F[i] - dataset.Fbox[i]) / dataset.sigma[i], 2);
      }
      else // use perturbed fluxes based on real data & error, for MonteCarlo
        xi = pow((dataset.Falt[i] - dataset.Fbox[i]) / dataset.sigma[i], 2);
      x += xi;
    }
  }
  
  // deal with the other cores if necessary
  #if OPEN_MPI_ == ENABLED_
    MPI_Bcast(&x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
  #endif // OPEN_MPI_ == ENABLED_
  
  chi_sqr = x;
  
  if (myid == 0 and p_chi_sqr != NULL)
  {
    print_fit_result(p_chi_sqr);
  }
  
  return x;
}

////////////////////////////////////////////////////////////////////////////////

void c_boxfit :: initialize_simplex()
{
  int i_sim, i_var;

  for (i_sim = 0; i_sim <= thawed; i_sim++)
  {
    for (i_var = 0; i_var < thawed; i_var++)
    {
      initial_simplex[i_sim][i_var] = fitvar_initial[fitsubentry[i_var]];
      if (i_var == i_sim - 1)
        initial_simplex[i_sim][i_var] = simplex_max[fitsubentry[i_var]];
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

int c_boxfit :: fit_dataset()
{
  FILE* p_fit = NULL; // pointer to fit variable evolution file
  int imcount = 0; // number of intermediate files stored
  int i, j; // dummy loop indices
  
  char savefilename[200]; // name of intermediate and final best fit files
  FILE* p_file; // pointer to intermediate and final best fit files
  int success = 0; // track error flags
  
  //----------------------------------------------------------------------------
  
  // allocate memory, open file
  if (myid == 0 and !quietfit) 
  { 
    p_fit = fopen("fit.txt", "w");
    print_opening_message(p_fit);
    print_settings(p_fit);
    fprintf(p_fit, "i, temp, theta_0, E, n, theta_obs, p, epsilon_B, "
      "epsilon_E, ksi_N, chi^2, red. chi^2\n");
    fflush(stdout);
    fflush(p_fit);
    p_chi_sqr = p_fit;
  }
  
  //----------------------------------------------------------------------------

  if (start_from_simplex)
  {
    simplex.load(simplex_filename);
    for(i = 0; i <= thawed; i++)
    {
      for(j = 0; j < thawed; j++)
      {
        initial_simplex[i][j] = simplex.x[i][j];
      }
      simplex.evaluate(i);
    }
  }
  else
  {
    for(i = 0; i <= thawed; i++)
    {
      for(j = 0; j < thawed; j++)
      {
        simplex.x[i][j] = initial_simplex[i][j];
      }
      simplex.evaluate(i);
    }
  }

  if (myid == 0 and !quietfit) 
  {
    fprintf(p_fit, "# initial simplex done\n");
    fflush(p_fit);
  }

  //----------------------------------------------------------------------------

  // Downhill simplex plus annealing at temperature > 0
  if (temp > 0.0)
  {
    simplex.max_iter = iter;
    simplex.temp = temp;
    
    do
    {
      simplex.minimize();
      if (myid == 0) simplex.save("currentsimplex.txt");

      // store best fit result in fit_var;
      for (i = 0; i < thawed; i++)
        fitvar[fitsubentry[i]] = simplex.x[simplex.i_l][i];
      chi_sqr = simplex.y[simplex.i_l];
  
      if (myid == 0 and !quietfit) 
      { 
        printf("# run at temp %e completed. Best fit: \n", simplex.temp);
        print_fit_result(stdout);
        print_fit_result(p_fit);
      }
      
      if (myid == 0 and !quietfit)
      {
        if (save_intermediate)
        {
          sprintf(savefilename, "intermediate%04d.txt", imcount);
          p_file = fopen(savefilename, "w");
          print_opening_message(p_file);
          print_settings(p_file);
          fprintf(p_file, "# Current temperature: %e\n", simplex.temp);
          fprintf(p_file, "# number of light curves calculated: %d\n", 
            lccounter);
          success = dataset.print(p_file);
          fclose(p_file);
          imcount++;
        }
      }
      
      #if OPEN_MPI_ == ENABLED_
        MPI_Barrier(MPI_COMM_WORLD); // wait until all cores have caught up
        MPI_Bcast(&success, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD); // wait until all cores have caught up
      #endif
      if (success != 0) 
      {
        if (myid == 0) fclose(p_fit);
        return success;
      }

      simplex.temp = simplex.temp * temp_factor;

    } while(simplex.temp >= temp_lowest);
  }
  
  //----------------------------------------------------------------------------
  
  // finalize, straightforward downhill simplex at t = 0, iter iterations
  simplex.max_iter = iter;
  simplex.temp = 0.0;
  simplex.minimize();
  if (myid == 0) simplex.save("currentsimplex.txt");
  
  // store best fit result in fit_var;
  for (i = 0; i < thawed; i++)
    fitvar[fitsubentry[i]] = simplex.x[simplex.i_l][i];
  chi_sqr = simplex.y[simplex.i_l];
  
  if (myid == 0) 
  { 
    if (!quietfit) 
    {
      printf("# Full annealing completed. Best fit: \n");
      print_fit_result(p_fit);
    }
    print_fit_result(stdout);
  
    sprintf(savefilename, "bestfit.txt");
    p_file = fopen(savefilename, "w");
    print_opening_message(p_file);
    print_settings(p_file);
    fprintf(p_file, "# Current temperature: %e\n", simplex.temp);
    fprintf(p_file, "# number of light curves calculated: %d\n", lccounter);
    dataset.print(p_file);
    fclose(p_file);
  }
  
  //----------------------------------------------------------------------------
  
  // release memory, close file
  if (myid == 0 and !quietfit)
  {
    fclose(p_fit);
  }
  
  return success;
}

////////////////////////////////////////////////////////////////////////////////

double c_boxfit :: get_average_error()
{
  int i;
  double err = 0.0;
  
  for (i = 0; i < dataset.n; i++)
    err += dataset.sigma[i] / dataset.F[i];
    
  return err / (double) dataset.n;
}

////////////////////////////////////////////////////////////////////////////////

int c_boxfit :: setup_derivatives()
{
  int i_dat, i_var, i_var2; // dummy loop indices
  double dx; // auxiliary variables to calculate derivatives
  FILE *p_file; // pointer to output file
  double Wup = 1.0, Wdown = 1.0; 
    // weight of derivatives on each side of datapoint. If one of two drops out,
    // Wup or Wdown will be set to 0.5 instead, resulting in a final Delta X
    // (with X the fit variable) that is only half times as big. The situation
    // where both Wup and Wdown are 0.5 should not occur: the allowed range for
    // the fit variable is then too small
  
  //----------------------------------------------------------------------------
  
  // set the average fractional error on the dataset
  dataset.ferr = get_average_error();
  if (myid == 0) 
  {
    printf("# average error on dataset = %e\n", dataset.ferr);
    fflush(stdout);
  }
  
  // calculate the light curves at the dataset values.
  // it is assumed that the fit values provided by boxfitsettings.txt are the
  // best fit values.
  for (i_var = 0; i_var < no_fit_vars_; i_var++)
    fitvar[i_var] = fitvar_initial[i_var];
  calculate();
  
  // store result
  if (myid == 0)
  {
    p_file = fopen("bestfit.txt", "w");
    print_opening_message(p_file);
    print_settings(p_file);
    dataset.print(p_file);
    fclose(p_file);
  }
  
  // store their logs in lnF
  for (i_dat = 0; i_dat < dataset.n; i_dat++)
    dataset.logF[i_dat] = log10(dataset.Fbox[i_dat]);
  
  // clean out the partial derivatives for all variables, including frozen
  for (i_var = 0; i_var < no_fit_vars_; i_var++)
    for (i_dat = 0; i_dat < dataset.n; i_dat++)
      dataset.dlogFdx[i_dat][i_var] = 0.0;
  
  for (i_var = 0; i_var < thawed; i_var++)
  {
    // reset all fit variables, before altering one
    for (i_var2 = 0; i_var2 < no_fit_vars_; i_var2++)
      fitvar[i_var2] = fitvar_initial[i_var2];
    
    // set stepsize from parameter file
    dx = 0.5 * 
      (simplex_max[fitsubentry[i_var]] - fitvar_initial[fitsubentry[i_var]]);
      
    // move down if possible within allowed parameter range
    if (fitvar_initial[fitsubentry[i_var]] - 2.0 * dx > 
      varmin[fitsubentry[i_var]])
    {
      // two steps down
      fitvar[fitsubentry[i_var]] = fitvar_initial[fitsubentry[i_var]] - 
        2.0 * dx;
      calculate();
    
      for (i_dat = 0; i_dat < dataset.n; i_dat++)
        dataset.dlogFdx[i_dat][fitsubentry[i_var]] += 
          -2.0 / 10.0 * log10(dataset.Fbox[i_dat]);
    
      // one step down
      fitvar[fitsubentry[i_var]] = fitvar_initial[fitsubentry[i_var]] 
        - 1.0 * dx;
      calculate();

      for (i_dat = 0; i_dat < dataset.n; i_dat++)
        dataset.dlogFdx[i_dat][fitsubentry[i_var]] += 
          -1.0 / 10.0 * log10(dataset.Fbox[i_dat]);
    }
    else
    {
      Wdown = 0.5;
    }

    // move up if possible within allowed parameter range
    if (fitvar_initial[fitsubentry[i_var]] + 2.0 * dx < 
      varmax[fitsubentry[i_var]])
    {
      // one step up
      fitvar[fitsubentry[i_var]] = fitvar_initial[fitsubentry[i_var]] 
        + 1.0 * dx;
      calculate();

      for (i_dat = 0; i_dat < dataset.n; i_dat++)
        dataset.dlogFdx[i_dat][fitsubentry[i_var]] += 
          1.0 / 10.0 * log10(dataset.Fbox[i_dat]);
    
      // two steps up
      fitvar[fitsubentry[i_var]] = fitvar_initial[fitsubentry[i_var]] 
        + 2.0 * dx;
      calculate();

      for (i_dat = 0; i_dat < dataset.n; i_dat++)
        dataset.dlogFdx[i_dat][fitsubentry[i_var]] += 
          2.0 / 10.0 * log10(dataset.Fbox[i_dat]);
    }
    else
    {
      Wup = 0.5;
    }
    
    // error check: if neither room to move up or down, bail out
    if (Wup < 0.9 and Wdown < 0.9)
    {
      if (myid == 0)
      {
        printf("ERROR: no room to determine partial derivative for fit "
          "parameter %d\n", fitsubentry[i_var]);
        printf("Please check varmin / varmax range to simplex_max / initial\n");
        fflush(stdout);
        return 1; // return failure
      }
    }
    
    // take stepsize into account
    for (i_dat = 0; i_dat < dataset.n; i_dat++)
      dataset.dlogFdx[i_dat][fitsubentry[i_var]] /= (Wdown * Wup * dx);
  }
  
  // store the resulting partial derivatives
  if (myid == 0) save_derivatives();
  
  return 0; // return success
}

////////////////////////////////////////////////////////////////////////////////

void c_boxfit :: load_derivatives()
{
  FILE *p_file;
  int i_dat;
  int success; // checks number of succesfully scanned variables from file
  
  p_file = fopen("partialdevs.txt", "r");
  
  for (i_dat = 0; i_dat < dataset.n; i_dat++)
  {
    success = fscanf(p_file, "%d, %le, %le, %le, %le, %le, %le, %le, %le, "
      "%le, %le, %le\n", &i_dat, &dataset.t[i_dat], &dataset.nu[i_dat], 
      &dataset.logF[i_dat], &(dataset.dlogFdx[i_dat][fit_theta_0_]),
      &(dataset.dlogFdx[i_dat][fit_E_]),
      &(dataset.dlogFdx[i_dat][fit_n_]),
      &(dataset.dlogFdx[i_dat][fit_theta_obs_]),
      &(dataset.dlogFdx[i_dat][fit_p_]),
      &(dataset.dlogFdx[i_dat][fit_epsilon_B_]),
      &(dataset.dlogFdx[i_dat][fit_epsilon_E_]),
      &(dataset.dlogFdx[i_dat][fit_ksi_N_]));
    if (success < 12)
    {
      printf("ERROR reading partialdevs.txt, could only read %d variables out"
        " of 12\n", success);
      fflush(stdout);
    }
  }
  fclose(p_file);
}

////////////////////////////////////////////////////////////////////////////////

void c_boxfit :: save_derivatives()
{
  FILE *p_file;
  int i_dat;
  
  p_file = fopen("partialdevs.txt", "w");
  for (i_dat = 0; i_dat < dataset.n; i_dat++)
  {
    fprintf(p_file, "%d, %e, %e, %e, %1.3e, %1.3e, %1.3e, %1.3e, %1.3e, "
      "%1.3e, %1.3e, %1.3e\n", i_dat, dataset.t[i_dat], dataset.nu[i_dat], 
      dataset.logF[i_dat], dataset.dlogFdx[i_dat][fit_theta_0_],
      dataset.dlogFdx[i_dat][fit_E_],
      dataset.dlogFdx[i_dat][fit_n_],
      dataset.dlogFdx[i_dat][fit_theta_obs_],
      dataset.dlogFdx[i_dat][fit_p_],
      dataset.dlogFdx[i_dat][fit_epsilon_B_],
      dataset.dlogFdx[i_dat][fit_epsilon_E_],
      dataset.dlogFdx[i_dat][fit_ksi_N_]);
  }
  fclose(p_file);
}

////////////////////////////////////////////////////////////////////////////////

void c_boxfit :: store_new_MC()
{
  int slot = 0; // new place in the list
  int i, j;
  
  // find place of new result in the ordered list
  while(slot < MCi and MCchisqr[slot] < simplex.y[simplex.i_l])
    slot++;
  
  // shift later entries downward
  for (i = MCi; i > slot; i--)
  {
    MCchisqr[i] = MCchisqr[i-1];
    for (j = 0; j < no_fit_vars_; j++) 
      MCfitvar[i][j] = MCfitvar[i-1][j];
  }
  
  // add the new entry
  MCchisqr[slot] = simplex.y[simplex.i_l];
  
  // add all the model variables, including those not included in the fit
  for (j = 0; j < no_fit_vars_; j++)
  {
    if (j == fit_p_ or j == fit_theta_0_ or j == fit_theta_obs_)
      MCfitvar[slot][j] = fitvar[j];
    else
      MCfitvar[slot][j] = pow(10.0, fitvar[j]);
  }
  // make sure that the entries for the variables fitted for are the best fit
  for (j = 0; j < thawed; j++)
  {
    if (fitsubentry[j] == fit_p_ or fitsubentry[j] == fit_theta_0_ or 
      fitsubentry[j] == fit_theta_obs_)
      MCfitvar[slot][fitsubentry[j]] = simplex.x[simplex.i_l][j];
    else
      MCfitvar[slot][fitsubentry[j]] = pow(10.0, simplex.x[simplex.i_l][j]);
  }
  
  // every once in a while store results
  if (MCi % 100 == 0 and myid == 0)
    save_MC_result();
}

////////////////////////////////////////////////////////////////////////////////

void c_boxfit :: save_MC_result()
{
  int i, j;
  FILE *p_MCfile;
  
  p_MCfile = fopen("MCresult.txt", "w");
  
  print_opening_message(p_MCfile);
  print_settings(p_MCfile);
  
  for (i = 0; i < MCi; i++)
  {
    fprintf(p_MCfile, "%d, ", MCi);
    for (j = 0; j < no_fit_vars_; j++)
      fprintf(p_MCfile, "%e, ", MCfitvar[i][j]);
    fprintf(p_MCfile, "%e, %e\n", MCchisqr[i], 
      MCchisqr[i] / (dataset.n - thawed));
  }
  fclose(p_MCfile);
}

////////////////////////////////////////////////////////////////////////////////

int c_boxfit :: montecarlo()
{
  int i, j;
  int MCfraction; // fraction of total MC runs that counts towards errors on
    // the fit parameters.
  
  double ferr; // average fractional error in dataset
  FILE *p_MCout = NULL; // pointer Monte Carlo results output file
  FILE *p_synthetic; // pointer to file containing results for first run
  int success = 0; // track error flags
  
  //----------------------------------------------------------------------------
  // set up variables
  
  MCfitvar = array_2d<double>(MC_runs, no_fit_vars_);
  MCchisqr = array_1d<double>(MC_runs);
  MCvarmax = array_1d<double>(no_fit_vars_);
  MCvarmin = array_1d<double>(no_fit_vars_);
 
  usealt = true; // set to use alternative, perturbed fluxes
  quietfit = true; // disable output during fit procedure
  
  ferr = get_average_error();
  if (myid == 0) printf("# average fractional error in dataset: %e\n", ferr);
  
  // setup output files
  if (myid == 0) 
  {
    p_MCout = fopen("MCout.txt", "w");
    print_opening_message(p_MCout);
    print_settings(p_MCout);
    fprintf(p_MCout, "# average fractional error in dataset: %e\n", ferr);
    
    p_synthetic = fopen("synthetic.txt", "w");
    print_opening_message(p_synthetic);
    print_settings(p_synthetic);
  }

  //----------------------------------------------------------------------------
  
  for (MCi = 0; MCi < MC_runs; MCi++)
  {
    // set up new values for Falt
    for (i = 0; i < dataset.n; i++)
      do
      { 
        dataset.Falt[i] = dataset.F[i] + random.globalGaussian() * 
          dataset.sigma[i];
      } while (dataset.Falt[i] < 0.0); // repeat until positive flux value

    success = fit_dataset();
    if (success != 0)
    {
      if (myid == 0)
      {
        fclose(p_MCout);
      }
  
      delete_array_2d(MCfitvar);
      delete_array_1d(MCchisqr);
      delete_array_1d(MCvarmax);
      delete_array_1d(MCvarmin);

      return success;
    }
    
    // store results
    store_new_MC();
  }

  //----------------------------------------------------------------------------
  // calculate the results
  
  save_MC_result();
  
  for (i = 0; i < no_fit_vars_; i++)
  {
    MCvarmax[i] = 0.0;
    MCvarmin[i] = 1e300;
  }
  
  MCfraction = (int) (0.682689492 * MC_runs);
  for (i = 0; i < MCfraction; i++)
  {
    for (j = 0; j < thawed; j++)
    {
      if (MCvarmax[fitsubentry[j]] < MCfitvar[i][fitsubentry[j]]) 
        MCvarmax[fitsubentry[j]] = MCfitvar[i][fitsubentry[j]];
      if (MCvarmin[fitsubentry[j]] > MCfitvar[i][fitsubentry[j]]) 
        MCvarmin[fitsubentry[j]] = MCfitvar[i][fitsubentry[j]];
    }
  }

  // print the best results:
  if (myid == 0)
  {
    printf("MC results, %d drawn from total %d runs:\n", MCfraction, MC_runs);
    printf("MC lower, input best fit, MC mean, MC upper\n");
    
    for (i = 0; i < thawed; i++)
    {
      if (fitsubentry[i] == fit_p_ or fitsubentry[i] == fit_theta_0_ or 
        fitsubentry[i] == fit_theta_obs_)
        printf("%d, %e, %e, %e\n", fitsubentry[i], MCvarmin[fitsubentry[i]], 
          fitvar_initial[fitsubentry[i]], MCvarmax[fitsubentry[i]]);
      else
        printf("%d, %e, %e, %e\n", fitsubentry[i], MCvarmin[fitsubentry[i]], 
          pow(10.0, fitvar_initial[fitsubentry[i]]), MCvarmax[fitsubentry[i]]);
    }
  }
  
  //----------------------------------------------------------------------------
  // free memory, close files
  
  if (myid == 0)
  {
    fclose(p_MCout);
  }
  
  delete_array_2d(MCfitvar);
  delete_array_1d(MCchisqr);
  delete_array_1d(MCvarmax);
  delete_array_1d(MCvarmin);
  
  return success;
}

////////////////////////////////////////////////////////////////////////////////

int c_boxfit :: initialize(int argc, char* argv[])
{
  int i, j; // dummy loop indices
  int box0 = -1; // first box entry number on disc
  int box1 = -1; // last box entry number on disc
  int success = 0; // nonzero if failure in loading or internal checks
  double tempdouble; // temporary storage when reading double from par file
  int tempint, tempint2; // temporary storage when reading int from par file
  bool read_status, read_status2; // true if parameter succesfully read from
    // parameter file or from command line
  
  //----------------------------------------------------------------------------
  // bail out if only one processor used in parallel mode

  #if OPEN_MPI_ == ENABLED_
    if (numprocs == 1)
    {
      printf("ERROR: only 1 core on parallel run. "
        "No cores left for data analysis.\n");
      fflush(stdout);
      return(1);
    }
  #endif // OPEN_MPI_

  //----------------------------------------------------------------------------
  // allocate memory

  fitvar = array_1d<double>(no_fit_vars_);
  fitvar_allocated = true;
  
  fitvar_initial = array_1d<double>(no_fit_vars_);
  fitvar_initial_allocated = true;
  
  varmax = array_1d<double>(no_fit_vars_);
  varmax_allocated = true;
  
  varmin = array_1d<double>(no_fit_vars_);
  varmin_allocated = true;
  
  varfrozen = array_1d<bool>(no_fit_vars_);
  varfrozen_allocated = true;
  
  fitsubentry = array_1d<int>(no_fit_vars_);
  fitsubentry_allocated = true;
  
  initial_simplex = array_2d<double>(no_fit_vars_ + 1, no_fit_vars_ + 1);
  initial_simplex_allocated = true;

  simplex_max = array_1d<double>(no_fit_vars_);
  simplex_max_allocated = true;
  
  // dummy string allocations, required by parameter file I/O functions 
  data_filename = new char[1];
  data_filename_allocated = true;

  simplex_filename = new char[1];
  simplex_filename_allocated = true;

  box_filename_base = new char[1];
  box_filename_base_allocated = true;

  //----------------------------------------------------------------------------
  // set pointers to class containing observer settings
  
  flux.set_observer_pointer(&Obs);

  //----------------------------------------------------------------------------
  // set variables to standard values
  
  // initialize the random number generator, different for each core
  random.initialize(34763 * (myid + 1));
  
  lccounter = 0; // set number of light curves calculated to zero
  usealt = false; // initialize to using standard fluxes F, not Falt
  fromderivatives = false; // don't use derivates approach for light curve
  p_chi_sqr = NULL; // disable printing of each chi^2 result
  
  simplex.p_f = funk_wrapper; // set simplex pointer to wrapper function
  simplex.p_random = &random; // set simplex pointer to random number generator
  
  quietfit = false; // assume were not running MC procedure
  
  //----------------------------------------------------------------------------
  // read the parameter file

  what_to_do = int_from_parfile(parfilename, "what_to_do");
  read_status = parse_int("-what_to_do=", tempint, argc, argv);
  if (read_status)
    what_to_do = tempint;
  
  // read IO related parameters
  if (what_to_do != -1 and what_to_do != 1 and what_to_do != 2)
  {
    string_from_parfile(parfilename, "data_filename", data_filename);
    if (parse("-data_filename=", argc, argv))
      parse_string("-data_filename=", data_filename, argc, argv);
  }

  read_status = read_parfile_string(parfilename, "box_filename_base", 
    box_filename_base);
  if (read_status == false and !parse("-box_filename_base=", argc, argv))
  {
    if (myid == 0)
    {
      printf("ERROR: box_filename_base not given in parameter file or command "
        "line.\n");
      fflush(stdout);
    }
    return 1;
  }
  
  if (parse("-box_filename_base=", argc, argv)) // command line overrules
    parse_string("-box_filename_base=", box_filename_base, argc, argv);

  if (int_from_parfile(parfilename, "save_intermediate") == 1)
    save_intermediate = true;
  else
    save_intermediate = false;
  read_status = parse_int("-save_intermediate=", tempint, argc, argv);
  if (read_status)
  {
    if (tempint == 0)
      save_intermediate = false;
    else
      save_intermediate = true;
  }

  // read intermediate log file settings
  read_status = read_parfile_int(parfilename, "save_emission_profile", tempint);
  read_status2 = parse_int("-save_emission_profile=", tempint2, argc, argv);
  if (what_to_do == 1 or what_to_do == 2) 
  {
    if (read_status == false and read_status2 == false and myid == 0)
    {
      printf("# WARNING, save_emission_profile not set, assuming 0.\n");
      fflush(stdout);
    }
    if (read_status == true)
    {
      if (tempint == 1) flux.save_emission_profile = true;
    }
    if (read_status2 == true) // command line overrules
    {
      if (tempint2 == 1)
        flux.save_emission_profile = true;
      else
        flux.save_emission_profile = false;
    }
  }

  // read intermediate image-on-sky file settings
  read_status = read_parfile_int(parfilename, "save_image", tempint);
  read_status2 = parse_int("-save_image=", tempint2, argc, argv);
  if (what_to_do == 1 or what_to_do == 2) 
  {
    if (read_status == false and read_status2 == false and myid == 0)
    {
      // feature is 'hidden' for now, while under development
      //printf("# WARNING, save_emission_profile not set, assuming 0.\n");
      //fflush(stdout);
    }
    if (read_status == true)
    {
      if (tempint == 1) flux.save_image = true;
    }
    if (read_status2 == true) // command line overrules
    {
      if (tempint2 == 1)
        flux.save_image = true;
      else
        flux.save_image = false;
    }
  }

  // read frequencies and times
  t_0 = double_from_parfile(parfilename, "t_0");
  if (parse("-t_0=", argc, argv)) // command line overrules
    parse_double("-t_0=", t_0, argc, argv);
  
  t_1 = double_from_parfile(parfilename, "t_1");
  if (parse("-t_1=", argc, argv)) // command line overrules
    parse_double("-t_1=", t_1, argc, argv);

  nu_0 = double_from_parfile(parfilename, "nu_0");
  if (parse("-nu_0=", argc, argv)) // command line overrules
    parse_double("-nu_0=", nu_0, argc, argv);

  nu_1 = double_from_parfile(parfilename, "nu_1");
  if (parse("-nu_1=", argc, argv)) // command line overrules
    parse_double("-nu_1=", nu_1, argc, argv);

  // read physics related parameters
  Obs.dL = double_from_parfile(parfilename, "d_L");
  if (parse("-d_L=", argc, argv)) // command line overrules
    parse_double("-d_L=", Obs.dL, argc, argv);

  Obs.z = double_from_parfile(parfilename, "z");
  if (parse("-z=", argc, argv)) // command line overrules
    parse_double("-z=", Obs.z, argc, argv);

  // read initial fit parameters. Note that all except p and angles are 
  // translated into log space
  fitvar[fit_theta_0_] = double_from_parfile(parfilename, "theta_0");
  if (parse("-theta_0=", argc, argv)) // command line overrules
    parse_double("-theta_0=", fitvar[fit_theta_0_], argc, argv);

  fitvar[fit_E_] = log10(double_from_parfile(parfilename, "E"));
  if (parse("-E=", argc, argv)) // command line overrules
  {
    parse_double("-E=", fitvar[fit_E_], argc, argv);
    fitvar[fit_E_] = log10(fitvar[fit_E_]);
  }

  fitvar[fit_n_] = log10(double_from_parfile(parfilename, "n"));
  if (parse("-n=", argc, argv)) // command line overrules
  {
    parse_double("-n=", fitvar[fit_n_], argc, argv);
    fitvar[fit_n_] = log10(fitvar[fit_n_]);
  }

  fitvar[fit_theta_obs_] = double_from_parfile(parfilename, "theta_obs");
  if (parse("-theta_obs=", argc, argv)) // command line overrules
    parse_double("-theta_obs=", fitvar[fit_theta_obs_], argc, argv);
  
  fitvar[fit_p_] = double_from_parfile(parfilename, "p");
  if (parse("-p=", argc, argv)) // command line overrules
    parse_double("-p=", fitvar[fit_p_], argc, argv);

  fitvar[fit_epsilon_B_] = 
    log10(double_from_parfile(parfilename, "epsilon_B"));
  if (parse("-epsilon_B=", argc, argv)) // command line overrules
  {
    parse_double("-epsilon_B=", fitvar[fit_epsilon_B_], argc, argv);
    fitvar[fit_epsilon_B_] = log10(fitvar[fit_epsilon_B_]);
  }

  fitvar[fit_epsilon_E_] = 
    log10(double_from_parfile(parfilename, "epsilon_E"));
  if (parse("-epsilon_E=", argc, argv)) // command line overrules
  {
    parse_double("-epsilon_E=", fitvar[fit_epsilon_E_], argc, argv);
    fitvar[fit_epsilon_E_] = log10(fitvar[fit_epsilon_E_]);
  }

  fitvar[fit_ksi_N_] = 
    log10(double_from_parfile(parfilename, "ksi_N"));
  if (parse("-ksi_N=", argc, argv)) // command line overrules
  {
    parse_double("-ksi_N=", fitvar[fit_ksi_N_], argc, argv);
    fitvar[fit_ksi_N_] = log10(fitvar[fit_ksi_N_]);
  }
    
  for (i = 0; i < no_fit_vars_; i++) fitvar_initial[i] = fitvar[i];
  
  //----------------------------------------------------------------------------
  
  varmax[fit_theta_0_] = double_from_parfile(parfilename, "theta_0_max");
  if (parse("-theta_0_max=", argc, argv)) // command line overrules
    parse_double("-theta_0_max=", varmax[fit_theta_0_], argc, argv);

  varmax[fit_E_] = log10(double_from_parfile(parfilename, "E_max"));
  if (parse("-E_max=", argc, argv)) // command line overrules
  {
    parse_double("-E_max=", varmax[fit_E_], argc, argv);
    varmax[fit_E_] = log10(varmax[fit_E_]);
  }

  varmax[fit_n_] = log10(double_from_parfile(parfilename, "n_max"));
  if (parse("-n_max=", argc, argv)) // command line overrules
  {
    parse_double("-n_max=", varmax[fit_n_], argc, argv);
    varmax[fit_n_] = log10(varmax[fit_n_]);
  }

  varmax[fit_theta_obs_] = double_from_parfile(parfilename, "theta_obs_max");
  if (parse("-theta_obs_max=", argc, argv)) // command line overrules
    parse_double("-theta_obs_max=", varmax[fit_theta_obs_], argc, argv);

  varmax[fit_p_] = double_from_parfile(parfilename, "p_max");
  if (parse("-p_max=", argc, argv)) // command line overrules
    parse_double("-p_max=", varmax[fit_p_], argc, argv);

  varmax[fit_epsilon_B_] = 
    log10(double_from_parfile(parfilename, "epsilon_B_max"));
  if (parse("-epsilon_B_max=", argc, argv)) // command line overrules
  {
    parse_double("-epsilon_B_max=", varmax[fit_epsilon_B_], argc, argv);
    varmax[fit_epsilon_B_] = log10(varmax[fit_epsilon_B_]);
  }

  varmax[fit_epsilon_E_] = 
    log10(double_from_parfile(parfilename, "epsilon_E_max"));
  if (parse("-epsilon_E_max=", argc, argv)) // command line overrules
  {
    parse_double("-epsilon_E_max=", varmax[fit_epsilon_E_], argc, argv);
    varmax[fit_epsilon_E_] = log10(varmax[fit_epsilon_E_]);
  }

  varmax[fit_ksi_N_] = 
    log10(double_from_parfile(parfilename, "ksi_N_max"));
  if (parse("-ksi_N_max=", argc, argv)) // command line overrules
  {
    parse_double("-ksi_N_max=", varmax[fit_ksi_N_], argc, argv);
    varmax[fit_ksi_N_] = log10(varmax[fit_ksi_N_]);
  }
  
  //----------------------------------------------------------------------------

  varmin[fit_theta_0_] = double_from_parfile(parfilename, "theta_0_min");
  if (parse("-theta_0_min=", argc, argv)) // command line overrules
    parse_double("-theta_0_min=", varmin[fit_theta_0_], argc, argv);

  varmin[fit_E_] = log10(double_from_parfile(parfilename, "E_min"));
  if (parse("-E_min=", argc, argv)) // command line overrules
  {
    parse_double("-E_min=", varmin[fit_E_], argc, argv);
    varmin[fit_E_] = log10(varmin[fit_E_]);
  }

  varmin[fit_n_] = log10(double_from_parfile(parfilename, "n_min"));
  if (parse("-n_min=", argc, argv)) // command line overrules
  {
    parse_double("-n_min=", varmin[fit_n_], argc, argv);
    varmin[fit_n_] = log10(varmin[fit_n_]);
  }

  varmin[fit_theta_obs_] = double_from_parfile(parfilename, "theta_obs_min");
  if (parse("-theta_obs_min=", argc, argv)) // command line overrules
    parse_double("-theta_obs_min=", varmin[fit_theta_obs_], argc, argv);

  varmin[fit_p_] = double_from_parfile(parfilename, "p_min");
  if (parse("-p_min=", argc, argv)) // command line overrules
    parse_double("-p_min=", varmin[fit_p_], argc, argv);

  varmin[fit_epsilon_B_] = 
    log10(double_from_parfile(parfilename, "epsilon_B_min"));
  if (parse("-epsilon_B_min=", argc, argv)) // command line overrules
  {
    parse_double("-epsilon_B_min=", varmin[fit_epsilon_B_], argc, argv);
    varmin[fit_epsilon_B_] = log10(varmin[fit_epsilon_B_]);
  }

  varmin[fit_epsilon_E_] = 
    log10(double_from_parfile(parfilename, "epsilon_E_min"));
  if (parse("-epsilon_E_min=", argc, argv)) // command line overrules
  {
    parse_double("-epsilon_E_min=", varmin[fit_epsilon_E_], argc, argv);
    varmin[fit_epsilon_E_] = log10(varmin[fit_epsilon_E_]);
  }

  varmin[fit_ksi_N_] = 
    log10(double_from_parfile(parfilename, "ksi_N_min"));
  if (parse("-ksi_N_min=", argc, argv)) // command line overrules
  {
    parse_double("-ksi_N_min=", varmin[fit_ksi_N_], argc, argv);
    varmin[fit_ksi_N_] = log10(varmin[fit_ksi_N_]);
  }

  //----------------------------------------------------------------------------
  // read fit settings

  if (int_from_parfile(parfilename, "theta_0_frozen") == 1)
    varfrozen[fit_theta_0_] = true;
  else
    varfrozen[fit_theta_0_] = false;
  read_status = parse_int("-theta_0_frozen=", tempint, argc, argv);
  if (read_status)
  {
    if (tempint == 0)
      varfrozen[fit_theta_0_] = false;
    else
      varfrozen[fit_theta_0_] = true;
  }
  
  if (int_from_parfile(parfilename, "E_frozen") == 1)
    varfrozen[fit_E_] = true;
  else
    varfrozen[fit_E_] = false;
  read_status = parse_int("-E_frozen=", tempint, argc, argv);
  if (read_status)
  {
    if (tempint == 0)
      varfrozen[fit_E_] = false;
    else
      varfrozen[fit_E_] = true;
  }

  if (int_from_parfile(parfilename, "n_frozen") == 1)
    varfrozen[fit_n_] = true;
  else
    varfrozen[fit_n_] = false;
  read_status = parse_int("-n_frozen=", tempint, argc, argv);
  if (read_status)
  {
    if (tempint == 0)
      varfrozen[fit_n_] = false;
    else
      varfrozen[fit_n_] = true;
  }

  if (int_from_parfile(parfilename, "theta_obs_frozen") == 1)
    varfrozen[fit_theta_obs_] = true;
  else
    varfrozen[fit_theta_obs_] = false;
  read_status = parse_int("-theta_obs_frozen=", tempint, argc, argv);
  if (read_status)
  {
    if (tempint == 0)
      varfrozen[fit_theta_obs_] = false;
    else
      varfrozen[fit_theta_obs_] = true;
  }

  if (int_from_parfile(parfilename, "p_frozen") == 1)
    varfrozen[fit_p_] = true;
  else
    varfrozen[fit_p_] = false;
  read_status = parse_int("-p_frozen=", tempint, argc, argv);
  if (read_status)
  {
    if (tempint == 0)
      varfrozen[fit_p_] = false;
    else
      varfrozen[fit_p_] = true;
  }

  if (int_from_parfile(parfilename, "epsilon_B_frozen") == 1)
    varfrozen[fit_epsilon_B_] = true;
  else
    varfrozen[fit_epsilon_B_] = false;
  read_status = parse_int("-epsilon_B_frozen=", tempint, argc, argv);
  if (read_status)
  {
    if (tempint == 0)
      varfrozen[fit_epsilon_B_] = false;
    else
      varfrozen[fit_epsilon_B_] = true;
  }

  if (int_from_parfile(parfilename, "epsilon_E_frozen") == 1)
    varfrozen[fit_epsilon_E_] = true;
  else
    varfrozen[fit_epsilon_E_] = false;
  read_status = parse_int("-epsilon_E_frozen=", tempint, argc, argv);
  if (read_status)
  {
    if (tempint == 0)
      varfrozen[fit_epsilon_E_] = false;
    else
      varfrozen[fit_epsilon_E_] = true;
  }

  if (int_from_parfile(parfilename, "ksi_N_frozen") == 1)
    varfrozen[fit_ksi_N_] = true;
  else
    varfrozen[fit_ksi_N_] = false;
  read_status = parse_int("-ksi_N_frozen=", tempint, argc, argv);
  if (read_status)
  {
    if (tempint == 0)
      varfrozen[fit_ksi_N_] = false;
    else
      varfrozen[fit_ksi_N_] = true;
  }

  //----------------------------------------------------------------------------
  // read initial simplex width from file
  
  simplex_max[fit_theta_0_] = 
    double_from_parfile(parfilename, "theta_0_simplex_max");
  if (parse("-theta_0_simplex_max=", argc, argv)) // command line overrules
    parse_double("-theta_0_simplex_max=", simplex_max[fit_theta_0_],
      argc, argv);

  simplex_max[fit_E_] = 
    log10(double_from_parfile(parfilename, "E_simplex_max"));
  if (parse("-E_simplex_max=", argc, argv)) // command line overrules
  {
    parse_double("-E_simplex_max=", simplex_max[fit_E_],
      argc, argv);
    simplex_max[fit_E_] = log10(simplex_max[fit_E_]);
  }

  simplex_max[fit_n_] = 
    log10(double_from_parfile(parfilename, "n_simplex_max"));
  if (parse("-n_simplex_max=", argc, argv)) // command line overrules
  {
    parse_double("-n_simplex_max=", simplex_max[fit_n_],
      argc, argv);
    simplex_max[fit_n_] = log10(simplex_max[fit_n_]);
  }
  
  simplex_max[fit_theta_obs_] = 
    double_from_parfile(parfilename, "theta_obs_simplex_max");
  if (parse("-theta_obs_simplex_max=", argc, argv)) // command line overrules
    parse_double("-theta_obs_simplex_max=", simplex_max[fit_theta_obs_],
      argc, argv);

  simplex_max[fit_p_] = double_from_parfile(parfilename, "p_simplex_max");
  if (parse("-p_simplex_max=", argc, argv)) // command line overrules
    parse_double("-p_simplex_max=", simplex_max[fit_p_], argc, argv);

  simplex_max[fit_epsilon_B_] = 
    log10(double_from_parfile(parfilename, "epsilon_B_simplex_max"));
  if (parse("-epsilon_B_simplex_max=", argc, argv)) // command line overrules
  {
    parse_double("-epsilon_B_simplex_max=", simplex_max[fit_epsilon_B_],
      argc, argv);
    simplex_max[fit_epsilon_B_] = log10(simplex_max[fit_epsilon_B_]);
  }

  simplex_max[fit_epsilon_E_] = 
    log10(double_from_parfile(parfilename, "epsilon_E_simplex_max"));
  if (parse("-epsilon_E_simplex_max=", argc, argv)) // command line overrules
  {
    parse_double("-epsilon_E_simplex_max=", simplex_max[fit_epsilon_E_],
      argc, argv);
    simplex_max[fit_epsilon_E_] = log10(simplex_max[fit_epsilon_E_]);
  }

  simplex_max[fit_ksi_N_] = 
    log10(double_from_parfile(parfilename, "ksi_N_simplex_max"));
  if (parse("-ksi_N_simplex_max=", argc, argv)) // command line overrules
  {
    parse_double("-ksi_N_simplex_max=", simplex_max[fit_ksi_N_],
      argc, argv);
    simplex_max[fit_ksi_N_] = log10(simplex_max[fit_ksi_N_]);
  }

  //----------------------------------------------------------------------------
  // read calculation settings from parfile
  
  if (int_from_parfile(parfilename, "self_absorption") == 1) usessa = true;
    else usessa = false;
  read_status = parse_int("-self_absorption=", tempint, argc, argv);
  if (read_status)
  {
    if (tempint == 0)
      usessa = false;
    else
      usessa = true;
  }

  if (int_from_parfile(parfilename, "electron_cooling") == 1) usecooling = true;
    else usecooling = false;
  read_status = parse_int("-electron_cooling=", tempint, argc, argv);
  if (read_status)
  {
    if (tempint == 0)
      usecooling = false;
    else
      usecooling = true;
  }

  if (int_from_parfile(parfilename, "include_box") == 1) usebox = true;
    else usebox = false;
  read_status = parse_int("-include_box=", tempint, argc, argv);
  if (read_status)
  {
    if (tempint == 0)
      usebox = false;
    else
      usebox = true;
  }

  if (int_from_parfile(parfilename, "include_BM") == 1) useBM = true;
    else useBM = false;
  read_status = parse_int("-include_BM=", tempint, argc, argv);
  if (read_status)
  {
    if (tempint == 0)
      useBM = false;
    else
      useBM = true;
  }

  //----------------------------------------------------------------------------
  // simplex I/O related
  
  if (int_from_parfile(parfilename, "start_from_simplex") == 1) 
    start_from_simplex = true;
    else start_from_simplex = false;
  read_status = parse_int("-start_from_simplex=", tempint, argc, argv);
  if (read_status)
  {
    if (tempint == 0)
      start_from_simplex = false;
    else
      start_from_simplex = true;
  }

  //----------------------------------------------------------------------------
  // read resolution related settings from parfile
  
  flux.res = int_from_parfile(parfilename, "fluid_res");
  if (parse("-fluid_res=", argc, argv)) // command line overrules
    parse_int("-fluid_res=", flux.res, argc, argv);
  
  flux.eds.ur_rays = int_from_parfile(parfilename, "eds_r_res");
  if (parse("-eds_r_res=", argc, argv)) // command line overrules
    parse_int("-eds_r_res=", flux.eds.ur_rays, argc, argv);

  if (flux.eds.ur_rays % 10 != 0)
  {
    if (myid == 0)
    {
      printf("ERROR: eds_r_res should be divisible by 10.\n");
      fflush(stdout);
    }
    return 1;
  }

  flux.eds.uphi_rays = int_from_parfile(parfilename, "eds_phi_res");
  if (parse("-eds_phi_res=", argc, argv)) // command line overrules
    parse_int("-eds_phi_res=", flux.eds.uphi_rays, argc, argv);

  if (flux.eds.uphi_rays != 1 and flux.eds.uphi_rays % 10 != 0)
  {
    if (myid == 0)
    {
      printf("ERROR: eds_phi_res should be 1 or divisible by 10.\n");
      fflush(stdout);
    }
    return 1;
  }

  // Because we integrate flux using a 6 point closed interval integrator, we 
  // add one point for off-axis (n points means n plus one boundaries)
  if (flux.eds.uphi_rays != 1)
    flux.eds.uphi_rays++;

  temp = double_from_parfile(parfilename, "temp");
  read_status = parse_double("-temp=", tempdouble, argc, argv);
  if (read_status)
    temp = tempdouble;
   
  temp_factor = double_from_parfile(parfilename, "temp_factor");
  read_status = parse_double("-temp_factor=", tempdouble, argc, argv);
  if (read_status)
    temp_factor = tempdouble;
  
  temp_lowest = double_from_parfile(parfilename, "temp_lowest");
  read_status = parse_double("-temp_lowest=", tempdouble, argc, argv);
  if (read_status)
    temp_lowest = tempdouble;

  iter = int_from_parfile(parfilename, "iter");
  read_status = parse_int("-iter=", tempint, argc, argv);
  if (read_status)
    iter = tempint;

  simplex.max_iter = iter;

  MC_runs = int_from_parfile(parfilename, "MC_runs");
  read_status = parse_int("-MC_runs=", tempint, argc, argv);
  if (read_status)
    MC_runs = tempint;
  
  read_status = read_parfile_int(parfilename, "box0", tempint);
  if (read_status)
    box0 = tempint;
  if (parse("-box0=", argc, argv)) // command line overrules
    parse_int("-box0=", box0, argc, argv);

  if (box0 < 0) // default value unchanged, so nothing was provided
  {
    if (myid == 0)
    {
      printf("ERROR: no box0 value provided.\n");
      fflush(stdout);
    }
    return 1;
  }

  read_status = read_parfile_int(parfilename, "box1", tempint);
  if (read_status)
    box1 = tempint;
  if (parse("-box1=", argc, argv)) // command line overrules
  {
    parse_int("-box1=", box1, argc, argv);
  }
  if (box1 < 0) // default value unchanged, so nothing was provided
  {
    if (myid == 0)
    {
      printf("ERROR: no box1 value provided.\n");
      fflush(stdout);
    }
    return 1;
  }
  
  if (box1 < box0)
  {
    if (myid == 0)
    {
      printf("ERROR: box1 < box0: %d < %d\n", box1, box0);
      fflush(stdout);
    }
    return 1;
  }
  
  no_boxes = box1 - box0 + 1;
  if (no_boxes <= 1)
  {
    if (myid == 0)
    {
      printf("ERROR: need to load at least two BOX files.\n");
      fflush(stdout);
    }
    return 1;
  }

  no_points = int_from_parfile(parfilename, "no_points");
  read_status = parse_int("-no_points=", tempint, argc, argv);
  if (read_status)
    no_points = tempint;

  flux.lfac_initial = double_from_parfile(parfilename, "BM_start");
  read_status = parse_double("-BM_start=", tempdouble, argc, argv);
  if (read_status)
    flux.lfac_initial = tempdouble;
  
  flux.lfac_final = double_from_parfile(parfilename, "BM_stop");
  read_status = parse_double("-BM_stop=", tempdouble, argc, argv);
  if (read_status)
    flux.lfac_final = tempdouble;

  read_status = read_parfile_double(parfilename, "eds_ur_max", tempdouble);
  if (read_status) flux.ur_overrule = tempdouble;
  read_status = parse_double("-eds_ur_max=", tempdouble, argc, argv);
  if (read_status)
    flux.ur_overrule = tempdouble;

  read_status = read_parfile_double(parfilename, "t_e_0", tempdouble);
  if (read_status) flux.t0_overrule = tempdouble;
  read_status = parse_double("-t_e_0=", tempdouble, argc, argv);
  if (read_status)
    flux.t0_overrule = tempdouble;
  
  read_status = read_parfile_double(parfilename, "t_e_1", tempdouble);
  if (read_status) flux.t1_overrule = tempdouble;
  read_status = parse_double("-t_e_1=", tempdouble, argc, argv);
  if (read_status)
    flux.t1_overrule = tempdouble;
  
  //----------------------------------------------------------------------------
  // Initialize the flux class, including multibox subclass

  flux.mbox.no_boxes = no_boxes;
  flux.mbox.initialize(); // allocate array of box classes. No file I/O yet.

  //----------------------------------------------------------------------------
  // check user provided settings dealing with BOX I/O

  box_filename = new char[strlen(box_filename_base) + 10];
  box_filename_allocated = true;

  if (no_boxes <= 0) 
  { 
    printf("ERROR: no_boxes <= 0\n");
    fflush(stdout);
    return 1;
  }

  for (i = box0; i <= box1; i++)
  {
    if (box_filename_base[strlen(box_filename_base)-1] == '_')
      sprintf(box_filename, "%s%02d.h5", box_filename_base, i);
    else
      sprintf(box_filename, "%s_%02d.h5", box_filename_base, i);

    if (!file_exists(box_filename))
    {
      printf("ERROR: no_boxes set to %d, while file %s does not exist.\n", 
        no_boxes, box_filename);
      fflush(stdout);
      return 1;
    }
  }

  #if OPEN_MPI_ == ENABLED_
    MPI_Barrier(MPI_COMM_WORLD); // wait until all cores have caught up
  #endif

  //----------------------------------------------------------------------------
  // load the boxes from disc

  if (myid == 0)
  {
    printf("# Loading %d BOXES ...", no_boxes); fflush(stdout);
  }

  j = 0;
  for (i = box0; i <= box1; i++)
  {
    // prepare the box
    if (box_filename_base[strlen(box_filename_base)-1] == '_')
      sprintf(box_filename, "%s%02d.h5", box_filename_base, i);
    else
      sprintf(box_filename, "%s_%02d.h5", box_filename_base, i);
    flux.mbox.box[j].load(box_filename);
    j++;
  }

  #if OPEN_MPI_ == ENABLED_
    MPI_Barrier(MPI_COMM_WORLD); // wait until all cores have caught up
  #endif

  if (myid == 0)
  {
    printf("done. First box has version %1.2f\n", flux.mbox.box[0].version); 
    fflush(stdout);
  }
 
  // check if angles in boxes make sense
  if (myid == 0)
  {
    for (i = 0; i < no_boxes; i++)
    {
      for (j = 0; j < no_boxes; j++)
      {
        // check if no repeated initial jet opening angles
        if (fabs(flux.mbox.box[i].theta_0 - flux.mbox.box[j].theta_0) < 1e-8 and
          i != j)
        {
          printf("ERROR: box %d and box %d have the same initial opening angle "
            "theta_0 = %1.5f.\n", i, j, flux.mbox.box[i].theta_0);
          fflush(stdout);
          return 1;
        }
        
        // check if opening angles in correct order
        if (flux.mbox.box[i].theta_0 > flux.mbox.box[j].theta_0 and i <= j)
        {
          printf("ERROR: box filenames in wrong order. Box %d has angle %1.5f, "
            "while box %d has angle %1.5f. Box %d should have the larger "
            "angle.\n", i, flux.mbox.box[i].theta_0, j, 
            flux.mbox.box[j].theta_0, j);
          fflush(stdout);
          return 1;
        }
      }
    }
  }

  #if OPEN_MPI_ == ENABLED_
    MPI_Barrier(MPI_COMM_WORLD); // wait until all cores have caught up
  #endif

  //----------------------------------------------------------------------------
  // jet and counterjet settings, default setting is BOTH jets, no counterjet
  // is possible for a boosted frame BOX without explicit counter-jet data

  read_status = read_parfile_string(parfilename, "jet", temp_string);
  if (read_status == false and !parse("-box_filename_base=", argc, argv))
  {
    if (myid == 0)
    {
      printf("# WARNING: no 'jet' setting provided by user, defaulting to"
        " both forward and receding jet.\n");
      fflush(stdout);
    }
  }
  else
  {
    // overrule with command line if command line setting is provided
    if (parse("-jet=", argc, argv))
      parse_string("-jet=", temp_string, argc, argv);
  
    if (strcmp(temp_string, "forward") == 0)
    {
      flux.mbox.box[0].jet = FORWARD_;
    }
    else if (strcmp(temp_string, "receding") == 0)
      flux.mbox.box[0].jet = RECEDING_;
    else if (strcmp(temp_string, "both") != 0)
    {
      if (myid == 0)
      {
        printf("ERROR: invalid string provided for 'jet' setting: %s\n",
          temp_string);
        fflush(stdout);
      }
      return 1;
    }
    
    for (i = 1; i < no_boxes; i++)
      flux.mbox.box[i].jet = flux.mbox.box[0].jet;
  }

  // overrule jet settings if boosted frame without counter-jet data
  #if BOOST_ == ENABLED_

    if (flux.mbox.box[0].counterjet == false)
    {
      // make sure to warn user if parameter file implies that emission
      // from more than / other than the forward jet is requested
      if (flux.mbox.box[0].jet != FORWARD_ and myid == 0)
      {
        printf("# WARNING: User settings indicate something other than forward "
          "jet only is requested, but the current Lorentz-boosted BOX files"
          " only include a forward jet.\n#   Forward jet only emission will be"
          " what is computed.\n");
        printf("#   ATTEMPTED jet setting: ");
        if (flux.mbox.box[0].jet == RECEDING_)
          printf("receding\n");
        else
          printf("both\n");
        fflush(stdout);
      }

      for (i = 0; i < no_boxes; i++)
        flux.mbox.box[i].jet = FORWARD_;
    }

  #endif

  //----------------------------------------------------------------------------

  // initialize the dataset
  if (what_to_do != -1 and what_to_do != 1 and what_to_do != 2)
    if (dataset.load(data_filename) != 0) return 1;

  //----------------------------------------------------------------------------
  // set up thawed fit variable indicators
  
  thawed = 0;
  for (i = 0; i < no_fit_vars_; i++)
  {
    if (!varfrozen[i])
    {
      fitsubentry[thawed] = i;
      thawed++;
    }
  }

  //----------------------------------------------------------------------------
  // set thawed subsets of fit parameters
  
  simplex.initialize(thawed);
  
  subvarmin = array_1d<double>(thawed);
  subvarmin_allocated = true;
  
  subvarmax = array_1d<double>(thawed);
  subvarmax_allocated = true;

  for (i = 0; i < thawed; i++)
  {
    subvarmin[i] = varmin[fitsubentry[i]];
    subvarmax[i] = varmax[fitsubentry[i]];
    simplex.x_min[i] = varmin[fitsubentry[i]];
    simplex.x_max[i] = varmax[fitsubentry[i]];
  }

  //----------------------------------------------------------------------------
  // error-check the user-provided settings

  if (myid == 0)
  {
    success = check_input_parameters();
  }

  #if OPEN_MPI_ == ENABLED_
    MPI_Barrier(MPI_COMM_WORLD); // wait until all cores have caught up
    MPI_Bcast(&success, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD); // wait until all cores have caught up
  #endif

  return success;
}

////////////////////////////////////////////////////////////////////////////////

void c_boxfit :: prepare()
{
  // prepare observer information
  Obs.theta = fitvar[fit_theta_obs_];
  Obs.set_sincos();
    
  // prepare flux class
  flux.mbox.E = pow(10.0, fitvar[fit_E_]);
  flux.mbox.n = pow(10.0, fitvar[fit_n_]);
  flux.mbox.set_theta_0(fitvar[fit_theta_0_]);

  // prepare the emab
  if (usecooling) flux.emab.electron_cooling = true;
  else flux.emab.electron_cooling = false;
  if (usessa) flux.emab.absorption = true;
  else flux.emab.absorption = false;
  
  flux.emab.p_synch = fitvar[fit_p_];
  flux.emab.epsilon_B = pow(10.0, fitvar[fit_epsilon_B_]);
  flux.emab.epsilon_E = pow(10.0, fitvar[fit_epsilon_E_]);
  flux.emab.ksi_N = pow(10.0, fitvar[fit_ksi_N_]);

  // prepare the eds
  flux.eds.initialize();

  flux.eds.p_emab = &(flux.emab);
  
  flux.mbox.set_BOXBM(usebox, useBM);
}

////////////////////////////////////////////////////////////////////////////////

int c_boxfit :: check_input_parameters()
{
  // check light curve and spectrum settings
  if (what_to_do == 2) 
  {
    if (nu_0 > nu_1)
    { 
      printf("ERROR: nu_0 > nu_1\n");
      fflush(stdout);
      return 1;
    }
    
    if (nu_1 <= 0.0) 
    { 
      printf("ERROR: nu_1 < 0\n"); 
      fflush(stdout); 
      return 1;
    }
  }
  
  if (what_to_do == 2 or what_to_do == 1)
  {
    if (nu_0 <= 0.0)
    { 
      printf("ERROR: nu_0 < 0\n");
      fflush(stdout);
      return 1;
    }
    
    if (t_0 <= 0.0)
    { 
      printf("ERROR: t_0 < 0\n");
      fflush(stdout);
      return 1;
    }
    
    if (no_points <= 0) 
    { 
      printf("ERROR: no_points <= 0\n");
      fflush(stdout);
      return 1;
    }
  }
  
  if (what_to_do == 3)
  {
    if (t_0 > t_1)
    { 
      printf("ERROR: t_0 > t_1\n"); 
      fflush(stdout);
      return 1;
    }
    
    if (t_1 <= 0.0 )
    { 
      printf("ERROR: t_1 < 0\n");
      fflush(stdout);
      return 1;
    }
  }

  // check redshift and distance
  if (Obs.dL <= 0.0)
  { 
    printf("ERROR: d_L < 0.0\n");
    fflush(stdout);
    return 1;
  }
  
  if (Obs.z < 0.0)
  { 
    printf("ERROR: z < 0.0\n");
    fflush(stdout);
    return 1;
  }
  
  // check what to include
  if (!usebox and !useBM)
  { 
    printf("ERROR: excluding both box and BM\n");
    fflush(stdout);
    return 1;
  }
  
  if (usebox and useBM and flux.lfac_final > 0)
  {
    printf("# Warning. BM final Lorentz factor BM_Stop is set to %e, but will"
      " be ignored since BOX files are also used\n", flux.lfac_final);
    fflush(stdout);
  }
  
  //----------------------------------------------------------------------------
  // check resolution settings

  if (!usebox and flux.lfac_final < 1)
  {
    printf("ERROR: if no BOX data included, lfac_final has to be set >= 1. It"
      " is now set to %e\n", flux.lfac_final);
    fflush(stdout);
    return 1;
  }
  
  if (flux.ur_overrule > 0. and flux.ur_overrule < 1e10)
  {
    printf("# WARNING EDS ur_max overruled to very small value at %e\n",
      flux.ur_overrule);
  }

  if (flux.t0_overrule > 0. and flux.t1_overrule > 0. and
    flux.t0_overrule > flux.t1_overrule)
  {
    printf("ERROR: t_e_0 > t_e_1, with %e > %e\n", flux.t0_overrule, 
      flux.t1_overrule);
    fflush(stdout);
    return 1;
  }
  
  if (flux.eds.ur_rays < 1)
  { 
    printf("ERROR: eds_r_res too low: %d\n", flux.eds.ur_rays);
    fflush(stdout);
    return 1;
  }
  
  if (flux.eds.ur_rays < 100)
  {
    printf("# WARNING: eds_r_res very low: %d\n", flux.eds.ur_rays);
    fflush(stdout);
  }
  
  if (flux.eds.uphi_rays < 1)
  { 
    printf("ERROR: eds_phi_res too low: %d\n", flux.eds.uphi_rays);
    fflush(stdout);
    return 1;
  }
  
  if (flux.eds.uphi_rays == 1 and fitvar[fit_theta_obs_] > 0.0)
  {
    printf("ERROR: eds_phi_res set to 1, but blast viewed off-axis\n");
    fflush(stdout);
    return 1;
  }
  
  if (flux.eds.uphi_rays > 1 and fitvar[fit_theta_obs_] <= 0.0 and
    varfrozen[fit_theta_obs_])
  {
    printf("# WARNING: inefficient resolution settings. When fitting on-axis,"
      " use eds_phi_res = 1\n");
    fflush(stdout);
  }
  
  if (flux.res < 1)
  {
    printf("ERROR: fluid_res too low: %d\n", flux.res);
    fflush(stdout);
    return 1;
  }
  
  if (flux.res < 100)
  {
    printf("# WARNING: fluid_res very low: %d\n", flux.res);
  }

  // Simulated annealing temperature settings
  if (temp < 0.0)
  {
    printf("ERROR: unphysical annealing temperature: %e\n", temp);
    fflush(stdout);
    return 1;
  }

  if (temp_factor >= 1.0)
  {
    printf("ERROR: annealing downscale factor too high: %e\n", temp_factor);
    fflush(stdout);
    return 1;
  }
  
  if (temp_lowest > temp)
  {
    printf("ERROR: Final annealing temperature %e > Initial temperature %e.\n",
      temp_lowest, temp);
    fflush(stdout);
    return 1;
  }

  // check initial fit parameters
  if (fitvar[fit_theta_0_] < flux.mbox.box[0].theta_0)
  {
    printf("ERROR: fit value for theta_0 too small to be covered by boxes\n");
    fflush(stdout);
    return 1;
  }

  if (fitvar[fit_theta_0_] > flux.mbox.box[no_boxes-1].theta_0)
  {
    printf("ERROR: fit value for theta_0 too large to be covered by boxes\n");
    printf("       fit value is %e\n", fitvar[fit_theta_0_]);
    printf("       largest covered value is %e\n", 
      flux.mbox.box[no_boxes-1].theta_0);
    fflush(stdout);
    return 1;
  }

  if (fitvar[fit_theta_obs_] < 0.0 or fitvar[fit_theta_obs_] > 0.5*PI)
  {
    printf("ERROR: unphysical fit value observer angle: %e\n", 
      fitvar[fit_theta_obs_]);
    fflush(stdout);
    return 1;
  }
  
  if (fitvar[fit_p_] <= 2.0)
  {
    printf("ERROR: domain error for fit variable p: %e\n", fitvar[fit_p_]);
    fflush(stdout);
    return 1;
  }

  if (fitvar[fit_epsilon_E_] > log10(1.0))
  {
    printf("ERROR: fit value epsilon_E > 1.0: %e\n",
      pow(10.0, fitvar[fit_epsilon_E_]));
    fflush(stdout);
    return 1;
  }

  if (fitvar[fit_epsilon_B_] > log10(1.0))
  {
    printf("ERROR: fit value epsilon_E > 1.0: %e\n",
      pow(10.0, fitvar[fit_epsilon_B_]));
    fflush(stdout);
    return 1;
  }

  if (fitvar[fit_ksi_N_] > log10(1.0))
  {
    printf("ERROR: fit value epsilon_E > 1.0: %e\n",
      pow(10.0, fitvar[fit_ksi_N_]));
    fflush(stdout);
    return 1;
  }
  
  // check if sufficient parameters are thawed
  if (thawed == 0)
  {
    printf("ERROR: all fit variables are frozen\n");
    fflush(stdout);
    return 1;
  }
  
  // check if no degeneracy exists due to xi_N
  if (!varfrozen[fit_ksi_N_] and !varfrozen[fit_E_] and !varfrozen[fit_n_]
    and !varfrozen[fit_epsilon_E_] and !varfrozen[fit_epsilon_B_])
  {
    printf("# WARNING: degeneracy between xi_N and E, n, epsilon_E, epsilon_B");
    printf(".\n# Do not keep all these parameters thawed.\n");
    printf("# See Eichler & Waxman 2005, ApJ 627, 861 for more information.\n");
    fflush(stdout);
  }
  
  // check fit parameter ranges, for any run involving fitting
  if (what_to_do == 3 or what_to_do == 4 or what_to_do == 6)
  {
    // check fit variable minima
    if (varmin[fit_theta_0_] < flux.mbox.box[0].theta_0)
    {
      printf("ERROR: min value for theta_0 too small to be covered by boxes\n");
      fflush(stdout);
      return 1;
    }

    if (varmin[fit_theta_0_] > flux.mbox.box[no_boxes-1].theta_0)
    {
      printf("ERROR: min value for theta_0 too large to be covered by boxes\n");
      fflush(stdout);
      return 1;
    }

    if (varmin[fit_theta_obs_] < 0.0 or varmin[fit_theta_obs_] > 0.5*PI)
    {
      printf("ERROR: unphysical min value observer angle: %e\n", 
        varmin[fit_theta_obs_]);
      fflush(stdout);
      return 1;
    }

    if (varmin[fit_p_] <= 2.0)
    {
      printf("ERROR: domain error for min variable p: %e\n", varmin[fit_p_]);
      fflush(stdout);
      return 1;
    }

    if (varmin[fit_epsilon_E_] > log10(1.0))
    {
      printf("ERROR: min value epsilon_E > 1.0: %e\n",
        pow(10.0, varmin[fit_epsilon_E_]));
      fflush(stdout);
      return 1;
    }

    if (varmin[fit_epsilon_B_] > log10(1.0))
    {
      printf("ERROR: min value epsilon_E > 1.0: %e\n",
        pow(10.0, varmin[fit_epsilon_B_]));
      fflush(stdout);
      return 1;
    }

    if (varmin[fit_ksi_N_] > log10(1.0))
    {
      printf("ERROR: min value epsilon_E > 1.0: %e\n",
        pow(10.0, varmin[fit_ksi_N_]));
      fflush(stdout);
      return 1;
    }
    
    // check fit variable maxima
    if (varmax[fit_theta_0_] < flux.mbox.box[0].theta_0)
    {
      printf("ERROR: max value for theta_0 too small to be covered by boxes\n");
      fflush(stdout);
      return 1;
    }

    if (varmax[fit_theta_0_] > flux.mbox.box[no_boxes-1].theta_0)
    {
      printf("ERROR: max value for theta_0 too large to be covered by boxes\n");
      fflush(stdout);
      return 1;
    }

    if (varmax[fit_theta_obs_] < 0.0 or varmax[fit_theta_obs_] > 0.5 * PI)
    {
      printf("ERROR: unphysical max value observer angle: %e\n", 
        varmax[fit_theta_obs_]);
      fflush(stdout);
      return 1;
    }

    if (varmax[fit_p_] <= 2.0)
    {
      printf("ERROR: domain error for max variable p: %e\n", varmax[fit_p_]);
      fflush(stdout);
      return 1;
    }

    if (varmax[fit_epsilon_E_] > log10(1.0))
    {
      printf("ERROR: max value epsilon_E > 1.0: %e\n",
        pow(10.0, varmax[fit_epsilon_E_]));
      fflush(stdout);
      return 1;
    }

    if (varmax[fit_epsilon_B_] > log10(1.0))
    {
      printf("ERROR: max value epsilon_E > 1.0: %e\n",
        pow(10.0, varmax[fit_epsilon_B_]));
      fflush(stdout);
      return 1;
    }

    if (varmax[fit_ksi_N_] > log10(1.0))
    {
      printf("ERROR: max value epsilon_E > 1.0: %e\n",
        pow(10.0, varmax[fit_ksi_N_]));
      fflush(stdout);
      return 1;
    }
    
    // check whether initial simplex max within allowed range
    if (simplex_max[fit_theta_0_] < flux.mbox.box[0].theta_0 and 
      !varfrozen[fit_theta_0_])
    {
      printf("ERROR: simplex max value for theta_0 too small to be covered by "
        "boxes\n");
      fflush(stdout);
      return 1;
    }

    if (simplex_max[fit_theta_0_] > flux.mbox.box[no_boxes-1].theta_0 and
      !varfrozen[fit_theta_0_])
    {
      printf("ERROR: simplex max value for theta_0 too large to be covered by "
        "boxes\n");
      fflush(stdout);
      return 1;
    }

    if ((simplex_max[fit_theta_obs_] < 0.0 or 
      simplex_max[fit_theta_obs_] > 0.5 * PI) and !varfrozen[fit_theta_obs_])
    {
      printf("ERROR: unphysical simplex max value observer angle: %e\n", 
        simplex_max[fit_theta_obs_]);
      fflush(stdout);
      return 1;
    }

    if (simplex_max[fit_p_] <= 2.0 and !varfrozen[fit_p_])
    {
      printf("ERROR: domain error for simplex max variable p: %e\n", 
       simplex_max[fit_p_]);
      fflush(stdout);
      return 1;
    }

    if (simplex_max[fit_epsilon_E_] > log10(1.0) and !varfrozen[fit_epsilon_E_])
    {
      printf("ERROR: simplex max value epsilon_E > 1.0: %e\n",
        pow(10.0, simplex_max[fit_epsilon_E_]));
      fflush(stdout);
      return 1;
    }

    if (simplex_max[fit_epsilon_B_] > log10(1.0) and 
      !varfrozen[fit_epsilon_B_])
    {
      printf("ERROR: simplex max value epsilon_E > 1.0: %e\n",
        pow(10.0, simplex_max[fit_epsilon_B_]));
      fflush(stdout);
      return 1;
    }

    if (simplex_max[fit_ksi_N_] > log10(1.0) and 
      !varfrozen[fit_ksi_N_])
    {
      printf("ERROR: simplex max value epsilon_E > 1.0: %e\n",
        pow(10.0, simplex_max[fit_ksi_N_]));
      fflush(stdout);
      return 1;
    }
    
    // check whether initial simplex is ok
    if (!varfrozen[fit_theta_obs_] and 
      simplex_max[fit_theta_obs_] <= fitvar[fit_theta_obs_])
    {
      printf("ERROR: simplex max for theta_obs should be larger than "
        "initial fit variable value\n");
      fflush(stdout);
      return 1;
    }

    if (!varfrozen[fit_theta_0_] and 
      simplex_max[fit_theta_0_] <= fitvar[fit_theta_0_])
    {
      printf("ERROR: simplex max for theta_0 should be larger than "
        "initial fit variable value\n");
      fflush(stdout);
      return 1;
    }

    if (!varfrozen[fit_E_] and simplex_max[fit_E_] <= fitvar[fit_E_])
    {
      printf("ERROR: simplex max for E should be larger than "
        "initial fit variable value\n");
      fflush(stdout);
      return 1;
    }

    if (!varfrozen[fit_n_] and simplex_max[fit_n_] <= fitvar[fit_n_])
    {
      printf("ERROR: simplex max for n should be larger than "
        "initial fit variable value\n");
      fflush(stdout);
      return 1;
    }

    if (!varfrozen[fit_epsilon_E_] and 
      simplex_max[fit_epsilon_E_] <= fitvar[fit_epsilon_E_])
    {
      printf("ERROR: simplex max for epsilon_E should be larger than "
        "initial fit variable value\n");
      fflush(stdout);
      return 1;
    }

    if (!varfrozen[fit_epsilon_B_] and 
      simplex_max[fit_epsilon_B_] <= fitvar[fit_epsilon_B_])
    {
      printf("ERROR: simplex max for epsilon_B should be larger than "
        "initial fit variable value\n");
      fflush(stdout);
      return 1;
    }

    if (!varfrozen[fit_p_] and simplex_max[fit_p_] <= fitvar[fit_p_])
    {
      printf("ERROR: simplex max for p should be larger than "
        "initial fit variable value\n");
      fflush(stdout);
      return 1;
    }

  }

  // Further safety checks can be added here in the future
  
  return 0; // success
}

////////////////////////////////////////////////////////////////////////////////

void c_boxfit :: release_memory()
{
  int i;
  
  // close the boxes
  for (i = 0; i < no_boxes; i++)
    flux.mbox.box[i].close();  
  
  if (fitvar_allocated) 
  {
    delete_array_1d(fitvar);
    fitvar_allocated = false;
  }
  
  if (varmax_allocated)
  {
    delete_array_1d(varmax);
    varmax_allocated = false;
  }
  
  if (varmin_allocated)
  {
    delete_array_1d(varmin);
    varmin_allocated = false;
  }

  if (subvarmin_allocated)
  {
    delete_array_1d(subvarmin);
    subvarmin_allocated = false;
  }
  
  if (subvarmax_allocated)
  {
    delete_array_1d(subvarmax);
    subvarmax_allocated = false;
  }
  
  if (varfrozen_allocated)
  {
    delete_array_1d(varfrozen);
    varfrozen_allocated = false;
  }
  
  if (fitsubentry_allocated)
  {
    delete_array_1d(fitsubentry);
    fitsubentry_allocated = false;
  }
  
  if (fitvar_initial_allocated)
  {
    delete_array_1d(fitvar_initial);
    fitvar_initial_allocated = false;
  }
  
  if (initial_simplex_allocated)
  {
    delete_array_2d(initial_simplex);
    initial_simplex_allocated = false;
  }
  
  if (simplex_max_allocated)
  {
    delete_array_1d(simplex_max);
    simplex_max_allocated = false;
  }
  
  if (data_filename_allocated)
  {
    delete[] data_filename;
    data_filename_allocated = false;
  }
  
  if (box_filename_base_allocated)
  {
    delete[] box_filename_base;
    box_filename_base_allocated = false;
  }
  
  if (simplex_filename_allocated)
  {
    delete[] simplex_filename;
    simplex_filename_allocated = false;
  }
  
  if (parfilename_allocated)
  {
    delete[] parfilename;
    parfilename_allocated = false;
  }
  
  if (box_filename_allocated)
  {
    delete[] box_filename;
    box_filename_allocated = false;
  }
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
  int success;
  
  // Initialize MPI if needed
  #if OPEN_MPI_ == ENABLED_

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Barrier(MPI_COMM_WORLD); // wait until all cores have caught up

    // get node-specific information
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, myid,
      MPI_INFO_NULL, &nodecom);

    MPI_Comm_size(nodecom, &nodesize);
    MPI_Comm_rank(nodecom, &noderank);
  
  #endif // PARALELL_
  #if OPEN_MPI_ == DISABLED_

    numprocs = 1; myid = 0; noderank = 0;

  #endif // PARALELL_

  c_boxfit boxfit;
  p_boxfit = &boxfit;
  
  success = boxfit.main(argc, argv);
  
  #if OPEN_MPI_ == ENABLED_

    MPI_Barrier(MPI_COMM_WORLD); // wait until all cores have caught up
    MPI_Finalize();

  #endif // OPEN_MPI_ == ENABLED_

  return success;
}

