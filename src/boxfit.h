////////////////////////////////////////////////////////////////////////////////
//
// boxfit.h
//
// Created Jun 28 2011, by HJvE
// Last Modified July 29, 2016 by HJvE
//
// Reference: van Eerten, van der Horst & MacFadyen 2012, ApJ 749, 44
//
////////////////////////////////////////////////////////////////////////////////

#include "environment.h"

#ifndef BOX_FIT_H_
#define BOX_FIT_H_

#include "extramath.h"
#include "observer.h"
#include "flux_from_box.h"
#include "numerical.h"
#include "param_IO.h"
#include "parse.h"

#if OPEN_MPI_ == ENABLED_
  #include "mpi.h"
#endif // OPEN_MPI_ == ENABLED_

////////////////////////////////////////////////////////////////////////////////
// define labels of fit variables

// main set of fit variables, full model fits
#define fit_theta_0_   0 // original jet half opening angle (radians)
#define fit_E_         1 // isotropic equivalent energy (Erg)
#define fit_n_         2 // circumburst particle number density (cm^-3)
#define fit_theta_obs_ 3 // observer angle (radians)
#define fit_p_         4 // synchrotron slope p
#define fit_epsilon_B_ 5 // magnetic field energy fraction
#define fit_epsilon_E_ 6 // accelerated particle energy fraction
#define fit_ksi_N_     7 // accelerated particle number density fraction

// second set of fit variables, power law fits for self-test
#define fit_S_         0 // scale factor single or broken power law fit
#define fit_alpha_1_   1 // pre-break slope
#define fit_t_break_   2 // break time
#define fit_alpha_2_   3 // post-break slope

#define no_fit_vars_   8 // total number of fit variables

////////////////////////////////////////////////////////////////////////////////
// Define wrapper functions and pointers

double funk_wrapper(double fv[]); // external wrapper for fit function
double selftest_funk_wrapper(double fv[]);

////////////////////////////////////////////////////////////////////////////////

class c_dataset 
{
  public:
  
  int n; // number of data points
  bool initialized; // if true, then memory as been allocated
  double * F; // flux as read from the dataset file
  double * Falt; // alternative fluxes, calculated for Monte Carlo procedure
  double * Fbox; // flux as calculated from box
  double * sigma; // one sigma error on the flux
  double * t; // observer time in (s)
  double * nu; // observer frequency in Hertz
  
  // Monte Carlo error determination related
  double *logF; // logarithm of best fit fluxes,
  double **dlogFdx; // partial derivatives
  double ferr; // average fractional error of data
  
  c_dataset();
  ~c_dataset();
  void initialize(); // allocate memory
  void free_memory(); // free memory
  int load(const char * filename); // load dataset from disc, return 0 if 
    // call successful
  int print(FILE *p_file); // print all datapoints, return 0 if call successfull
};

////////////////////////////////////////////////////////////////////////////////

class c_boxfit
{
  public:
  
    c_boxfit();
    ~c_boxfit();

    int main(int argc, char* argv[]);
    int selftest(); // return 0 if successful

    double funk(double fv[]); // internal wrapper for fit function
    double bkn_power(double fv[]); // wrapper for broken power law self-test

  //----------------------------------------------------------------------------
  
  protected:
  
    c_dataset dataset; // class containing all obervational and synthetic data
    c_flux_box flux; // class containing flux calculation driver routine
    c_random random; // class containing random number generation functions
    c_simplex simplex; // class containing simplex functions 

    int what_to_do; // what to calculate: light curve, iterative fit, etc.

    // light curve and spectrum settings & variables
    c_observer Obs; // class containing observer position, measurement time and
      // frequency
    double nu_0, nu_1; // starting and stopping frequency (Hz)
    double t_0, t_1; // starting and stopping time (s)
    int no_points; // how many data points in the light curve / spectrum
    double *fluxes; // array containing the results from light curve / spectrum

    // file I/O variables
    char* parfilename; // name of parameter file
    char *data_filename; // name of data file
    char *box_filename_base; // path to box files
    char *box_filename; // full box filename, used when loading boxes
    char *simplex_filename; // name of simplex file in case start from simplex
    char *temp_string; // all-purpose string storage variable for disc I/O
    bool save_intermediate; // switch save after each temp cycle
    bool start_from_simplex; // if true, read starting simplex from disc
    bool quietfit; // prevents too much output during MC runs

    // radiation calculation switches
    bool usealt; // specify whether to use F (if false) or Falt (if true)
    bool fromderivatives; // specify whether to use quick partial der. method
    bool usessa; // switch synchrotron self-absorption
    bool usecooling; // switch electron cooling
    bool usebox; // switch use box based fluid states
    bool useBM; // switch use Blandford-McKee based fluid states
    bool radial_interpolation; // switch whether to use radial interpolation
      // within a BOX to obtain local fluid variables

    // fit variable settings    
    double *fitvar; // current fit values, modified after each fit
    double *fitvar_initial; // separately store the initial fit values
    double *varmax; // upper range for each fit value
    double *varmin; // lower range for each fit value
    double *subvarmax, *subvarmin; // subsets of ranges, thawed variables only
    bool *varfrozen; // which fit values are kept fixed during fit procedure
    int *fitsubentry; // which subset entries are linked to which fitvars
    int thawed; // how many active fitvars
    double chi_sqr; // chi^2

    // Monte Carlo results
    double** MCfitvar; // stores the results for the different MCruns
    double* MCchisqr; // stores the separate MC chi^2 results
    double* MCvarmax; // stores the upper limit on MC results for each variable
    double* MCvarmin; // stores the lower limit on MC results for each variable
    
    // resolution settings
    double temp; // simulated annealing start temperature
    double temp_factor; // fractional drop in temperature after each cycle
    double temp_lowest; // lowest temperature before drop to zero temperature
    double *simplex_max; // range of initial simplex
    int iter; // number of iterations at single temperature
    int MC_runs; // number of light curves calculated during Monte Carlo run
    int no_boxes; // number of box data files (running 0 .. no_boxes - 1)
    int MCi; // current MC iteration / number of currently stored MC results
    double BM_start, BM_stop; // Blandford-McKee starting and stopping Lorentz
      // factors as set by the parfile.

    int lccounter; // counts how often light curve is calculated

    // NR based simulated annealing method variables
    double** initial_simplex; // initial simplex

    // flags tracking memory allocation
    bool fitvar_allocated;
    bool varmax_allocated;
    bool varmin_allocated;
    bool subvarmin_allocated;
    bool subvarmax_allocated;
    bool varfrozen_allocated;
    bool fitsubentry_allocated;
    bool fitvar_initial_allocated;
    bool initial_simplex_allocated;
    bool simplex_max_allocated;

    bool data_filename_allocated;
    bool box_filename_base_allocated;
    bool simplex_filename_allocated;
    bool parfilename_allocated;
    
    bool box_filename_allocated;

    //--------------------------------------------------------------------------
    
    // general initialization and shutdown related
    int initialize(int argc, char* argv[]); // general initialization, returns
      // 0 if succesfull.
    void prepare(); // prepare box, eds etc. for new values of fit parameters
    int check_input_parameters(); // check if boxfitsettings are consistent,
      // returns 0 if successful.
    double get_average_error(); // returns average fractional error of dataset
    void release_memory(); // deallocate memory
    
    // single calculation
    void initialize_simplex(); // set up simplex of initial fit attempts
    void calculate(); // calculate the box values
    void calculate_from_derivatives(); // calculate box values using part. der.
    double set_chisqr(); // assuming datapoints are updated, calculates chisqr
    
    // light curve and spectra related
    void spectrum(); // calculate a spectrum
    void lightcurve(); // calculate a light curve

    // single and multiple fit routines
    int fit_dataset(); // perform a dataset fit, return 0 if successful
    int montecarlo(); // perform a number of fits using Monte Carlo approach,
      // return 0 if successful
    int setup_derivatives(); // set up partial derivatives at best fit values,
      // return 0 if successful
    void store_new_MC(); // store latest MC iteration results in ordered list

    // printing and I/O
    void print_settings(FILE *p_where); // print settings from parameter file
    void print_opening_message(FILE *p_where); // print program name, version
    void print_fluxes(FILE *p_where); // print light curve or spectrum
    void print_fit_result(FILE *p_where); // prints fit vars and chi^2 value
    FILE *p_chi_sqr; // if not NULL print chi^2 result here whenever chi^2 calc.
    void save_MC_result();
    void save_derivatives();
    void load_derivatives();
};

extern c_boxfit *p_boxfit;

#endif // BOX_FIT_H_
