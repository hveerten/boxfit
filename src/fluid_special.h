////////////////////////////////////////////////////////////////////////////////
//
// fluid_BM.h
//
// Created June 29, 2011 by HJvE
// Last modified September 8, 2011 by HJvE
//
// Routines to interact with analytical solutions for the fluid dynamics
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FLUID_BM_H_
#define FLUID_BM_H_

#include "environment.h"
#include "physics.h"
#include "BM.h"
#include "fluid.h"
#include "coords.h"
#include "param_IO.h"
#include "parse.h"

////////////////////////////////////////////////////////////////////////////////

class c_fluid_special : public c_fluid
{
  public:

    c_BM BM;
    
    double lfac_final, lfac_initial; // FLUID Lorentz factors begin and end
    #if BOOST_ == ENABLED_
      double boost, boostsqr, beta_sim; 
    #endif // BOOST_
    
    // initial and final times are derived from lfac's
    double t_e_final; // final BM time (but usually first snapshots time)
    double t_e_initial;
//    double t_e; // current emission time. 

    double theta_max;
  
    c_fluid_special(); // constructor, set boost to 1.0
    void load_settings(char *parfilename, int argc, char* argv[]);
    void print_log(FILE *p_logfile);
    void init_times(double toverz);
    void set_global();
    void set_local(s_coordinates cor);
    double get_var(int vartype, s_coordinates cor);
};

#endif // FLUID_BM_H_
