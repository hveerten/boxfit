////////////////////////////////////////////////////////////////////////////////
//
// fluid.h
// 
// Created Jun 27, 2011 by HJvE
// Last modified Nov 24, 2014 by HJvE
//
// contains base class fluid, plus labels for different fluid quantities
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FLUID_H_
#define FLUID_H_

#include <string.h>
#include "coords.h"

#define rho_         0 // comoving density
#define eint_        1 // comoving internal energy density, rest mass excluded
#define v_x_         2 // GENERALIZED x coordinate
#define v_y_         3 // GENERALIZED y coordinate
#define v_           4 // lab frame velocity (absolute value)
#define lfac_        5 // Lorentz factor (lab frame)
#define D_           6 // lab frame density
#define p_           7 // pressure (comoving frame)
#define tau_         8 // lab frame energy density tau
#define N_           9 // comoving number density
#define lfacbeta_   10 // gamma * beta
#define dpdr_       11 // radial derivative of pressure
#define drhodr_     12 // radial derivative of comoving density
#define c_s_        13 // local sound speed
#define extra_1_    14 // additional advected quantity 1

#define fluid_vars_ 15 // number of fluid variables
#define var_string_max_length_ 12 // longest variable name string including \n 

////////////////////////////////////////////////////////////////////////////////

void sprintvar(char* varstring, int varno);
  // prints the associated string into varstring, for a given varno
  
int sscanvar(char* varstring);
  // returns the var number for a given varstring. Returns -1 if no match found

class c_fluid
{
  public:
  
  virtual ~c_fluid();
  
  double local[fluid_vars_]; // array containing local fluid values
  double t_e;
  int pos_flag; // position flag, set to nonzero value for coordinate queries
    // landing on special coordinates (e.g. out of grid boundaries), depending
    // on choices made by child class

  // virtual function to set global quantities for given time
  virtual void set_global() = 0;
  // virtual function to set all fluid variables at once  
  virtual void set_local(s_coordinates cor) = 0;
  // virtual function to query a single fluid variable.
  virtual double get_var(int vartype, s_coordinates cor) = 0;
};

#endif // FLUID_H_
