////////////////////////////////////////////////////////////////////////////////
//
// flux_from_box.h
//
// Created, Jun 23 2011, by HJvE
// Last Modified, July 29, 2016 by HJvE
//
////////////////////////////////////////////////////////////////////////////////

#include "environment.h"

#ifndef FLUX_FROM_BOX_H_
#define FLUX_FROM_BOX_H_

#include "box.h"
#include "fluid_special.h"
#include "flux.h"

////////////////////////////////////////////////////////////////////////////////

class c_flux_box : public c_flux
{
  public:
  
  c_flux_box();
  ~c_flux_box();
  
  // pointers to classes
  c_multibox mbox; // probes both the box files and BM solution, if needed
  
  // parameters
  double F; // flux

  // paremeters between special grid and box
  double lfac_final, lfac_initial; // FLUID Lorentz factors begin and end
  
  int res;
  
  int counter; // data point counter, useful for debugging
  bool save_emission_profile; // true if log files enabled
  bool save_image; // true if image is to be saved to disc as well.

  double tfirst, tlast; // first and last received BOX times from which flux
   // is received will be reported here
  
  double ur0_overrule; 
  double ur_overrule, t0_overrule, t1_overrule; // if > 0, these will be used
    // to set EDS ur_max and emission times. Otherwise these will be determined
    // within the set_flux routine
  
  int set_flux(); // Set flux for continuous th0. Return 0 on success
  int set_flux_th(int i_box, double &F_th); // set flux for tabulated theta_0.
    // return 0 on success. i_box is table theta_0 entry, F_th contains flux
};

#endif // FLUX_FROM_BOX_H_
