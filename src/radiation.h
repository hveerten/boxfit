////////////////////////////////////////////////////////////////////////////////
//
// radiation.h
//
// Created, Jun 28 2011, by HJvE
// Last Modified, Nov 24, 2014 by HJvE
//
// Routines to provide emission and absorption coefficients for a given
// fluid state.
//
// References:
//   van Eerten & Wijers 2009, MNRAS 394, 2164
//   van Eerten et al. 2010, MNRAS, 403, 300
//   van Eerten, Zhang & MacFadyen 2010, ApJ 722, 235
//
////////////////////////////////////////////////////////////////////////////////

#include "environment.h"

#ifndef RADIATION_H_
#define RADIATION_H_

#include "observer.h"
#include "coords.h"
#include "fluid.h"
#include "physics.h"

////////////////////////////////////////////////////////////////////////////////

// definitions for different jet types
#ifndef FORWARD_
  #define FORWARD_  0 // these settings determine whether to calculate the
  #define RECEDING_ 1 // contribution from the forward jet, the receding jet
  #define BOTH_     2 // or both.
#endif

////////////////////////////////////////////////////////////////////////////////

class c_emab
{
  public:
   
    c_emab();
    ~c_emab();
   
    c_fluid *p_fluid;

    // observer related
    c_observer *p_Obs; // pointer to observer settings class

    // radiation mechanism settings
    bool electron_cooling, absorption;
    //int jet; // remove completely. All fluid info moved to fluid src code
    double overrule_theta_max; // if positive, overrule grid or box value
    double p_synch;
    double ksi_N;
    double epsilon_B;
    double epsilon_E;

    void get_emab(s_coordinates cor, double &em, double &ab);

  protected:

    // return synchrotron shape for simplified connected power law model
    double powerlawsync(double nu_c, double invnu_m);
    double powerlawsync(double nu_c, double invnu_m, double invnu_M);
};

#endif // RADIATION_H_
