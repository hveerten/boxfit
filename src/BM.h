////////////////////////////////////////////////////////////////////////////////
//
// BM.h
//
// Created pre 18-dec-2008, HJvE
// Last modified July 29, 2016, by HJvE
//
// This contains all the math routines for the Blandford-McKee solution.
//
// reference: Blandford & McKee 1976, Physics of Fluids, 19, 1130
//
////////////////////////////////////////////////////////////////////////////////

#ifndef BM_H_
#define BM_H_

#include "physics.h"
#include "extramath.h"

////////////////////////////////////////////////////////////////////////////////

class c_BM
{
  public:

    // input parameters
    double E; // explosion energy in erg
    double n_ext; // external medium density at distance R_scale in 1/cm^3
    double k; // slope of external density distribution
    double R_scale; // initialized at 1e17 cm
    double GAMMA; // adiabatic index
    double c; // speed of light. Initially set to v_light, but can be
      // set equal to one
    
    double r, t; // radius and time
    double A; // proportionality between lfac^2 and t_e^(k-3)

    // output
    double D, lfac_shock, lfac_shock_sqrd, e_therm, p, xi, R_shock, rho;
    double n; // number density in comoving frame
    double lfac, lfac_sqrd, beta, tau, S_x;
    
    // fluid quantities at the front, same time
    double D_front, p_front, rho_front;
    double lfac_sqrd_front, lfac_front, beta_front;
    double e_therm_front, tau_front, S_x_front;
    
    double rho_ahead; // density just ahead of shock front

    //--------------------------------------------------------------------------

    c_BM(); // constructor
    
    void set_A(); // sets the proportionality between lfac^2 and t_e^(k-3),
      // which is a very useful auxiliary quantity
    void set_global_from_time( double a_t );
    void set_global_from_lfac( double a_lfac_shock );
    void set_global_from_radius( double a_R_shock );

    void set_local( double a_r ); // assumes set_global_xxxx already called !
 
  protected:
  
    void set_front(); // set variables at the shock front, called from both
      // set_global routines after they fixed lfac_shock / lab frame time
};

////////////////////////////////////////////////////////////////////////////////
// BEGIN CODE NOT PUBLIC_

class c_BM_cooling : public c_BM
{
  // synchtron paremeters
  double p_synch, frac_E, frac_N, frac_B;

  // cooling related fluid quantities
  double D_cross, rho_cross, lfac_cross, e_therm_cross, t_cross;
  double lfac_M, lfac_m, lfac_m_cross;

  void set_local_cooling(); // assumes set_local already called !
};

// END CODE NOT PUBLIC

#endif // BM_H_
