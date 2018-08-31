////////////////////////////////////////////////////////////////////////////////
//
// BM.h
//
// Created pre 18-dec-2008, HJvE
// Last modified July 29, 2016 by HJvE
//
// This contains all the math routines for the Blandford-McKee solution with
// impulsive injection of energy
//
////////////////////////////////////////////////////////////////////////////////

#include "BM.h"

////////////////////////////////////////////////////////////////////////////////

c_BM :: c_BM()
{
  // constructor, initialize default values

  R_scale = 1.0e17; // circumburst density values refer to this radius if k != 0
  c = v_light;
  GAMMA = 4.0 / 3.0; // ultra-relativistic equation of state
}

////////////////////////////////////////////////////////////////////////////////

void c_BM :: set_A()
{
  // set proportionality between lfac^2 and t_e^(k-3), useful auxiliary variable
  A = E * (17.0 - 4.0*k ) / ( 8.0 * PI*n_ext * m_p * pow( c, 5.0 - k)
    * pow( R_scale, k ) );
}

////////////////////////////////////////////////////////////////////////////////

void c_BM :: set_global_from_time( double a_t )
{
  // set shock characteristics at given time a_t

  t = a_t;
  rho_ahead = m_p * n_ext * pow( c * t / R_scale , -k );
    // initial guess for rho_ahead, assuming shock front moves at light speed
    // we use this to calculate the Lorentz factor within the BM accuracy
  lfac_shock_sqrd = dmax( E * (17.0-4.0*k) / ( 8.0 * PI * rho_ahead
    * pow( c, 2.0 ) * pow( c * t, 3.0 ) ), 1.0 );
  lfac_shock = sqrt( lfac_shock_sqrd );
  R_shock = c * t * (1.0 - 1.0 / ( 2.0 * (4.0 - k ) * lfac_shock_sqrd ) );
  rho_ahead = n_ext * m_p * pow( R_shock / R_scale, -k );
  set_front();
}

////////////////////////////////////////////////////////////////////////////////

void c_BM :: set_global_from_lfac( double a_lfac_shock )
{
  // set shock characteristics from given shock Lorentz factor

  lfac_shock = a_lfac_shock;
  lfac_shock_sqrd = lfac_shock * lfac_shock;
  t = pow(E * (17.0 - 4.0 * k) / (8.0 * PI * n_ext * m_p * pow(c, 5.0 - k)
    * pow( R_scale, k ) * lfac_shock_sqrd ), 1.0 / ( 3.0 - k ) );
  R_shock = c * t * (1.0 - 1.0 / ( 2.0 * (4.0 - k ) * lfac_shock_sqrd ) );
  rho_ahead = n_ext * m_p * pow( R_shock / R_scale, -k );
  set_front();
}

////////////////////////////////////////////////////////////////////////////////

void c_BM :: set_global_from_radius( double a_R_shock )
{
  // set shock characteristics from given shock radius

  R_shock = a_R_shock;
  rho_ahead = n_ext * m_p * pow( R_shock / R_scale, -k );
  lfac_shock_sqrd = E * ( 17.0 - 4.0 * k ) / ( 8.0 * PI * rho_ahead * c * c
    * pow( R_shock, 3.0 ) );
  lfac_shock = sqrt( lfac_shock_sqrd );
  t = R_shock / c / (1.0 - 1.0 / ( 2.0 * ( 4.0 - k ) * lfac_shock_sqrd ) );
  set_front();
}

////////////////////////////////////////////////////////////////////////////////

void c_BM :: set_front()
{
  // set conditions immediately behind the shock front

  D_front = 2.0 * rho_ahead * lfac_shock_sqrd;
  p_front = 2.0 * onethird * rho_ahead * lfac_shock_sqrd * c * c;
  e_therm_front = p_front / ( GAMMA - 1.0 );
  lfac_sqrd_front = dmax( 1.0, 0.5 * lfac_shock_sqrd );
  lfac_front = sqrt( lfac_sqrd_front );
  beta_front = sqrt( dmax( 1.0 - 1.0 / lfac_sqrd_front , 0.0 ) );
  rho_front = D_front / lfac_front;
  tau_front = lfac_sqrd_front
    * ( rho_front * c * c + e_therm_front + p_front )
    - p_front - rho_front * lfac_front * c * c;
  S_x_front = lfac_sqrd_front 
    * ( rho_front * c * c + e_therm_front + p_front )
    * beta_front;
}

////////////////////////////////////////////////////////////////////////////////

void c_BM :: set_local( double a_r )
{
  // note that this assumes that one of the set_global routines has been called
  // sets the local fluid state at given radius
  
  r = a_r;
  if ( r > R_shock ) // outside of the shock
  {
    n = n_ext * pow( r / R_scale, - k );
    D = m_p * n;
    rho = D;
    p = 1e-5 * D * c * c; // sufficiently small not to influence dynamics
    e_therm = p / ( GAMMA - 1.0 );
    lfac_sqrd = 1.0; lfac = 1.0; beta = 0.0;
    tau = e_therm;
    S_x = 0.0;
  }
  else // inside of the shock
  {
    // set self similar variable. Note that this is less accurate than eq. 27
    // from Blandford & McKee '76. This way we use the analytical error
    // ( ~ lfac_shock^-2 ) inherent in the BM solution to our advantage. We end
    // up being slightly closer to the shock front and still able to resolve
    // the hot region (of importance if cooling is enabled) once the size of
    // this region (expressed in xi) has become of the order of the error on xi.
    xi = ( 2.0 * ( 4.0 - k ) * lfac_shock_sqrd ) * ( 1.0 - r * invv_light / t );

    // lab frame density
    D = D_front * pow(xi, -(7.0 - 2.0 * k) / (4.0 - k));
    D = dmax( D, rho_ahead * 1e-6 ); // prevent underflow
    // pressure and thermal energy
    p = p_front * pow( xi, - ( 17.0 - 4.0 * k ) / ( 12.0 - 3.0 * k ) );
    p = dmax( p, rho_ahead * 1e-10 ); // prevent underflow
    e_therm = p / ( GAMMA - 1.0 );
    // Lorentz factor and velocity
    lfac_sqrd = lfac_sqrd_front / xi + 1.0; // modified BM for smooth transit
    //lfac_sqrd = lfac_sqrd_front / xi;
    lfac = sqrt ( lfac_sqrd );
    beta = sqrt( dmax( 1.0 - 1.0 / lfac_sqrd, 0.0 ) );
    // comoving density
    rho = D / lfac;
    n = rho / m_p;
    // energy momentum tensor quantities
    tau = lfac_sqrd * ( rho * c * c + e_therm + p ) - p - rho * lfac * c * c;
    S_x = lfac_sqrd * ( rho * c * c + e_therm + p ) * beta * c;
  }
}

////////////////////////////////////////////////////////////////////////////////
// BEGIN CODE NOT PUBLIC

void c_BM_cooling :: set_local_cooling()
// assumes that set_local has been called before
{
  if ( r > R_shock )
  {
    lfac_m = 1.0; lfac_M = 1.0; // i.e. no accelerated particles
  }
  else
  {
    D_cross = D * pow( xi, (10-2*k)/(4-k) );
    e_therm_cross = e_therm * pow( xi, 2.0 * (13.0 - 2.0*k ) / (12.0 - 3.0*k) );
    lfac_cross = lfac * pow( xi, (7-2*k)/(8-2*k) );
    rho_cross = D_cross / lfac_cross;
    t_cross = pow( xi, -1.0/(4.0-k) ) * t;
    
    lfac_m_cross = (2.0 - p_synch) / (1.0 - p_synch) * frac_E * e_therm_cross 
      / ( rho_cross * invm_p * frac_N * m_e * c * c );
    
    lfac_M = 2.0 * ( 19.0 - 2.0*k ) * PI * m_e * c * lfac_cross 
      / ( sigma_T * 8.0 * PI * frac_B * e_therm_cross * t_cross )
      * pow( xi, ( 25.0 - 2.0 * k ) / ( 6.0 * ( 4.0 - k ) ) )
       / ( pow( xi, ( 19.0 - 2.0*k ) / ( 3.0 * ( 4.0 - k ) ) ) - 1.0 );
       
     if ( lfac_M != lfac_M or xi == 1.0 ) lfac_M = 1e30; 
       // upper cut-off effectively Inf

    lfac_m = lfac_m_cross 
      / ( pow( xi, (13.0-2.0*k)/(6.0*(4.0-k)) ) + lfac_m_cross / lfac_M );  
  }
}

// END CODE NOT PUBLIC
