////////////////////////////////////////////////////////////////////////////////
//
// physics.cpp
//
// Created pre July 28, 2010 by HJvE
// Last modified September 8, 2010 by HJvE
//
////////////////////////////////////////////////////////////////////////////////

#include "physics.h"

////////////////////////////////////////////////////////////////////////////////
// BEGIN CODE NOT PUBLIC

double Planck_w(double w, double T)
{
  double y = h_bar * w / (k_B * T);
  // the double datatype resolves only about 16 digits difference between
  // two numbers, so if y<<1 then we take the 1st order Taylor expansion
  // to prevent 1.00...0000x - 1.000....0000 -> 0
  if (y < 1e-13) return h_planck / h_planck / (v_light * v_light) * 2 * 
    pow(inv2PI, 4) * pow(w,3) / y;
    else return h_planck / (v_light * v_light) * 2 * pow(inv2PI, 4) * pow(w, 3)
      / (exp(h_bar * w / (k_B * T)) - 1.0);
} 

////////////////////////////////////////////////////////////////////////////////

double Planck_nu(double nu, double T)
{
  double y=h_planck*nu/(k_B*T);
  // the double datatype resolves only about 16 digits difference between
  // two numbers, so if y<<1 then we take the 1st order Taylor expansion
  // to prevent 1.00...0000x - 1.000....0000 -> 0
  if (y < 1e-13) return 2 * nu * nu * invv_light * invv_light * k_B * T; 
    // Rayleigh-Jeans
    else return h_planck / (v_light * v_light) * 2 * pow(inv2PI, 4) * pow(nu, 3)
      / (exp(h_planck * nu / (k_B * T)) - 1.0);
}

// END CODE NOT PUBLIC
