////////////////////////////////////////////////////////////////////////////////
//
// physics.h
//
// additional generic routines and constants for physics calculations
//
// Created pre July 28, 2010 by HJvE
// Last modified January 16, 2012 by HJvE
//
////////////////////////////////////////////////////////////////////////////////

#ifndef PHYSICS_H_
#define PHYSICS_H_

#include <cmath>
#include "extramath.h"

#define v_light     2.99792458e10   // speed of light in cm / s
#define invv_light  3.335640952e-11 // inverse speed of light s / cm
#define M_sun       1.98892e33      // solar mass in g
#define R_sun       6.9599e10       // solar radius in cm
#define m_e         9.1093897e-28   // electron mass in g
#define m_p         1.6726231e-24   // proton mass in g
#define invm_e      1.097768383e27  // inverse electron mass in 1/g
#define invm_p      5.978633202e23  // inverse proton mass in 1/g
#define h_planck    6.6260755e-27   // Planck's constant in erg * s
#define h_bar       1.05457266e-27  // Planck's constant / 2 PI in erg /s
#define k_B         1.380658e-16    // Boltzmann's constant in erg / K
#define e_e         4.803e-10       // electron charge in Gaussian cgs units
#define sigma_T     6.6524e-25      // Thomson cross section free electron cm^2
#define cgs2mJy     1e26            // quantity in cgs to mJy
#define mJy2cgs     1e-26           // quantity in mJy to cgs
#define deg2rad     0.017453292     // quantity in degrees to radians
#define rad2deg     57.29577951     // quantity in radians to degrees
#define sec2day     0.000011574     // quantity in seconds to days
#define day2sec     86400.          // quantity in days to seconds
#define parsec      3.0857e18       // one parsec in cm
#define eV2Hz       2.417661e14     // quantity in electron volt to Hertz

////////////////////////////////////////////////////////////////////////////////
// BEGIN CODE NOT PUBLIC

double Planck_w(double w, double T); // Planck's function (arg in rad/s)
double Planck_nu(double nu, double T); // Planck's function 

// END CODE NOT PUBLIC

#endif // PHYSICS_H_
