////////////////////////////////////////////////////////////////////////////////
//
// observer.h
//
// Created July 2, 2014 by HJvE
// Last modified July 2, 2014 by HJvE
//
// provides observer class, which contains all information on the observer
// position
//
////////////////////////////////////////////////////////////////////////////////

#ifndef OBSERVER_H_
#define OBSERVER_H_

#include <math.h>
#include "extramath.h"

class c_observer
{
  public:
  
    // observer position
    
    double dL;    // observer luminosity distance (cm)
    double z;     // observer redshift
    double theta; // observer angle (rad)
    
    double sintheta, costheta; // auxiliary functions for observer angle
    
    // observer measurement
    
    double nu;    // observer frame frequency (Hz)
    double t;     // observer frame time (s)
    
    void set_sincos(); // set sintheta, costheta based on theta value
};

// Note: an obvious extension would be to include cosmology constants and 
// routines calculating z from dL and vice versa in this class as well

#endif // OBSERVER_H_
