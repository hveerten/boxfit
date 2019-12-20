////////////////////////////////////////////////////////////////////////////////
//
// flux.h
//
// Created, July 2 2014, by HJvE
// Last Modified, July 2, 2014 by HJvE
//
// - Provides base class structure for flux calculation routines. This is then
// inherited by flux_from_box, flux_from_grid etc, which include the additional
// layer with specifics for given sources of fluid information.
//
////////////////////////////////////////////////////////////////////////////////

#include "environment.h"

#ifndef FLUX_H_
#define FLUX_H_

#include "observer.h"
#include "eds_2D_regular.h"
#include "radiation.h"
#include "fluid.h"

////////////////////////////////////////////////////////////////////////////////

class c_flux
{
  public:

    virtual ~c_flux(); // class has virtual functions, so needs virtual
      // destructor. This is a safeguard in case an inherited class based on
      // c_flux is defined through a pointer not to that inherited class
      // but to the base class. If such an instance of the inherited class
      // gets deleted via a delete[] statement acting on the pointer to that
      // instance, the delete[] command would in principle call the base
      // destructor, not the inherited class destructor, which could potentially
      // lead to memory leaks. The 'virtual' statement ensures that the 
      // inherited class destructor is not ignored.

    // classes used by the flux class

    c_eds eds;
    c_emab emab;
    c_observer *p_Obs;

    double F; // monochromatic flux value in cgs units (erg cm^-2 s^-1 Hz^-1)
    
    void set_observer_pointer(c_observer *p_Obs_arg);
    virtual int set_flux() = 0; // specific implementation for each inheriting
      // class
};

#endif // FLUX_H_
