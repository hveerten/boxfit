#ifndef _EDS_2D_REGULAR_H_
#define _EDS_2D_REGULAR_H_

////////////////////////////////////////////////////////////////////////////////
//
// EDS_2D_REGULAR
//
// Created ~ May 2011 by HJvE
// Last Modified January 12, 2012 by HJvE
//
// EDS (EquiDistant Surface)
//
////////////////////////////////////////////////////////////////////////////////

#include "environment.h"
#include "extramath.h"
#include "arraytools.h"
#include "radiation.h"
#include "observer.h"
#include "hdf5.h"

////////////////////////////////////////////////////////////////////////////////

struct s_ray
{
  double em, ab; // local emission and absorption coefficients

  #if LOCAL_SELF_ABSORPTION_ == DISABLED_
  
    double emdr, abdr; // cumulative emission and absorption coefficients
    double emdr_lores, abdr_lores; // low time resolution versions of same
  
  #else
  
    double I, k1, k2, k3, k4; // intensity I and intermediate Runge-Kutta I's
    double I_lores; // low time resolution I

  #endif
  
  double ur, uphi; // cell coordinates on EDS plane, polar
  double ux, uy;   // cell coordinates on EDS plane, Cartesian
};

////////////////////////////////////////////////////////////////////////////////

class c_eds
{
  // the EquiDistant Surface (EDS) lies parallel to the y-axis, since the
  // observer is located in the x-z plane (the jet is always along
  // the z-axis).

  public:
  c_eds();
  virtual ~c_eds();
  
  // location of the origin of the EDS, with respect to the grid origin.
  // All vectors are in grid frame Cartesian coordinates.
  double O[3]; // position of origin on grid in Cartesian coordinates
  // EDS origin distance to fluid grid origin
  double O_r;
  
  #if BOOST_ == ENABLED_
    double boost, boostsqr, beta_sim;
  #endif // BOOST

  double Dr; // travel distance between snapshots, i.e. ray stepsize
  double t_sim;

  // additional quantities for log purpose
  double F; // flux
  
  // observer information
  c_observer *p_Obs;
  
  int ur_rays; // number of rays logarithmically space in r direction
  int uphi_rays; // number of rays linearly spaced in phi direction
  bool memory_assigned; // false initially. Memory is assigned in initialize()

  double ur_max; // upper r boundary of the eds. 'u' denotes coord in eds frame
  double ur_min; // lower boundary of eds. 
  
  // related to shape of blast wave
  double R_100, R_99, R_95, R_75, R_50;
    
  //double ulnrImax; // max radius with nonzero intensity. Helps to determine
    // plot size for spatially resolved plots.

  c_emab *p_emab;

  static s_ray** ray; // 2D array of rays

  // unit vectors on the EDS
  double eux[3], euy[3]; // if O_theta is zero these are the same as x and y

  // core functions
  void initialize(); // assign memory for rays
  
  #if LOCAL_SELF_ABSORPTION_ == DISABLED_
  
    void prepare_update(); // calculates emissivities and absorptions but 
      // does not yet calculate new integrated values
    void finalize_update(double dt_sim);  // updates the eds
    void finalize_lores_update(double dt_sim);  // low times resolution update
  
  #else

    void prepare_update(int RK, double dt_sim);
    void finalize_RK_update();
  
  #endif

  // Disc IO related function
  void save_image(const char * filename);

  // output functions  
  double get_F_annulus(int iur);
    // returns the total flux in the annulus at EDS iur'th radial position,
    // for the full 0 .. 2PI angle domain of u_phi
  double get_F_annulus_lores_phi(int iur); // low resolution in uphi
  double get_F_annulus_lores_t(int iur); // low resolution in time
  double get_total_flux(); 
    // returns total flux, observed at distance R.
  double get_F_phi_error();
    // returns measure of error due to angular resolution
  double get_F_r_error();
    // returns measure of error due to radial resolution
  double get_F_t_error();
    // returns measure of error due to time step
  void set_R(); // set 1D measures of object size
  void reset(); // clears all values
  void set_coordinates(double a_t_sim);

};

#endif /* EDS_2D_REGULAR_H_ */
