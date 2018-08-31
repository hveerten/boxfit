////////////////////////////////////////////////////////////////////////////////
//
// box.h
//
// Created Jun 15, 2011 by HJvE
// Last modified July 29, 2016 by HJvE
//
// defines BOX structures and classes
//
// Reference: van Eerten, van der Horst & MacFadyen 2012, ApJ 749, 44
//
////////////////////////////////////////////////////////////////////////////////

#include "environment.h"

#ifndef BOX_H_
#define BOX_H_

#include "arraytools.h"
#include <stdlib.h>
#include "hdf5.h"
#include "physics.h"
#include <algorithm>
#include "coords.h"
#include "fluid.h"
#include "fluid_special.h"

////////////////////////////////////////////////////////////////////////////////

// definitions for different jet types

#ifndef FORWARD_
  #define FORWARD_  0 // these settings determine whether to calculate the
  #define RECEDING_ 1 // contribution from the forward jet, the receding jet
  #define BOTH_     2 // or both.
#endif

////////////////////////////////////////////////////////////////////////////////

class c_box : public c_fluid
{
  // The BOX class contains all the information for a given theta_0 BOX
  
  public:
  
  c_box();
  ~c_box();

  c_BM BM; // Blandford-McKee analytical solution
  
  double version; // version number of BOX data

  double ***dens, ***eint, ***velr; // fluid quantities at each R, theta, time
  double ***veltheta, ***pres;
  double *theta_max; // outer theta boundaries at each time
  double **R_max; // outer radius boundaries at each theta, time
  double **R_peak; // peak positions at each theta, time
  double ***r; // lower radial cell boundary at each radius, theta, time
  double ***dr; // cell widths
  double *t; // grid snapshot times

  #if BOOST_ == ENABLED_
  
    // fluid quantities for the counterjet, if included on the grid
    double ***dens_ctr, ***eint_ctr, ***velr_ctr;
    double ***veltheta_ctr, ***pres_ctr;

    double *theta_max_ctr; // outer theta boundaries at each time
    double **R_max_ctr; // outer radius boundaries at each theta, time
    double **R_peak_ctr; // peak tau radius at each theta, time
    double ***r_ctr; // radial positions
    double ***dr_ctr; // cell widths

  #endif

  int k; // circumburst medium slope k, either 0 or 2

  // actual physics settings, i.e. those used to generate the BOX
  double lfac_initial; // initial post-shock fluid Lorentz factor in BOX file
  double E_actual; // isotropic equivalent explosion energy before scaling
  double n_actual; // circumburst medium number density before scaling
  double theta_0; // half opening angle in radians
  
  // scaled physics settings requested by user
  double n, E; // desired values (after scaling)
  double Sr; // scale factor for radius and time
  double Sn; // scale factor for n, E
  
  #if BOOST_ == ENABLED_
    double boost;
    bool counterjet; // flag whether counterjet data explicitly in BOX. This
      // does not automatically mean counterjet data is used, just that it needs
      // to be loaded from disc separately (in practice, this applies to medium
      // boost cylindrical coordinate simulations)
  #endif

  char* BOX_filename; // string containing BOX filename

  // time scaling related
  double *r_scale0, *r_scale1;
  double *scale0, *scale1;
  double theta_scale0, theta_scale1;
  int t0, t1; // box snapshot times bounding actual time requested
  
  int i_t_cur; // current time, if snapshot times are used
  
  // time and angle interpolation related
  double fract;
  double theta_max_cur; // current jet opening angle at interpolated time
  double dtheta_cur; // current theta direction cell width

  #if BOOST_ == ENABLED_
    double theta_max_cur_ctr;
    double dtheta_cur_ctr;
  #endif
  
  bool BOX_enabled; // TRUE if BOX solution also used
  bool BM_enabled;  // TRUE if BM solution also used
  int jet;          // use both jets, or jet / counterjet only

  // BOX resolution
  int thetares;
  int tres;
  int Rres;

  //----------------------------------------------------------------------------
  // class subroutines.
  
  void load(const char* filename);
  void close();

  void set_global(); // set values valid throughout BOX at given simulation time

  void set_local(s_coordinates cor);
  void set_local_empty(s_coordinates &cor, double *a_local);
  double get_var(int vartype, s_coordinates cor);

  void set_cell_coords(int ir, int ith, double &a_r, double &a_th);

  //----------------------------------------------------------------------------
  // variables and routines not accesible from the outside

  protected:
  
  bool box_file_open; // true if the data file is open
  bool box_allocated; // true if memory for box is allocated
  bool BOX_filename_allocated; // true if memory for filename allocated

  hid_t h5_fid;
  double get_var(int vartype, int i_r, int i_theta, int i_t);

  void find_cell_mix(double a_r, double a_theta, int &cr, int &ctheta);
    // MIX VERSION. UPDATE COMMENT

  #if BOOST_ == ENABLED_
    void set_cell_coords_ctr(int ir, int ith, double &a_r, double &a_th);
    double get_var_ctr(int vartype, int i_r, int i_theta, int i_t);
    void find_cell_mix_ctr(double a_r, double a_theta, int &cr, int &ctheta);
  #endif

};

////////////////////////////////////////////////////////////////////////////////

class c_multibox
{
  // the multiBOX class contains information for a series of BOXes with a range
  // of theta_0 values
  
  public:
  
    c_multibox(); // sets initialization status to false
    ~c_multibox(); // releases memory
    
    static c_box *box;
    bool initialized; // true only if memory for boxes is assigned
    int no_boxes; // How many boxes in the multibox, which each box covering
      // a single initial opening angle

    double E, n; // requested energy and circumburst number density

    int i_th0; // index of box with initial opening angle on left of th0
    double th0; // requested opening angle
  
    //--------------------------------------------------------------------------
    
    void initialize(); // allocate array of class box
    void set_theta_0(double a_th0);
    
    // enable / disable BM solution and BOX fluid data for all boxes
    void set_BOXBM(bool BOX_setting, bool BM_setting);
    
    double get_largest_radius(); // returns scaled maximum radius in the boxes
    double fractheta0; // weight factor for interpolation between boxes
};

#endif // BOX_H_
