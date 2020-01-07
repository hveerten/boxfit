////////////////////////////////////////////////////////////////////////////////
// eds_2D.h
//
// Created November 17, 2010 by HJvE
// Last Modified January 16, 2012 by HJvE
//
////////////////////////////////////////////////////////////////////////////////

#include "eds_2D_regular.h"

////////////////////////////////////////////////////////////////////////////////
// Static member definitions

s_ray** c_eds :: ray;

////////////////////////////////////////////////////////////////////////////////
// constructor en destructor
////////////////////////////////////////////////////////////////////////////////

c_eds :: c_eds()
{
  // set up initialization values for parameters
  ur_rays = 750;
  uphi_rays = 31;
  memory_assigned = false;
  #if BOOST_ == ENABLED_
    boost = 1.0, boostsqr = 1.0, beta_sim = 0.0; // start with nonmoving grid
  #endif // BOOST_
}

////////////////////////////////////////////////////////////////////////////////

void c_eds :: initialize()
{
  // assign memory
    
  if (memory_assigned) delete_array_2d(ray);
  ray = array_2d<s_ray>(ur_rays, uphi_rays);
  memory_assigned = true;
}

////////////////////////////////////////////////////////////////////////////////

c_eds :: ~c_eds()
{
  if (memory_assigned) delete_array_2d(ray);
}

////////////////////////////////////////////////////////////////////////////////
// public functions
////////////////////////////////////////////////////////////////////////////////

// Prepare update routine for global approach to synchrotron self-absorption

#if LOCAL_SELF_ABSORPTION_ == DISABLED_

  void c_eds :: prepare_update()
  {
    // set local emission and absorption coefficients for all rays
  
    double r2;
    s_coordinates cor;
    #if BOOST_ == ENABLED_
      double Ex, Ey, Ez; // projections on EDS through origin
    #endif // BOOST_
    int iur, iuphi;

    // loop over all rays
    for (iur = 0; iur < ur_rays; iur++)
    {
      for (iuphi = 0; iuphi < uphi_rays; iuphi++)
      {
        #if BOOST_ == DISABLED_

          // get coordinates in grid frame
          cor.x = O[0] + ray[iur][iuphi].ux * eux[0] + 
            ray[iur][iuphi].uy * euy[0];
          cor.y = O[1] + ray[iur][iuphi].ux * eux[1] + 
            ray[iur][iuphi].uy * euy[1];
          cor.z = O[2] + ray[iur][iuphi].ux * eux[2] + 
            ray[iur][iuphi].uy * euy[2];
          cor.t = t_sim;
    
        #else // BOOST_

          // if boost, the EDS is curved
          Ex = ray[iur][iuphi].ux * eux[0] + ray[iur][iuphi].uy * euy[0];
          Ey = ray[iur][iuphi].ux * eux[1] + ray[iur][iuphi].uy * euy[1];
          Ez = ray[iur][iuphi].ux * eux[2] + ray[iur][iuphi].uy * euy[2];

          cor.z = (Ez + p_Obs->costheta * v_light * 
            (t_sim / boost - p_Obs->t / (p_Obs->z + 1.))) / 
            (1.0 - beta_sim * p_Obs->costheta);

          cor.t = t_sim / boost + beta_sim * cor.z * invv_light;

          cor.x = Ex + (cor.t - p_Obs->t / (p_Obs->z + 1.)) * 
            p_Obs->sintheta * v_light;

          cor.y = Ey;

        #endif // BOOST_
  
        // get emmission and absorption coefficients
        r2 = cor.x * cor.x + cor.y * cor.y + cor.z * cor.z;
        if (r2 >= v_light * v_light * cor.t * cor.t)
        {
          ray[iur][iuphi].em = 0.; ray[iur][iuphi].ab = 0.;
        }
        else
        {
          // get emission and absorption coefficients
          p_emab->get_emab(cor, ray[iur][iuphi].em, ray[iur][iuphi].ab);
        }
      }
    }
  }

#endif // LOCAL_SELF_ABSORPTION == DISABLED_

////////////////////////////////////////////////////////////////////////////////
// update routine for local approach to synchrotron self-absorption

#if LOCAL_SELF_ABSORPTION_ == ENABLED_

  void c_eds :: prepare_update(int RK, double dt_sim)
  {
    int iur, iuphi;
    double tau, S, expo;
    s_coordinates cor;
    #if BOOST_ == ENABLED_
      double Ex, Ey, Ez; // projections on EDS through origin
    #endif // BOOST_

    double r2;

    #if BOOST_ == ENABLED_
      Dr = dt_sim * v_light / (boost * (1.0 - beta_sim * p_Obs->costheta));
    #endif // BOOST_
    #if BOOST_ == DISABLED_
      Dr = dt_sim * v_light;
    #endif // BOOST_

    //--------------------------------------------------------------------------
    // loop over all rays to determine em and ab

    if (RK != 1 and RK != 3)
    { 
      for (iur = 0; iur < ur_rays; iur++)
      {
        for (iuphi = 0; iuphi < uphi_rays; iuphi++)
        {
          #if BOOST_ == DISABLED_

            // get coordinates in grid frame
            cor.x = O[0] + ray[iur][iuphi].ux * eux[0] + 
              ray[iur][iuphi].uy * euy[0];
            cor.y = O[1] + ray[iur][iuphi].ux * eux[1] + 
              ray[iur][iuphi].uy * euy[1];
            cor.z = O[2] + ray[iur][iuphi].ux * eux[2] + 
              ray[iur][iuphi].uy * euy[2];
            cor.t = t_sim;
        
          #else
  
            // if boost, the EDS is curved for off-axis observers
            Ex = ray[iur][iuphi].ux * eux[0] + ray[iur][iuphi].uy * euy[0];
            Ey = ray[iur][iuphi].ux * eux[1] + ray[iur][iuphi].uy * euy[1];
            Ez = ray[iur][iuphi].ux * eux[2] + ray[iur][iuphi].uy * euy[2];

            cor.z = (Ez + p_Obs->costheta * v_light * 
              (t_sim / boost - p_Obs->t / (p_Obs->z + 1.))) / 
              (1.0 - beta_sim * p_Obs->costheta);

            cor.t = t_sim / boost + beta_sim * cor.z * invv_light;

            cor.x = Ex + (cor.t - p_Obs->t / (p_Obs->z + 1.)) * 
              p_Obs->sintheta * v_light;

            cor.y = Ey;

          #endif // BOOST_

          r2 = cor.x * cor.x + cor.y * cor.y + cor.z * cor.z;
          if (r2 >= v_light * v_light * cor.t * cor.t)
          {
            ray[iur][iuphi].em = 0.; ray[iur][iuphi].ab = 0.;
          }
          else
          {
            // get emission and absorption coefficients
            p_emab->get_emab(cor, ray[iur][iuphi].em, ray[iur][iuphi].ab);
          }
        }  // loop over iuphi
      }  // loop over iur
    } // RK != 1 and RK != 3

    //--------------------------------------------------------------------------
    // apply em and ab
  
    for (iur = 0; iur < ur_rays; iur++)
    {
      for (iuphi = 0; iuphi < uphi_rays; iuphi++)
      {
        tau = ray[iur][iuphi].ab * Dr;
 
        // store in the approprate Runge-Kutta bin
        switch(RK)
        {
          case 1:
            if (tau > 1e-3)
            {
              expo = exp(-tau);
              S = ray[iur][iuphi].em / (ray[iur][iuphi].ab + 1e-200);
              ray[iur][iuphi].k1 = ray[iur][iuphi].I * (expo - 1.) + 
              S * (1.0 - expo);
            }
            else
            {
              ray[iur][iuphi].k1 = ray[iur][iuphi].I * tau * (.5 * tau - 1.) +
                ray[iur][iuphi].em * Dr * (1. - 0.5 * tau);
            }
          break;
          case 2:
            if (tau > 1e-3)
            {
              expo = exp(-tau);
              S = ray[iur][iuphi].em / (ray[iur][iuphi].ab + 1e-200);
              ray[iur][iuphi].k2 = (ray[iur][iuphi].I + 
                0.5 * ray[iur][iuphi].k1) * (expo - 1.) + S * (1. - expo);
            }
            else
            {
              ray[iur][iuphi].k2 = (ray[iur][iuphi].I + 
                0.5 * ray[iur][iuphi].k1) * tau * (0.5 * tau - 1.) +
                ray[iur][iuphi].em * Dr * (1. - 0.5 * tau);
            }
          break;
          case 3:
            if (tau > 1e-3)
            {
              expo = exp(-tau);
              S = ray[iur][iuphi].em / (ray[iur][iuphi].ab + 1e-200);
              ray[iur][iuphi].k3 = (ray[iur][iuphi].I + 
                0.5 * ray[iur][iuphi].k2) * (expo - 1.) + S * (1. - expo);
            }
            else
            {
              ray[iur][iuphi].k3 = (ray[iur][iuphi].I + 
                0.5 * ray[iur][iuphi].k2) * tau * (0.5 * tau - 1.) +
                ray[iur][iuphi].em * Dr * (1. - 0.5 * tau);
            }
          break;
          case 4:
            if (tau > 1e-3)
            {
              expo = exp(-tau);
              S = ray[iur][iuphi].em / (ray[iur][iuphi].ab + 1e-200);
              ray[iur][iuphi].k4 = (ray[iur][iuphi].I + 
                ray[iur][iuphi].k3) * (expo - 1.) + S * (1. - expo);
            }
            else
            {
              ray[iur][iuphi].k4 = (ray[iur][iuphi].I + 
                ray[iur][iuphi].k3) * tau * (0.5 * tau - 1.) + 
                ray[iur][iuphi].em * Dr * (1. - 0.5 * tau);
            }
          break;
        } // switch
      } // iuphi loop
    } // iur loop
  }

#endif // LOCAL_SELF_ABSORPTION_ == ENABLED_

////////////////////////////////////////////////////////////////////////////////

#if LOCAL_SELF_ABSORPTION_ == DISABLED_

  void c_eds :: finalize_update(double dt_sim) 
  { 
    int iur, iuphi;

    #if BOOST_ == ENABLED_
      Dr = dt_sim * v_light / (boost * (1.0 - beta_sim * p_Obs->costheta));
    #endif // BOOST_
    #if BOOST_ == DISABLED_
    Dr = dt_sim * v_light;
    #endif // BOOST_

    // loop over all rays
    for (iur = 0; iur < ur_rays; iur++)
    {
      for (iuphi = 0; iuphi < uphi_rays; iuphi++)
      {
        ray[iur][iuphi].emdr += ray[iur][iuphi].em * Dr;
        ray[iur][iuphi].abdr += ray[iur][iuphi].ab * Dr;
      }
    }
  }

#endif // LOCAL_SELF_ABSORPTION_ == DISABLED_

////////////////////////////////////////////////////////////////////////////////

#if LOCAL_SELF_ABSORPTION_ == DISABLED_

  void c_eds :: finalize_lores_update(double dt_sim)  
  {
    int iur, iuphi;

    #if BOOST_ == ENABLED_
      Dr = dt_sim * v_light / (boost * (1.0 - beta_sim * p_Obs->costheta));
    #endif // BOOST_
    #if BOOST_ == DISABLED_
      Dr = dt_sim * v_light;
    #endif // BOOST_

    // loop over all rays
    for (iur = 0; iur < ur_rays; iur++)
    {
      for (iuphi = 0; iuphi < uphi_rays; iuphi++)
      {
        ray[iur][iuphi].emdr_lores += ray[iur][iuphi].em * Dr;
        ray[iur][iuphi].abdr_lores += ray[iur][iuphi].ab * Dr;
      }
    }
  }

#endif // LOCAL_SELF_ABSORPTION_ == DISABLED_

////////////////////////////////////////////////////////////////////////////////

#if LOCAL_SELF_ABSORPTION_ == ENABLED_

  void c_eds :: finalize_RK_update()
  {
    int iur, iuphi;
    double Iold;

    for (iur = 0; iur < ur_rays; iur++)
    {
      for (iuphi = 0; iuphi < uphi_rays; iuphi++)
      {
        Iold = ray[iur][iuphi].I;
        ray[iur][iuphi].I = fmax((Iold + 
          ray[iur][iuphi].k1 / 6. +
          ray[iur][iuphi].k2 / 3. +
          ray[iur][iuphi].k3 / 3. +
          ray[iur][iuphi].k4 / 6.), 0.);

        // traced low resolution as well, for determining error due to time step
        Iold = ray[iur][iuphi].I_lores;
        ray[iur][iuphi].I_lores = fmax(Iold + ray[iur][iuphi].k1, 0.);
      } // iuphi
    } // iur
  }

#endif // LOCAL_SELF_ABSORPTION_ == ENABLED_

////////////////////////////////////////////////////////////////////////////////

void c_eds :: set_coordinates(double a_t_sim)
{
  // The coordinates of the EDS are as follows. In 1D the normal to the eds
  // lies along the z-axis of the grid. The x-direction on the EDS lies along
  // the grid x-axis, the y on the EDS along the grid y-axis. If the EDS is
  // tilted (i.e. the observer angle is nonzero), the normal to the eds surface
  // moves in the x-z plane. As a consequence the y direction on the EDS remains
  // unaffected. The x direction on the EDS no longer coincides with the x
  // direction on the grid.

  t_sim = a_t_sim;

  #if BOOST_ == DISABLED_
    double r = (t_sim - p_Obs->t / (p_Obs->z + 1.)) * v_light;
    O[0] = r * p_Obs->sintheta; O[1] = 0.0; O[2] = r * p_Obs->costheta;
    O_r = r; // O_r < 0 if grid origin between eds and observer
  #endif
  
  // set unit vectors on the EDS plane in the frame of the grid
  eux[0] = p_Obs->costheta; eux[1] = 0.0; eux[2] = -p_Obs->sintheta;
  euy[0] = 0.0; euy[1] = 1.0; euy[2] = 0.0;
}

////////////////////////////////////////////////////////////////////////////////

void c_eds :: save_image(const char *filename)
{
  hid_t h5_fid;
  hid_t dataset, dataspace;
  hsize_t dim_2d[2], dim_1d;
  double value; // used for storing single double values in hdf5 file
  int int_value; // used for storing integer values in hdf5 file
  double **I; // local version of I
  double *ur;
  double *uphi;
  int i_ur, i_uphi;
  
  #if LOCAL_SELF_ABSORPTION_ == DISABLED_
    double tau; // cumulative optical depth
    double S; // source function
  #endif

  I = array_2d<double>(ur_rays, uphi_rays);
  ur = array_1d<double>(ur_rays);
  uphi = array_1d<double>(uphi_rays);

  //----------------------------------------------------------------------------
  // prepare file
  
  h5_fid = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  //----------------------------------------------------------------------------
  // Store the header data
  
  value = F;
  dim_1d = 1;
  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "F", H5T_NATIVE_DOUBLE, dataspace,
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &value);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  value = p_Obs->t;
  dim_1d = 1;
  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "t_obs", H5T_NATIVE_DOUBLE, dataspace,
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &value);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  value = p_Obs->z;
  dim_1d = 1;
  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "z", H5T_NATIVE_DOUBLE, dataspace,
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &value);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  value = p_Obs->dL;
  dim_1d = 1;
  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "dL", H5T_NATIVE_DOUBLE, dataspace,
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &value);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  value = ur_min;
  dim_1d = 1;
  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "ur_min", H5T_NATIVE_DOUBLE, dataspace,
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &value);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  value = 0.; // this is always the same
  dim_1d = 1;
  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "uphi_min", H5T_NATIVE_DOUBLE, dataspace,
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &value);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  value = ur_max;
  dim_1d = 1;
  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "ur_max", H5T_NATIVE_DOUBLE, dataspace,
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &value);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  value = PI; // this is always the same
  dim_1d = 1;
  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "uphi_max", H5T_NATIVE_DOUBLE, dataspace,
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &value);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  value = R_50;
  dim_1d = 1;
  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "R_50", H5T_NATIVE_DOUBLE, dataspace,
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &value);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  value = R_75;
  dim_1d = 1;
  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "R_75", H5T_NATIVE_DOUBLE, dataspace,
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &value);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  value = R_95;
  dim_1d = 1;
  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "R_95", H5T_NATIVE_DOUBLE, dataspace,
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &value);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  value = R_99;
  dim_1d = 1;
  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "R_99", H5T_NATIVE_DOUBLE, dataspace,
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &value);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  value = R_100;
  dim_1d = 1;
  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "R_100", H5T_NATIVE_DOUBLE, dataspace,
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &value);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  int_value = ur_rays;
  dim_1d = 1;
  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "ur_rays", H5T_NATIVE_INT, dataspace,
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_INT, dataspace, dataspace,
    H5P_DEFAULT, &int_value);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  int_value = uphi_rays;
  dim_1d = 1;
  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "uphi_rays", H5T_NATIVE_INT, dataspace,
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_INT, dataspace, dataspace,
    H5P_DEFAULT, &int_value);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  //----------------------------------------------------------------------------
  // Store the main data
  
  dim_2d[0] = ur_rays;
  dim_2d[1] = uphi_rays;
  dataspace = H5Screate_simple(2, &dim_2d[0], NULL);
  dataset = H5Dcreate1(h5_fid, "I", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  // prepare emergent intensities for storing
  for(i_ur = 0; i_ur < ur_rays; i_ur++)
    for (i_uphi = 0; i_uphi < uphi_rays; i_uphi++)
    {
      #if LOCAL_SELF_ABSORPTION_ == ENABLED_
  
        I[i_ur][i_uphi] = ray[i_ur][i_uphi].I;
 
      #else

        tau = ray[i_ur][i_uphi].abdr;
        if (tau > 1e-3)
        {
          S = ray[i_ur][i_uphi].emdr / ray[i_ur][i_uphi].abdr;
          I[i_ur][i_uphi] = S * (1. - exp(-tau));
        }
        else
        {
          I[i_ur][i_uphi] = ray[i_ur][i_uphi].emdr * (1. - 0.5 * tau);
        }

      #endif
    }
  // store the intensities
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &I[0][0]);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  // coordinate information, r coordinates
  for (i_ur = 0; i_ur < ur_rays; i_ur++)
  {
    ur[i_ur] = ray[i_ur][0].ur;
    //printf("[ur[%d] = %e\n]", i_ur, ur[i_ur]);
  }
  dim_1d = ur_rays ;
  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "ur", H5T_NATIVE_DOUBLE, dataspace,
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, ur);
  H5Dclose(dataset);
  H5Sclose(dataspace);
  
  // coordinate information, phi coordinates
  for (i_uphi = 0; i_uphi < uphi_rays; i_uphi++)
    uphi[i_uphi] = ray[0][i_uphi].uphi;
  dim_1d = uphi_rays;
  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "uphi", H5T_NATIVE_DOUBLE, dataspace,
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, uphi);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  //----------------------------------------------------------------------------

  // deallocate memory and close the file
  H5Fclose(h5_fid);
  
  delete_array_2d(I);
  delete_array_1d(ur);
  delete_array_1d(uphi);
}

////////////////////////////////////////////////////////////////////////////////
// internal ( protected ) functions c_eds_2D
////////////////////////////////////////////////////////////////////////////////

double c_eds :: get_F_annulus(int iur)
{
  #if LOCAL_SELF_ABSORPTION_ == ENABLED_

    // on-axis observer
    if (uphi_rays == 1) 
      return 2. * PI * ray[iur][0].I * ray[iur][0].ur * ray[iur][0].ur;
        // 2PI is the domain size in Phi
        // ur^2 is the Jacobian for polar log r coordinates

    // off-axis observer
    double h = PI / (uphi_rays - 1);
    double F;
    int iuphi;

    // The loop implements a closed 6-point Newton-Cotes formula between 0 .. PI
    // formula reference: http://mathworld.wolfram.com/Newton-CotesFormulas.html
    F = -19. * ray[iur][0].I; // correct for overcounting
    
    for (iuphi = 0; iuphi < uphi_rays - 1; iuphi += 5)
    {
      F += 38. * ray[iur][iuphi].I +
        75. * ray[iur][iuphi + 1].I +
        50 * ray[iur][iuphi + 2].I
        50 * ray[iur][iuphi + 3].I
        75 * ray[iur][iuphi + 4].I;
    }
    F += 19. * ray[iur][iuphi]; // add outer closed boundary

    F = F * 5. / 288. * h

    return 2 * F * ray[iur][0].ur * ray[iur][0].ur;
      // factor 2 to account for the u_phi half from PI to 2 PI

  #else // LOCAL_SELF_ABSORPTION_ == DISABLED_

    // on-axis observer
    if (uphi_rays == 1) 
    {
      double S, tau;
  
      if (ray[iur][0].abdr > 1e-3)
      {
        S = ray[iur][0].emdr / ray[iur][0].abdr;
        tau = ray[iur][0].abdr;
        return S * (1. - exp(-tau)) * 2 * PI * ray[iur][0].ur * ray[iur][0].ur;
          // 2PI is the domain size in Phi
          // ur^2 is the Jacobian for polar log r coordinates
      }
      else
        return ray[iur][0].emdr * (1. - 0.5 * ray[iur][0].abdr) * 2 * PI * 
          ray[iur][0].ur * ray[iur][0].ur;
          // 2PI is the domain size in Phi
          // ur^2 is the Jacobian for polar log r coordinates
    }

    // off-axis observer
    double F, S, tau;
    double h = PI / (uphi_rays - 1);
    int iuphi;

    // The loop implements a closed 6-point Newton-Cotes formula between 0..PI
    // formula ref: http://mathworld.wolfram.com/Newton-CotesFormulas.html

    // correct for overshooting
    tau = ray[iur][0].abdr;
    if (tau > 1e-3)
    {
      S = ray[iur][0].emdr / ray [iur][0].abdr;
      F = - 19. * S * (1 - exp(-tau));
    }
    else
    {
      F = -19. * ray[iur][0].emdr * (1. - 0.5 * tau);
    }

    for (iuphi = 0; iuphi < uphi_rays - 1; iuphi += 5)
    {
      // point 0
      tau = ray[iur][iuphi].abdr;
      if (tau > 1e-3)
      {
        S = ray[iur][iuphi].emdr / ray [iur][iuphi].abdr;
        F += 38. * S * (1 - exp(-tau));
      }
      else
      {
        F += 38. * ray[iur][iuphi].emdr * (1. - 0.5 * tau);
      }

      // point 1
      tau = ray[iur][iuphi + 1].abdr;
      if (tau > 1e-3)
      {
        S = ray[iur][iuphi + 1].emdr / ray [iur][iuphi + 1].abdr;
        F += 75. * S * (1 - exp(-tau));
      }
      else
      {
        F += 75. * ray[iur][iuphi + 1].emdr * (1. - 0.5 * tau);
      }
      
      // point 2
      tau = ray[iur][iuphi + 2].abdr;
      if (tau > 1e-3)
      {
        S = ray[iur][iuphi + 2].emdr / ray [iur][iuphi + 2].abdr;
        F += 50. * S * (1 - exp(-tau));
      }
      else
      {
        F += 50. * ray[iur][iuphi + 2].emdr * (1. - 0.5 * tau);
      }

      // point 3
      tau = ray[iur][iuphi + 3].abdr;
      if (tau > 1e-3)
      {
        S = ray[iur][iuphi + 3].emdr / ray [iur][iuphi + 3].abdr;
        F += 50. * S * (1 - exp(-tau));
      }
      else
      {
        F += 50. * ray[iur][iuphi + 3].emdr * (1. - 0.5 * tau);
      }

      // point 4
      tau = ray[iur][iuphi + 4].abdr;
      if (tau > 1e-3)
      {
        S = ray[iur][iuphi + 4].emdr / ray [iur][iuphi + 4].abdr;
        F += 75. * S * (1 - exp(-tau));
      }
      else
      {
        F += 75. * ray[iur][iuphi + 4].emdr * (1. - 0.5 * tau);
      }
    }
    
    // add the closed end boundary
    tau = ray[iur][iuphi].abdr;
    if (tau > 1e-3)
    {
      S = ray[iur][iuphi].emdr / ray [iur][iuphi].abdr;
      F += 19. * S * (1 - exp(-tau));
    }
    else
    {
      F += 19. * ray[iur][iuphi].emdr * (1. - 0.5 * tau);
    }

    F = F * 5 / 288. * h;      
      
    // include a factor 2 to add the other half (phi = PI .. 2PI), which mirrors
    // the first half (phi = 0.. PI)
    return 2. * F * ray[iur][0].ur * ray[iur][0].ur;

  #endif
}

////////////////////////////////////////////////////////////////////////////////

double c_eds :: get_F_annulus_lores_t(int iur)
{
  #if LOCAL_SELF_ABSORPTION_ == ENABLED_

    // on-axis observer
    if (uphi_rays == 1) 
      return ray[iur][0].I_lores * 2 * PI * ray[iur][0].ur * ray[iur][0].ur;
        // 2PI is the domain size in Phi
        // ur^2 is the Jacobian for polar log r coordinates

    // off-axis observer
    double h = PI / (uphi_rays - 1);
    double F = 0.;
    int iuphi;

    // The loop implements a closed 6-point Newton-Cotes formula between 0 .. PI
    // formula reference: http://mathworld.wolfram.com/Newton-CotesFormulas.html
    F = -19. * ray[iur][0].I_lores; // correct for overcounting
    
    for (iuphi = 0; iuphi < uphi_rays - 1; iuphi += 5)
    {
      F += 38. * ray[iur][iuphi].I_lores +
        75. * ray[iur][iuphi + 1].I_lores +
        50 * ray[iur][iuphi + 2].I_lores
        50 * ray[iur][iuphi + 3].I_lores
        75 * ray[iur][iuphi + 4].I_lores;
    }
    F += 19. * ray[iur][iuphi]; // add outer closed boundary

    F = F * 5. / 288. * h

    return 2 * F * ray[iur][0].ur * ray[iur][0].ur;
      // factor 2 to account for the u_phi half from PI to 2 PI

  #else // LOCAL_SELF_ABSORPTION_ == DISABLED_



    // on-axis observer
    if (uphi_rays == 1) 
    {
      double S, tau;
  
      if (ray[iur][0].abdr_lores > 1e-3)
      {
        S = ray[iur][0].emdr_lores / ray[iur][0].abdr_lores;
        tau = ray[iur][0].abdr_lores;
        return S * (1. - exp(-tau)) * 2 * PI * ray[iur][0].ur * ray[iur][0].ur;
          // 2PI is the domain size in Phi
          // ur^2 is the Jacobian for polar log r coordinates
      }
      else
        return ray[iur][0].emdr_lores * (1. - 0.5 * ray[iur][0].abdr_lores) * 
          2 * PI * ray[iur][0].ur * ray[iur][0].ur;
          // 2PI is the domain size in Phi
          // ur^2 is the Jacobian for polar log r coordinates
    }

    // off-axis observer
    double F, S, tau;
    double h = PI / (uphi_rays - 1);
    int iuphi;

    // The loop implements a closed 6-point Newton-Cotes formula between 0..PI
    // formula ref: http://mathworld.wolfram.com/Newton-CotesFormulas.html

    // correct for overshooting
    tau = ray[iur][0].abdr_lores;
    if (tau > 1e-3)
    {
      S = ray[iur][0].emdr_lores / ray [iur][0].abdr_lores;
      F = - 19. * S * (1 - exp(-tau));
    }
    else
    {
      F = -19. * ray[iur][0].emdr_lores * (1. - 0.5 * tau);
    }

    for (iuphi = 0; iuphi < uphi_rays - 1; iuphi += 5)
    {
      // point 0
      tau = ray[iur][iuphi].abdr_lores;
      if (tau > 1e-3)
      {
        S = ray[iur][iuphi].emdr_lores / ray [iur][iuphi].abdr_lores;
        F += 38. * S * (1 - exp(-tau));
      }
      else
      {
        F += 38. * ray[iur][iuphi].emdr_lores * (1. - 0.5 * tau);
      }

      // point 1
      tau = ray[iur][iuphi + 1].abdr_lores;
      if (tau > 1e-3)
      {
        S = ray[iur][iuphi + 1].emdr_lores / ray [iur][iuphi + 1].abdr_lores;
        F += 75. * S * (1 - exp(-tau));
      }
      else
      {
        F += 75. * ray[iur][iuphi + 1].emdr_lores * (1. - 0.5 * tau);
      }
      
      // point 2
      tau = ray[iur][iuphi + 2].abdr;
      if (tau > 1e-3)
      {
        S = ray[iur][iuphi + 2].emdr_lores / ray [iur][iuphi + 2].abdr_lores;
        F += 50. * S * (1 - exp(-tau));
      }
      else
      {
        F += 50. * ray[iur][iuphi + 2].emdr_lores * (1. - 0.5 * tau);
      }

      // point 3
      tau = ray[iur][iuphi + 3].abdr_lores;
      if (tau > 1e-3)
      {
        S = ray[iur][iuphi + 3].emdr_lores / ray [iur][iuphi + 3].abdr_lores;
        F += 50. * S * (1 - exp(-tau));
      }
      else
      {
        F += 50. * ray[iur][iuphi + 3].emdr_lores * (1. - 0.5 * tau);
      }

      // point 4
      tau = ray[iur][iuphi + 4].abdr_lores;
      if (tau > 1e-3)
      {
        S = ray[iur][iuphi + 4].emdr_lores / ray [iur][iuphi + 4].abdr_lores;
        F += 75. * S * (1 - exp(-tau));
      }
      else
      {
        F += 75. * ray[iur][iuphi + 4].emdr_lores * (1. - 0.5 * tau);
      }
    }
    
    // add the closed end boundary
    tau = ray[iur][iuphi].abdr_lores;
    if (tau > 1e-3)
    {
      S = ray[iur][iuphi].emdr_lores / ray [iur][iuphi].abdr_lores;
      F += 19. * S * (1 - exp(-tau));
    }
    else
    {
      F += 19. * ray[iur][iuphi].emdr_lores * (1. - 0.5 * tau);
    }

    F = F * 5 / 288. * h;      
      
    // include a factor 2 to add the other half (phi = PI .. 2PI), which mirrors
    // the first half (phi = 0.. PI)
    return 2. * F * ray[iur][0].ur * ray[iur][0].ur;

  #endif
}

////////////////////////////////////////////////////////////////////////////////

double c_eds :: get_F_annulus_lores_phi(int iur)
{
  #if LOCAL_SELF_ABSORPTION_ == ENABLED_

    // on-axis observer
    if (uphi_rays == 1) 
      return 2. * PI * ray[iur][0].I * ray[iur][0].ur * ray[iur][0].ur;
        // 2PI is the domain size in Phi
        // ur^2 is the Jacobian for polar log r coordinates

    // off-axis observer
    double h = 2 * PI / (uphi_rays - 1); // note increased stepsize
    double F;
    int iuphi;

    // The loop implements a closed 6-point Newton-Cotes formula between 0 .. PI
    // formula reference: http://mathworld.wolfram.com/Newton-CotesFormulas.html
    F = -19. * ray[iur][0].I; // correct for overcounting
    
    for (iuphi = 0; iuphi < uphi_rays - 1; iuphi += 10)
    {
      F += 38. * ray[iur][iuphi].I +
        75. * ray[iur][iuphi + 2].I +
        50 * ray[iur][iuphi + 4].I
        50 * ray[iur][iuphi + 6].I
        75 * ray[iur][iuphi + 8].I;
    }
    F += 19. * ray[iur][iuphi]; // add outer closed boundary

    F = F * 5. / 288. * h

    return 2 * F * ray[iur][0].ur * ray[iur][0].ur;
      // factor 2 to account for the u_phi half from PI to 2 PI

  #else // LOCAL_SELF_ABSORPTION_ == DISABLED_

    // on-axis observer
    if (uphi_rays == 1) 
    {
      double S, tau;
  
      if (ray[iur][0].abdr > 1e-3)
      {
        S = ray[iur][0].emdr / ray[iur][0].abdr;
        tau = ray[iur][0].abdr;
        return S * (1. - exp(-tau)) * 2 * PI * ray[iur][0].ur * ray[iur][0].ur;
          // 2PI is the domain size in Phi
          // ur^2 is the Jacobian for polar log r coordinates
      }
      else
        return ray[iur][0].emdr * (1. - 0.5 * ray[iur][0].abdr) * 2 * PI * 
          ray[iur][0].ur * ray[iur][0].ur;
          // 2PI is the domain size in Phi
          // ur^2 is the Jacobian for polar log r coordinates
    }

    // off-axis observer
    double F, S, tau;
    double h = 2 * PI / (uphi_rays - 1);
    int iuphi;

    // The loop implements a closed 6-point Newton-Cotes formula between 0..PI
    // formula ref: http://mathworld.wolfram.com/Newton-CotesFormulas.html

    // correct for overshooting
    tau = ray[iur][0].abdr;
    if (tau > 1e-3)
    {
      S = ray[iur][0].emdr / ray [iur][0].abdr;
      F = - 19. * S * (1 - exp(-tau));
    }
    else
    {
      F = -19. * ray[iur][0].emdr * (1. - 0.5 * tau);
    }

    for (iuphi = 0; iuphi < uphi_rays - 1; iuphi += 10)
    {
      // point 0
      tau = ray[iur][iuphi].abdr;
      if (tau > 1e-3)
      {
        S = ray[iur][iuphi].emdr / ray [iur][iuphi].abdr;
        F += 38. * S * (1 - exp(-tau));
      }
      else
      {
        F += 38. * ray[iur][iuphi].emdr * (1. - 0.5 * tau);
      }

      // point 1
      tau = ray[iur][iuphi + 2].abdr;
      if (tau > 1e-3)
      {
        S = ray[iur][iuphi + 2].emdr / ray [iur][iuphi + 2].abdr;
        F += 75. * S * (1 - exp(-tau));
      }
      else
      {
        F += 75. * ray[iur][iuphi + 2].emdr * (1. - 0.5 * tau);
      }
      
      // point 2
      tau = ray[iur][iuphi + 4].abdr;
      if (tau > 1e-3)
      {
        S = ray[iur][iuphi + 4].emdr / ray [iur][iuphi + 4].abdr;
        F += 50. * S * (1 - exp(-tau));
      }
      else
      {
        F += 50. * ray[iur][iuphi + 4].emdr * (1. - 0.5 * tau);
      }

      // point 3
      tau = ray[iur][iuphi + 6].abdr;
      if (tau > 1e-3)
      {
        S = ray[iur][iuphi + 6].emdr / ray [iur][iuphi + 6].abdr;
        F += 50. * S * (1 - exp(-tau));
      }
      else
      {
        F += 50. * ray[iur][iuphi + 6].emdr * (1. - 0.5 * tau);
      }

      // point 4
      tau = ray[iur][iuphi + 8].abdr;
      if (tau > 1e-3)
      {
        S = ray[iur][iuphi + 8].emdr / ray [iur][iuphi + 8].abdr;
        F += 75. * S * (1 - exp(-tau));
      }
      else
      {
        F += 75. * ray[iur][iuphi + 8].emdr * (1. - 0.5 * tau);
      }
    }
    
    // add the closed end boundary
    tau = ray[iur][iuphi].abdr;
    if (tau > 1e-3)
    {
      S = ray[iur][iuphi].emdr / ray [iur][iuphi].abdr;
      F += 19. * S * (1 - exp(-tau));
    }
    else
    {
      F += 19. * ray[iur][iuphi].emdr * (1. - 0.5 * tau);
    }

    F = F * 5 / 288. * h;      
      
    // include a factor 2 to add the other half (phi = PI .. 2PI), which mirrors
    // the first half (phi = 0.. PI)
    return 2. * F * ray[iur][0].ur * ray[iur][0].ur;

  #endif
}

////////////////////////////////////////////////////////////////////////////////

double c_eds :: get_total_flux()
{
  int iur; // , iuphi;
  double h = (log(ur_max) - log(ur_min)) / (double) (ur_rays - 1);

  // Compute the flux by summing over the annuli. Since the radial coordinates
  // are logarithmically spaced, we can only use a higher-order method for
  // equally spaced abscissas if we integrate over ln r instead of r. The
  // integral then becomes Int r^2 I d u_phi d ln _ur in polar coordinates on
  // the EDS.

  // The loop implements a closed 6-point Newton-Cotes formula. However, since
  // the EDS is assumed to cover more than the full area on the sky of the
  // source, the outer domain limit intensity value is assumed zero.
  // formula reference: http://mathworld.wolfram.com/Newton-CotesFormulas.html
  
  F = -19. * get_F_annulus(0); // start with correction for overcounting
  for (iur = 0; iur < ur_rays; iur += 5)
  {
    F += (38. * get_F_annulus(iur) +
      75 * get_F_annulus(iur + 1) +
      50 * get_F_annulus(iur + 2) +
      50 * get_F_annulus(iur + 3) +
      75 * get_F_annulus(iur + 4));
  }
  F = F * 5 / 288. * h; // apply scale factor and domain step size h
  
  // correct the flux for distance
  F = F * (1.0 + p_Obs->z) / (p_Obs->dL * p_Obs->dL);
  
  return F;
}

////////////////////////////////////////////////////////////////////////////////

double c_eds :: get_F_r_error()
{
  // returns measure of error due to radial resolution
  double F = get_total_flux();
  double Flores;
  int iur;
  double h = (log(ur_max) - log(ur_min)) * 2. / (double) (ur_rays - 1);

  Flores = -19. * get_F_annulus(0); // start with overcounting correction
  for (iur = 0; iur < ur_rays; iur += 10)
  {
    F += (38. * get_F_annulus(iur) +
      75 * get_F_annulus(iur + 2) +
      50 * get_F_annulus(iur + 3) +
      50 * get_F_annulus(iur + 6) +
      75 * get_F_annulus(iur + 8));
  }
  F = F * 5 / 288. * h; // apply scale factor and domain step size h
  
  // correct the flux for distance
  F = F * (1.0 + p_Obs->z) / (p_Obs->dL * p_Obs->dL);

  // Flores = (2.0 * Flores) * (1.0 + p_Obs->z) / (p_Obs->dL * p_Obs->dL);
    
  return 2. * fabs((Flores - F) / (Flores + F + 1e-200));
}

////////////////////////////////////////////////////////////////////////////////

double c_eds :: get_F_phi_error()
{
  // returns measure of error due to angular resolution
  double F = get_total_flux();

  int iur;
  double F_lores = 0.;
  double h = (log(ur_max) - log(ur_min)) / (double) (ur_rays - 1);

  // The loop implements a closed 6-point Newton-Cotes formula. However, since
  // the EDS is assumed to cover more than the full area on the sky of the
  // source, the outer domain limit intensity value is assumed zero.
  // formula reference: http://mathworld.wolfram.com/Newton-CotesFormulas.html
  
  F_lores = -19. * get_F_annulus_lores_phi(0); 
    // start with correction for overcounting
  for (iur = 0; iur < ur_rays; iur += 5)
  {
    F_lores += (38. * get_F_annulus_lores_phi(iur) +
      75 * get_F_annulus_lores_phi(iur + 1) +
      50 * get_F_annulus_lores_phi(iur + 2) +
      50 * get_F_annulus_lores_phi(iur + 3) +
      75 * get_F_annulus_lores_phi(iur + 4));
  }
  F_lores = F_lores * 5 / 288. * h; // apply scale factor and domain step size h
  
  // correct the flux for distance
  F_lores = F_lores * (1.0 + p_Obs->z) / (p_Obs->dL * p_Obs->dL);
  
  return 2. * fabs((F_lores - F) / (F_lores + F + 1e-200));
}

////////////////////////////////////////////////////////////////////////////////

double c_eds :: get_F_t_error()
{
  // returns measure of error due to temporal resolution
  double F = get_total_flux();

  int iur;
  double F_lores = 0.;
  double h = (log(ur_max) - log(ur_min)) / (double) (ur_rays - 1);

  // The loop implements a closed 6-point Newton-Cotes formula. However, since
  // the EDS is assumed to cover more than the full area on the sky of the
  // source, the outer domain limit intensity value is assumed zero.
  // formula reference: http://mathworld.wolfram.com/Newton-CotesFormulas.html
  
  F_lores = -19. * get_F_annulus_lores_t(0); 
    // start with correction for overcounting
  for (iur = 0; iur < ur_rays; iur += 5)
  {
    F_lores += (38. * get_F_annulus_lores_t(iur) +
      75 * get_F_annulus_lores_t(iur + 1) +
      50 * get_F_annulus_lores_t(iur + 2) +
      50 * get_F_annulus_lores_t(iur + 3) +
      75 * get_F_annulus_lores_t(iur + 4));
  }
  F_lores = F_lores * 5 / 288. * h; // apply scale factor and domain step size h
  
  // correct the flux for distance
  F_lores = F_lores * (1.0 + p_Obs->z) / (p_Obs->dL * p_Obs->dL);
  
  return 2. * fabs((F_lores - F) / (F_lores + F + 1e-200));  
}

////////////////////////////////////////////////////////////////////////////////

void c_eds :: set_R()
{
  int iur;
  double df;
  
  R_50 = -1.;
  R_75 = -1.;
  R_95 = -1.;
  R_99 = -1.;
  R_100 = -1;
  
  double F_total = 0.0, f = 0.0; // 'f' is fractional flux

  // get total
  for (iur = 0; iur < ur_rays; iur++)
    F_total += get_F_annulus(iur);
  
  // get fractional results
  R_99 = -1.0; R_95 = -1.0; R_75 = -1.0; R_50 = -1.0;
  
  for (iur = 0; iur < ur_rays; iur++)
  {
    df = get_F_annulus(iur);
    f += df;
    if (f > 0.5 * F_total and R_50 < 0.0) R_50 = ray[iur][0].ur;
    if (f > 0.75 * F_total and R_75 < 0.0) R_75 = ray[iur][0].ur;
    if (f > 0.95 * F_total and R_95 < 0.0) R_95 = ray[iur][0].ur;
    if (f > 0.99 * F_total and R_99 < 0.0) R_99 = ray[iur][0].ur;
    if (df > 0.) R_100 = ray[iur][0].ur;
  }
}

////////////////////////////////////////////////////////////////////////////////

void c_eds :: reset()
{
  int iur, iuphi;
  double dlnur, duphi;
  
  //----------------------------------------------------------------------------
  // set coordinates assumes that lnur_max and lnur_min are set
  if (uphi_rays == 1)
    duphi = PI;
  else 
    duphi = (PI - 0.0) / (double) (uphi_rays - 1);
  dlnur = (log(ur_max) - log(ur_min)) / (double) (ur_rays - 1);

  // set unit vectors on the EDS plane in the frame of the grid, 
  // assuming theta_obs related variables are set
  eux[0] = p_Obs->costheta; eux[1] = 0.0; eux[2] = -p_Obs->sintheta;
  euy[0] = 0.0; euy[1] = 1.0; euy[2] = 0.0;

  for (iur = 0; iur < ur_rays; iur++)
  {
    for (iuphi = 0; iuphi < uphi_rays; iuphi++)
    {
      ray[iur][iuphi].ur = exp(log(ur_min) + iur * dlnur);
      ray[iur][iuphi].uphi = iuphi * duphi;

      ray[iur][iuphi].ux = ray[iur][iuphi].ur * cos(ray[iur][iuphi].uphi); 
      ray[iur][iuphi].uy = ray[iur][iuphi].ur * sin(ray[iur][iuphi].uphi);
      
      // reset radiation
      ray[iur][iuphi].em = 0.0;
      ray[iur][iuphi].ab = 0.0;
      
      #if LOCAL_SELF_ABSORPTION_ == DISABLED_

        ray[iur][iuphi].emdr = 0.0;
        ray[iur][iuphi].abdr = 0.0;

        ray[iur][iuphi].emdr_lores = 0.0;
        ray[iur][iuphi].abdr_lores = 0.0;
      
      #else
      
        ray[iur][iuphi].I = 0.0;
        ray[iur][iuphi].I_lores = 0.0;
        
      #endif
      
    }
  }
  
  F = 0.;

  return;
}
