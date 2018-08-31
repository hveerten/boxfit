////////////////////////////////////////////////////////////////////////////////
//
// eds_IO_2D.cpp
//
// Created November 30, 2010 by HJvE
// last modified January 13, 2012 by HJvE
//
////////////////////////////////////////////////////////////////////////////////

#include "eds_IO.h"

////////////////////////////////////////////////////////////////////////////////

c_eds_IO :: c_eds_IO()
{
  save_EDS = false;
}

////////////////////////////////////////////////////////////////////////////////

c_eds_IO :: ~c_eds_IO()
{
}

////////////////////////////////////////////////////////////////////////////////

void c_eds_IO :: save(const char* filenamebase, int entry)
{
  char filename[strlen(filenamebase)+4+strlen(".h5")+1];
  hid_t h5_fid, dataset, dataspace;
  hsize_t dim_1d, dim_2d[2];
  //int ix, iy;
  //double dx, dy, x, y, r, phi;
  //double **Ireg;

  //----------------------------------------------------------------------------

  // create file
  sprintf(filename,"%s%04d.h5", filenamebase, entry);
  h5_fid = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  //----------------------------------------------------------------------------
  // save general info
  
  dim_1d = 1;
  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "nu_obs", H5T_NATIVE_DOUBLE, dataspace, 
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &p_Obs->nu);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "t_obs", H5T_NATIVE_DOUBLE, dataspace, 
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &p_Obs->t);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "theta_obs", H5T_NATIVE_DOUBLE, dataspace, 
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &p_Obs->theta);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "F", H5T_NATIVE_DOUBLE, dataspace, 
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &F);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "D_L", H5T_NATIVE_DOUBLE, dataspace, 
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &p_Obs->dL);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "z", H5T_NATIVE_DOUBLE, dataspace, 
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &p_Obs->z);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "ur_rays", H5T_NATIVE_INT, dataspace, 
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_INT, dataspace, dataspace,
    H5P_DEFAULT, &ur_rays);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "uphi_rays", H5T_NATIVE_INT, dataspace, 
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_INT, dataspace, dataspace,
    H5P_DEFAULT, &uphi_rays);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  /*
  // save eds extent
  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "ulnrImax", H5T_NATIVE_DOUBLE, dataspace, 
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &ulnrImax);
  H5Dclose(dataset);
  H5Sclose(dataspace);
  */
  
  //----------------------------------------------------------------------------
  // full dataset on polar grid

  dim_2d[0] = ur_rays; dim_2d[1] = uphi_rays;
  dataspace = H5Screate_simple(2, &dim_2d[0], NULL);
  dataset = H5Dcreate1(h5_fid, "I", H5T_NATIVE_DOUBLE, dataspace, 
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &ray[0][0]);
  H5Dclose(dataset);
  H5Sclose(dataspace);
  
  //----------------------------------------------------------------------------
  // dataset on low-resolution Cartesian preview grid


  /*
  COMMENTED OUT COMPLETELY BECAUSE IT IS WRONG! NO CORRECT WEIGHTS!!!!!


  // create Cartesian grid
  Ireg = array_2d<double>(ur_rays * 2 - 1, ur_rays * 2 - 1);
  
  x = -exp(ulnrImax); y = -exp(ulnrImax);
  dx = -2 * x / (ur_rays * 2);
  x += 0.5 * dx;
  dy = -2 * y / (ur_rays * 2);
  y += 0.5 * dy;
  for (ix = 0; ix < ur_rays * 2 - 1; ix++)
  {
    for (iy = 0; iy < ur_rays * 2 - 1; iy++)
    {
       r = sqrt(x * x + y * y);
       phi = atan(y / x);
       Ireg[ix][iy] = get_I(r, phi);
       y += dy;
    }
    x += dx;
  }

  // save regular grid
  dim_2d[0] = 2 * ur_rays - 1; dim_2d[1] = 2 * ur_rays - 1;
  dataspace = H5Screate_simple(2, &dim_2d[0], NULL);
  dataset = H5Dcreate1(h5_fid, "Ireg", H5T_NATIVE_DOUBLE, dataspace, 
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &Ireg[0][0]);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  // release regular grid memory
  delete_array_2d(Ireg);  */
  
  // close file
  H5Fclose(h5_fid);
}
