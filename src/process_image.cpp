#include <stdlib.h>
#include <stdio.h>
#include "hdf5.h"
#include "arraytools.h"
#include <cmath>
#include "extramath.h"

int main(int argc, char* argv[])
{
  hid_t dataset, dataspace;
  hid_t h5_fid;
  char filename[255] = "image0050.h5"; // some default number
  
  hsize_t dim_2d[2], dim_1d;
  
  // content already in the image
  double **I; // 2D array containing the intensities in radial coordinates
  double *ur; // 1D array containing the radial coordinates
  double *uphi; // 1D array containing phi coordinates
  int ur_rays; // number of rays in r direction
  int uphi_rays; // number of rays in phi direction
  double ur_max; // upper r boundary of the eds. 'u' denotes coord in eds frame
  double ur_min; // lower boundary of eds.
  double R_100;
  double t_obs;
  
  double value; // for writing scalars to new hdf5 file
  int intvalue; // same as above, but for integers
  
  // additional content
  double *Ix, *Iy; // projected intensities on axes
  double **Ixy; // full map of intensities in Cartesian coordinates
  int xy_res;
  double I_total; // total summed intentsity (not including area information,
    // so not equal to flux
  double Ix_left, Ix_right, Ix_total;
  double Iy_left, Iy_total, Iy_right;
  
  int i_ux, i_uy, i_ur, i_uphi;
  double *ux, *uy; // x,y-coordinates of Cartesian version
  double dux;
  double ur_local, uphi_local, lnur_local;

  // open the image file for reading AND writing
  if (argc > 1)
    sprintf(filename, "%s", argv[1]);
  else
  {
    printf("please provide a filename for processing, e.g. process_image "
      "image0000.h5\n");
    printf("Now picking image0093.h5 as default, but file might not exist.\n");
    sprintf(filename, "image0093.h5");
  }
  
  printf("Analyzing image file: %s\n", filename);
  h5_fid = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (h5_fid < 0)
  {
    printf("(process_image) Error. Could not open file %s\n",filename);
    fflush(stdout);
    exit(1);
  }

  //----------------------------------------------------------------------------
  // read the relevant content
  
  dataset = H5Dopen1(h5_fid, "t_obs");
  dataspace = H5Dget_space(dataset);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
    &t_obs);
  H5Dclose(dataset);
  H5Sclose(dataspace);
  
  
  dataset = H5Dopen1(h5_fid, "ur_rays");
  dataspace = H5Dget_space(dataset);
  H5Dread(dataset, H5T_NATIVE_INT, dataspace, dataspace, H5P_DEFAULT, 
    &ur_rays);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  dataset = H5Dopen1(h5_fid, "ur_max");
  dataspace = H5Dget_space(dataset);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
    &ur_max);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  dataset = H5Dopen1(h5_fid, "R_100");
  dataspace = H5Dget_space(dataset);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
    &R_100);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  dataset = H5Dopen1(h5_fid, "uphi_rays");
  dataspace = H5Dget_space(dataset);
  H5Dread(dataset, H5T_NATIVE_INT, dataspace, dataspace, H5P_DEFAULT, 
    &uphi_rays);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  dataset = H5Dopen1(h5_fid, "ur_min");
  dataspace = H5Dget_space(dataset);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT,
    &ur_min);
  H5Dclose(dataset);
  H5Sclose(dataspace);
  
  I = array_2d<double>(ur_rays, uphi_rays);
  dataset = H5Dopen1(h5_fid, "I");
  dataspace = H5Dget_space(dataset);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT,
    &I[0][0]);
  H5Dclose(dataset);
  H5Sclose(dataspace);
  
  ur = array_1d<double>(ur_rays);
  dataset = H5Dopen1(h5_fid, "ur");
  dataspace = H5Dget_space(dataset);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, ur);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  uphi = array_1d<double>(uphi_rays);
  dataset = H5Dopen1(h5_fid, "uphi");
  dataspace = H5Dget_space(dataset);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, uphi);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  //----------------------------------------------------------------------------
  // print file summary
  
  printf("t_obs = %e s = %e days\n", t_obs, t_obs / (24. * 60 * 60));
  printf("ur_max = %e\n", ur_max);
  printf("ur_rays = %d\n", ur_rays);
  printf("uphi_rays = %d\n", uphi_rays);
  printf("R_100 = %e\n", R_100);
  
  // temporary hack to zoom in on the actual image on the sky
  if (R_100 > 0.) ur_max = R_100;
  
  //----------------------------------------------------------------------------
  // generate additional content
  xy_res = 2 * ur_rays - 1; // pos and negative, not double-counting the 
    // central point

  //xy_res *= 5;

  Ix = array_1d<double>(xy_res);
  Iy = array_1d<double>(xy_res);
  ux = array_1d<double>(xy_res);
  uy = array_1d<double>(xy_res);
  Ixy = array_2d<double>(xy_res, xy_res);
  
  // clean up projected intensity arrays first
  for (i_ux = 0; i_ux < xy_res; i_ux++)
  {
    Ix[i_ux] = 0.;
    Iy[i_ux] = 0.;
  }
  
  // set up x,y coordinates
  dux = ur_max / (ur_rays - 1);
  for (i_ux = 0; i_ux < xy_res; i_ux++)
  {
    ux[i_ux] = -ur_max + i_ux * dux;
    uy[i_ux] = -ur_max + i_ux * dux;
    //printf("ux[%d] = %e\n", i_ux, ux[i_ux]);
  }
  
  for (i_ux = 0; i_ux < xy_res; i_ux ++)
  {
    //i_uy = xy_res / 2;
    for (i_uy = 0; i_uy < xy_res; i_uy ++)
    {
      ur_local = sqrt(ux[i_ux] * ux[i_ux] + uy[i_uy] * uy[i_uy]);
      
      // deal with first quadrant, including y = 0
      if (ux[i_ux] > 0 and uy[i_uy] >= 0)
      {
        if (uy[i_uy] < 1e-8) uphi_local = 0.0;
        else uphi_local = atan(uy[i_uy] / ux[i_ux]);
      }
      else // 2nd quadrant
      if (ux[i_ux] <= 0 and uy[i_uy] > 0)
      {
        if (fabs(ux[i_ux]) < 1e-8) uphi_local = 0.5 * PI;
        else uphi_local = PI - atan(-uy[i_uy] / ux[i_ux]);
      }
      else
      if (ux[i_ux] < 0 and uy[i_uy] <= 0)
      {
        if (fabs(uy[i_uy]) < 1e-8) uphi_local = PI;
        else uphi_local = PI + atan(uy[i_uy] / ux[i_ux]);
      }
      else
      {
        if (ux[i_ux] < 1e-8) uphi_local = 1.5 * PI;
        else uphi_local = 2.0 * PI - atan(-uy[i_uy] / ux[i_ux]);
      }      

      if (uphi_local < 0)
        uphi_local = -uphi_local;
      if (uphi_local > PI)
        uphi_local = PI - (uphi_local - PI);
      
      // find the appropriate r, phi entries
      i_uphi = (uphi_local / PI * (uphi_rays)); // type-casting to int drops
        // everything beyond the decimal point
      if (i_uphi == uphi_rays)
          i_uphi--; // deals with the case of being at PI exactly
      
      if (ur_local <= ur_min) 
        i_ur = 0;
      else
      {
        lnur_local = log(ur_local);
        i_ur = (lnur_local - log(ur_min)) / (log(ur_max) - log(ur_min)) * 
          ur_rays;
      }
      
      if (i_ur < ur_rays)
        Ixy[i_ux][i_uy] = I[i_ur][i_uphi];
      else
        Ixy[i_ux][i_uy] = 0;// Cartesian square includes parts not in EDS circle
      //if (Ixy[i_ux][i_uy] > 0.)
      //  printf("(x = %e, y= %e) -> (r = %e, phi = %e). i_r = %d, i_phi = %d, I = %e\n ", ux[i_ux], uy[i_uy], ur_local, uphi_local, i_ur, i_uphi, Ixy[i_ux][i_uy]);
    }
  }
  
  // project the intensity profiles on the x and y axes
  I_total = 0;
  for (i_ux = 0; i_ux < xy_res; i_ux++)
  {
    for (i_uy = 0; i_uy < xy_res; i_uy++)
    {
       Ix[i_ux] += Ixy[i_ux][i_uy];
    }
    I_total += Ix[i_ux];
  }

  for (i_uy = 0; i_uy < xy_res; i_uy++)
  {
    for (i_ux = 0; i_ux < xy_res; i_ux++)
    {
       Iy[i_uy] += Ixy[i_ux][i_uy];
    }
  }
  
  // compute the size both in x and y directions
  Ix_total = 0;
  Ix_right = xy_res - 1;
  for (i_ux = 0; i_ux < xy_res; i_ux++)
  {
    Ix_total += Ix[i_ux];
    if (Ix_total > 0.99 * I_total)
    {
      Ix_right = ux[i_ux];
      break;
    }    
  }

  Ix_total = 0;
  Ix_left = xy_res - 1;
  for (i_ux = xy_res - 1; i_ux >= 0; i_ux--)
  {
    Ix_total += Ix[i_ux];
    if (Ix_total > 0.99 * I_total)
    {
      Ix_left = ux[i_ux];
      break;
    }    
  }

  Iy_total = 0;
  Iy_right = xy_res - 1;
  for (i_uy = 0; i_uy < xy_res; i_uy++)
  {
    Iy_total += Iy[i_uy];
    if (Iy_total > 0.99 * I_total)
    {
      Iy_right = uy[i_uy];
      break;
    }    
  }

  Iy_total = 0;
  Iy_left = xy_res - 1;
  for (i_uy = xy_res - 1; i_uy >= 0; i_uy--)
  {
    Iy_total += Iy[i_uy];
    if (Iy_total > 0.99 * I_total)
    {
      Iy_left = uy[i_uy];
      break;
    }    
  }
  
  printf("x_left = %e, x_right = %e, size = %e\n", Ix_left, Ix_right, Ix_right - Ix_left);
  printf("y_left = %e, y_right = %e, size = %e\n", Iy_left, Iy_right, Iy_right - Iy_left);

  // close the input file
  H5Fclose(h5_fid);
  
  //----------------------------------------------------------------------------
  // write additional content to file

  // create a new file for output
  sprintf(filename, "out.h5");
  h5_fid = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (h5_fid < 0)
  {
    printf("(process_image) Error. Could not open file %s "
      "(did it exist already?)\n",filename);
    fflush(stdout);
    exit(1);
  }

  // the old stuff, now into the new file
  intvalue = ur_rays;
  dim_1d = 1;
  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "ur_rays", H5T_NATIVE_INT, dataspace,
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_INT, dataspace, dataspace,
    H5P_DEFAULT, &intvalue);
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

  value = t_obs;
  dim_1d = 1;
  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "t_obs", H5T_NATIVE_DOUBLE, dataspace,
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

  intvalue = uphi_rays;
  dim_1d = 1;
  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "uphi_rays", H5T_NATIVE_INT, dataspace,
    H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_INT, dataspace, dataspace,
    H5P_DEFAULT, &intvalue);
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

  dim_2d[0] = ur_rays;
  dim_2d[1] = uphi_rays;
  dataspace = H5Screate_simple(2, &dim_2d[0], NULL);
  dataset = H5Dcreate1(h5_fid, "I", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &I[0][0]);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  dim_1d = ur_rays;
  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "ur", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &ur[0]);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  dim_1d = uphi_rays;
  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "uphi", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &ur[0]);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  //----------------------------------------------------------------------------
  // new stuff  
  
  dim_2d[0] = xy_res;
  dim_2d[1] = xy_res;
  dataspace = H5Screate_simple(2, &dim_2d[0], NULL);
  dataset = H5Dcreate1(h5_fid, "Ixy", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &Ixy[0][0]);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  dim_1d = xy_res;
  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "Ix", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &Ix[0]);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  dim_1d = xy_res;
  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "Iy", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &Iy[0]);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  dim_1d = xy_res;
  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "ux", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &ux[0]);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  dim_1d = xy_res;
  dataspace = H5Screate_simple(1, &dim_1d, NULL);
  dataset = H5Dcreate1(h5_fid, "uy", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
    H5P_DEFAULT, &uy[0]);
  H5Dclose(dataset);
  H5Sclose(dataspace);
  
  // close the image file and allocated arrays
  
  delete_array_2d(I);
  delete_array_1d(ur);
  delete_array_1d(uphi);
  
  delete_array_1d(Ix);
  delete_array_1d(Iy);
  delete_array_1d(ux);
  delete_array_1d(uy);
  delete_array_2d(Ixy);
  
  H5Fclose(h5_fid);

  return 0;
}
