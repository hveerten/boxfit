////////////////////////////////////////////////////////////////////////////////
//
// box.h
//
// Created Jun 15, 2011 by HJvE
// Last modified July 29, 2016 by HJvE
//
// defines box structures and classes
//
////////////////////////////////////////////////////////////////////////////////

#include "box.h"

////////////////////////////////////////////////////////////////////////////////
// global variables

extern int noderank; // identity of single core within node, 0 if no MPI
#if OPEN_MPI_ == ENABLED_
  // these variables should have been declared as global variables in the main
  // program
  extern int nodesize;
  extern MPI_Comm nodecom;
#endif

////////////////////////////////////////////////////////////////////////////////
// static member declarations

c_box* c_multibox :: box;

////////////////////////////////////////////////////////////////////////////////

c_box :: c_box()
{
  BOX_enabled = true; // use BOX data at times covered by BOX
  BM_enabled = true; // use BM data. If BOX_enabled also TRUE, only use BM data
    // at times before those covered by BOX
  jet = BOTH_; // by default, use both jet and counterjet

  box_allocated = false; // no memory allocated on startup
  box_file_open = false; // no files have been opened on startup
  
  BOX_filename_allocated = false; // no filename string stored on startup
}

////////////////////////////////////////////////////////////////////////////////

c_box :: ~c_box()
{
  close(); // call the memory de-allocation routines again, in case of
    // unexpected shutdown
}

////////////////////////////////////////////////////////////////////////////////

void c_box :: load(const char* filename)
{
  hid_t dataset, dataspace;
  
  #if OPEN_MPI_ == ENABLED_
    MPI_Aint size;
    MPI_Aint localsize; // only noderank 0 gets to allocate the full array
    int size_disp;
  #endif
  
  #if BOOST_ == ENABLED_
    int i_counterjet; // counterjet flag. 1 or 0 in file, stored as bool in code 
  #endif

  //----------------------------------------------------------------------------
  // open the file
  
  h5_fid = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (h5_fid < 0)
  {
    printf("(c_box :: load) Error. Could not open file %s\n",filename);
    fflush(stdout);
    exit(1);
  }
  else
    box_file_open = true;

  BOX_filename = new char[strlen(filename) + 1];
  BOX_filename_allocated = true;
  sprintf(BOX_filename, "%s", filename);

  //----------------------------------------------------------------------------
  // Doublecheck consistency BOOST_ compiler setting with BOX file data
  
  #if BOOST_ == ENABLED_
  
    dataset = H5Dopen1(h5_fid, "boost");
    if (dataset < 0)
    {
      printf("(c_box :: load) Error. BOX file %s has no boost, while "
        "program compiled for boosted frames.\n", filename);
      fflush(stdout);
      exit(1);
    }
    else
    {
      H5Dclose(dataset);
    }

  #else

    if (H5Lexists(h5_fid, "boost", H5P_DEFAULT))
    {
      printf("(c_box :: load) Error. BOX file %s has a boost, while program "
        "compiled without boosted frame functionality.\n", filename);
      fflush(stdout);
      exit(1);
    }

  #endif

  //----------------------------------------------------------------------------
  // get the physics parameters

  dataset = H5Dopen1(h5_fid, "lfac_initial");
  dataspace = H5Dget_space(dataset);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
    &lfac_initial);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  dataset = H5Dopen1(h5_fid, "E");
  dataspace = H5Dget_space(dataset);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
    &E_actual);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  dataset = H5Dopen1(h5_fid, "n");
  dataspace = H5Dget_space(dataset);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
    &n_actual);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  dataset = H5Dopen1(h5_fid, "theta_0");
  dataspace = H5Dget_space(dataset);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
    &theta_0);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  if (H5Lexists(h5_fid, "k", H5P_DEFAULT))
  {
  dataset = H5Dopen1(h5_fid, "k");
  dataspace = H5Dget_space(dataset);
  H5Dread(dataset, H5T_NATIVE_INT, dataspace, dataspace, H5P_DEFAULT, 
    &k);
  H5Dclose(dataset);
  H5Sclose(dataspace);
  BM.k = k;
  }
  else // BOX version < 2
  {
    BM.k = 0.;
    k = 0.;
  }
  
  if (H5Lexists(h5_fid, "boxversion", H5P_DEFAULT))
  {
  dataset = H5Dopen1(h5_fid, "boxversion");
  dataspace = H5Dget_space(dataset);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
    &version);
  H5Dclose(dataset);
  H5Sclose(dataspace);
  }
  else // BOX version < 2
  {
    version = 1.;
  }
  
  
  #if BOOST_ == ENABLED_
  
    dataset = H5Dopen1(h5_fid, "boost");
    dataspace = H5Dget_space(dataset);
    H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
      &boost);
    H5Dclose(dataset);
    H5Sclose(dataspace);

    dataset = H5Dopen1(h5_fid, "counterjet");
    dataspace = H5Dget_space(dataset);
    H5Dread(dataset, H5T_NATIVE_INT, dataspace, dataspace, H5P_DEFAULT, 
      &i_counterjet);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    if (i_counterjet == 1)
      counterjet = true;
    else
      counterjet = false;
    
  #endif

  // at first set requested energy and density to default values
  E = E_actual;
  n = n_actual;
  Sr = 1.; // scale factor unity
  
  //----------------------------------------------------------------------------
  // get the box resolution
  
  dataset = H5Dopen1(h5_fid, "tres");
  dataspace = H5Dget_space(dataset);
  H5Dread(dataset, H5T_NATIVE_INT, dataspace, dataspace, H5P_DEFAULT, &tres);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  dataset = H5Dopen1(h5_fid, "thetares");
  dataspace = H5Dget_space(dataset);
  H5Dread(dataset, H5T_NATIVE_INT, dataspace, dataspace, H5P_DEFAULT, 
    &thetares);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  dataset = H5Dopen1(h5_fid, "Rres");
  dataspace = H5Dget_space(dataset);
  H5Dread(dataset, H5T_NATIVE_INT, dataspace, dataspace, H5P_DEFAULT, 
    &Rres);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  //----------------------------------------------------------------------------
  // assign memory, based on box resolution

  #if OPEN_MPI_ == ENABLED_
    if (noderank == 0)
    {
      size = Rres * thetares * tres;
      localsize = size;
    }
    else
      localsize = 0;

    // allocate shared window
    if (H5Lexists(h5_fid, "pres", H5P_DEFAULT))
      MPI_Win_allocate_shared(localsize*sizeof(double), sizeof(double),
        MPI_INFO_NULL, nodecom, &pres, &preswin);
    MPI_Win_allocate_shared(localsize*sizeof(double), sizeof(double),
      MPI_INFO_NULL, nodecom, &dens, &denswin);
    MPI_Win_allocate_shared(localsize*sizeof(double), sizeof(double),
      MPI_INFO_NULL, nodecom, &eint, &eintwin);
    MPI_Win_allocate_shared(localsize*sizeof(double), sizeof(double),
      MPI_INFO_NULL, nodecom, &velr, &velrwin);
    MPI_Win_allocate_shared(localsize*sizeof(double), sizeof(double),
      MPI_INFO_NULL, nodecom, &veltheta, &velthetawin);

    MPI_Win_allocate_shared(localsize*sizeof(double), sizeof(double),
      MPI_INFO_NULL, nodecom, &r, &rwin);
    MPI_Win_allocate_shared(localsize*sizeof(double), sizeof(double),
      MPI_INFO_NULL, nodecom, &dr, &drwin);

    // tell cores not of top rank in their respective node where to find data

    if (noderank != 0)
    {
      if (H5Lexists(h5_fid, "pres", H5P_DEFAULT))
        MPI_Win_shared_query(preswin, 0, &size, &size_disp, &pres);
  
      MPI_Win_shared_query(denswin, 0, &size, &size_disp, &dens);
      MPI_Win_shared_query(eintwin, 0, &size, &size_disp, &eint);
      MPI_Win_shared_query(velrwin, 0, &size, &size_disp, &velr);
      MPI_Win_shared_query(velthetawin, 0, &size, &size_disp, &veltheta);

      MPI_Win_shared_query(rwin, 0, &size, &size_disp, &r);
      MPI_Win_shared_query(drwin, 0, &size, &size_disp, &dr);
    }

  #endif
  #if OPEN_MPI_ == DISABLED_
    if (H5Lexists(h5_fid, "pres", H5P_DEFAULT))
      pres = new double[Rres * thetares * tres];
    
    dens = new double[Rres * thetares * tres];
    eint = new double[Rres * thetares * tres];
    velr = new double[Rres * thetares * tres];
    veltheta = new double[Rres * thetares * tres];
  
    r = new double[Rres * thetares * tres];
    dr = new double[Rres * thetares * tres];
  #endif


  // TO DO: ERROR CHECKING WHETHER THE RIGHT MPI CAPABILITIES ARE THERE
  // TO DO: CONFIRM NO MEMORY LEAKS WITH DEBUGGING TOOLS

  
  theta_max = array_1d<double>(tres); // outer theta boundaries at each time
  t = array_1d<double>(tres);
  R_max = array_2d<double>(thetares, tres);
  R_peak = array_2d<double>(thetares, tres);

  // assign scaling & interpolation related
  r_scale0 = array_1d<double>(thetares);
  r_scale1 = array_1d<double>(thetares);
  
  #if BOOST_ == ENABLED_
  
    if (counterjet)
    {
      // if there is a counterjet, assign memory
      
      #if OPEN_MPI_ == DISABLED_
        pres_ctr = new double[Rres * thetares * tres];
        dens_ctr = new double[Rres * thetares * tres];
        eint_ctr = new double[Rres * thetares * tres];
        velr_ctr = new double[Rres * thetares * tres];
        veltheta_ctr = new double[Rres * thetares * tres];

        r_ctr = new double[Rres * thetares * tres];
        dr_ctr = new double[Rres * thetares * tres];
      #endif
      #if OPEN_MPI_ == ENABLED_

        if (noderank == 0)
        {
          size = Rres * thetares * tres;
          localsize = size;
        }
        else
          localsize = 0;

        // allocate shared window
        MPI_Win_allocate_shared(localsize*sizeof(double), sizeof(double),
          MPI_INFO_NULL, nodecom, &pres_ctr, &pres_ctrwin);
        MPI_Win_allocate_shared(localsize*sizeof(double), sizeof(double),
          MPI_INFO_NULL, nodecom, &dens_ctr, &dens_ctrwin);
        MPI_Win_allocate_shared(localsize*sizeof(double), sizeof(double),
          MPI_INFO_NULL, nodecom, &eint_ctr, &eint_ctrwin);
        MPI_Win_allocate_shared(localsize*sizeof(double), sizeof(double),
          MPI_INFO_NULL, nodecom, &velr_ctr, &velr_ctrwin);
        MPI_Win_allocate_shared(localsize*sizeof(double), sizeof(double),
          MPI_INFO_NULL, nodecom, &veltheta_ctr, &veltheta_ctrwin);

        MPI_Win_allocate_shared(localsize*sizeof(double), sizeof(double),
          MPI_INFO_NULL, nodecom, &r_ctr, &r_ctrwin);
        MPI_Win_allocate_shared(localsize*sizeof(double), sizeof(double),
          MPI_INFO_NULL, nodecom, &dr_ctr, &dr_ctrwin);

        // tell cores not of top rank in their node where to find data

        if (noderank != 0)
        {
          MPI_Win_shared_query(pres_ctrwin, 0, &size, &size_disp, &pres_ctr);
          MPI_Win_shared_query(dens_ctrwin, 0, &size, &size_disp, &dens_ctr);
          MPI_Win_shared_query(eint_ctrwin, 0, &size, &size_disp, &eint_ctr);
          MPI_Win_shared_query(velr_ctrwin, 0, &size, &size_disp, &velr_ctr);
          MPI_Win_shared_query(veltheta_ctrwin, 0, &size, &size_disp,
            &veltheta_ctr);

          MPI_Win_shared_query(r_ctrwin, 0, &size, &size_disp, &r_ctr);
          MPI_Win_shared_query(dr_ctrwin, 0, &size, &size_disp, &dr_ctr);
        }

      #endif // OPEN_MPI_ == ENABLED_

      theta_max_ctr = array_1d<double>(tres);
      R_max_ctr = array_2d<double>(thetares, tres);
      R_peak_ctr = array_2d<double>(thetares, tres);
    }
    
  #endif
  
  box_allocated = true;

  //----------------------------------------------------------------------------
  // load grid times

  dataset = H5Dopen1(h5_fid, "t");
  dataspace = H5Dget_space(dataset);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
    &t[0]);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  //----------------------------------------------------------------------------
  // load theta_max
  
  dataset = H5Dopen1(h5_fid, "theta_max");
  dataspace = H5Dget_space(dataset);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
    &theta_max[0]);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  //----------------------------------------------------------------------------
  // load R_max
  
  dataset = H5Dopen1(h5_fid, "R_max");
  dataspace = H5Dget_space(dataset);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
    &R_max[0][0]);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  //----------------------------------------------------------------------------
  // load R_peak
  
  dataset = H5Dopen1(h5_fid, "R_peak");
  dataspace = H5Dget_space(dataset);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
    &R_peak[0][0]);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  //----------------------------------------------------------------------------
  
  if (noderank == 0)
  {
    // load dr
  
    dataset = H5Dopen1(h5_fid, "dr");
    dataspace = H5Dget_space(dataset);
    H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, dr);
    H5Dclose(dataset);
    H5Sclose(dataspace);

    //--------------------------------------------------------------------------
    // load r
  
    dataset = H5Dopen1(h5_fid, "r");
    dataspace = H5Dget_space(dataset);
    H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, r);
    H5Dclose(dataset);
    H5Sclose(dataspace);
  }

  //----------------------------------------------------------------------------
  
  #if BOOST_ == ENABLED_

    if (counterjet)
    {
      // load theta_max
  
      dataset = H5Dopen1(h5_fid, "theta_max_ctr");
      dataspace = H5Dget_space(dataset);
      H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
        &theta_max_ctr[0]);
      H5Dclose(dataset);
      H5Sclose(dataspace);

      // load R_max
  
      dataset = H5Dopen1(h5_fid, "R_max_ctr");
      dataspace = H5Dget_space(dataset);
      H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
        &R_max_ctr[0][0]);
      H5Dclose(dataset);
      H5Sclose(dataspace);

      // load R_peak
  
      dataset = H5Dopen1(h5_fid, "R_peak_ctr");
      dataspace = H5Dget_space(dataset);
      H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
        &R_peak_ctr[0][0]);
      H5Dclose(dataset);
      H5Sclose(dataspace);

      if (noderank == 0)
      {
        // load dr
  
        dataset = H5Dopen1(h5_fid, "dr_ctr");
        dataspace = H5Dget_space(dataset);
        H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
          dr_ctr);
        H5Dclose(dataset);
        H5Sclose(dataspace);

        // load r
  
        dataset = H5Dopen1(h5_fid, "r_ctr");
        dataspace = H5Dget_space(dataset);
        H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
          r_ctr);
        H5Dclose(dataset);
        H5Sclose(dataspace);
      }
    }

  #endif // BOOST_ == ENABLED_

  //----------------------------------------------------------------------------
  // load eint
    
  dataset = H5Dopen1(h5_fid, "eint");
  dataspace = H5Dget_space(dataset);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, eint);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  //----------------------------------------------------------------------------
  
  if (noderank == 0)
  {
    // load dens (comoving frame density)

    dataset = H5Dopen1(h5_fid, "dens");
    dataspace = H5Dget_space(dataset);
    H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, dens);
    H5Dclose(dataset);
    H5Sclose(dataspace);
  
    //--------------------------------------------------------------------------
    // load v_r
  
    dataset = H5Dopen1(h5_fid, "velr");
    dataspace = H5Dget_space(dataset);
    H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, velr);
    H5Dclose(dataset);
    H5Sclose(dataspace);

    //--------------------------------------------------------------------------
  // load v_theta
  
    dataset = H5Dopen1(h5_fid, "veltheta");
    dataspace = H5Dget_space(dataset);
    H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
      veltheta);
    H5Dclose(dataset);
    H5Sclose(dataspace);

    //--------------------------------------------------------------------------
    // load pres
  
    if (H5Lexists(h5_fid, "pres", H5P_DEFAULT))
    {
      dataset = H5Dopen1(h5_fid, "pres");
      dataspace = H5Dget_space(dataset);
      H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT,
        pres); 
      H5Dclose(dataset);
      H5Sclose(dataspace);
    }
  }

  //----------------------------------------------------------------------------
  
  #if BOOST_ == ENABLED_
  
    if (counterjet and noderank == 0)
    {
      // load pres
      
      dataset = H5Dopen1(h5_fid, "pres_ctr");
      dataspace = H5Dget_space(dataset);
      H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
        pres_ctr);
      H5Dclose(dataset);
      H5Sclose(dataspace);

      // load eint
    
      dataset = H5Dopen1(h5_fid, "eint_ctr");
      dataspace = H5Dget_space(dataset);
      H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
        eint_ctr);
      H5Dclose(dataset);
      H5Sclose(dataspace);

      // load dens (comoving frame density)
  
      dataset = H5Dopen1(h5_fid, "dens_ctr");
      dataspace = H5Dget_space(dataset);
      H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
        dens_ctr);
      H5Dclose(dataset);
      H5Sclose(dataspace);

      // load v_r
  
      dataset = H5Dopen1(h5_fid, "velr_ctr");
      dataspace = H5Dget_space(dataset);
      H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
        velr_ctr);
      H5Dclose(dataset);
      H5Sclose(dataspace);

      // load v_theta
  
      dataset = H5Dopen1(h5_fid, "veltheta_ctr");
      dataspace = H5Dget_space(dataset);
      H5Dread(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
        veltheta_ctr);
      H5Dclose(dataset);
      H5Sclose(dataspace);
    }
    
  #endif // BOOST_ == ENABLED_
}

////////////////////////////////////////////////////////////////////////////////

void c_box :: close()
{
  if (box_allocated)
  {
    // release memory
    #if OPEN_MPI_ == DISABLED_
      if (version >= 2.) delete[] pres;
      delete[] dens;
      delete[] eint;
      delete[] velr;
      delete[] veltheta;
      delete[] r;
      delete[] dr;
    #endif
    #if OPEN_MPI_ == ENABLED_
      if (version >= 1.999) MPI_Win_free(&preswin);
      MPI_Win_free(&denswin);
      MPI_Win_free(&eintwin);
      MPI_Win_free(&velrwin);
      MPI_Win_free(&velthetawin);
      MPI_Win_free(&rwin);
      MPI_Win_free(&drwin);
    #endif

    delete_array_1d(t);
  
    delete_array_1d(theta_max);
    delete_array_2d(R_max);
    delete_array_2d(R_peak);
  
    delete_array_1d(r_scale0);
    delete_array_1d(r_scale1);
 
    #if BOOST_ == ENABLED_
  
      if (counterjet)
      {
        #if OPEN_MPI_ == DISABLED_
          delete[] pres_ctr;
          delete[] dens_ctr;
          delete[] eint_ctr;
          delete[] velr_ctr;
          delete[] veltheta_ctr;
          delete[] r_ctr;
          delete[] dr_ctr;
        #else
          if (version >= 1.999) MPI_Win_free(&pres_ctrwin);
          MPI_Win_free(&dens_ctrwin);
          MPI_Win_free(&eint_ctrwin);
          MPI_Win_free(&velr_ctrwin);
          MPI_Win_free(&veltheta_ctrwin);
          MPI_Win_free(&r_ctrwin);
          MPI_Win_free(&dr_ctrwin);
        #endif
        
        delete_array_1d(theta_max_ctr);
        delete_array_2d(R_max_ctr);
        delete_array_2d(R_peak_ctr);
      }
        
    #endif
    
    box_allocated = false;
  }
  
  // close file
  if (box_file_open)
  {
    H5Fclose(h5_fid);
    box_file_open = false;
  }
  
  if (BOX_filename_allocated)
  {
    delete[] BOX_filename;
    BOX_filename_allocated = false;
  }
}

////////////////////////////////////////////////////////////////////////////////

double c_box :: get_var(int vartype, int i_r, int i_theta, int i_t)
{
  double aux, x;
  double Sn;
  int i = i_t + i_r * thetares * tres + i_theta * tres; // flattened index
  
  // compute rescale factor for mass densities (pressure and energy density
  // scale in the same way, since cm / s is not affected by the rescaling).
  if (k == 0)  
    Sn = n / n_actual; 
  else // assuming wind (k=2) as only alternative
    Sn = pow(n / n_actual, 3) / pow(E / E_actual, 2);
  
  switch( vartype )
  {
    case eint_: 
      if (i_r < 0) return 1e-10 * n * m_p * v_light * v_light;
      return eint[i] * Sn;
    case rho_: 
      if (i_r < 0) return n * m_p;
      return dens[i] * Sn;
    case p_:
      if (version < 2.) return -1.;
      if (i_r < 0) return 1e-10 * n * m_p * v_light * v_light;
      return pres[i] * Sn;    
    case v_x_: 
      if (i_r < 0) return 0.0;
      return velr[i];
    case v_y_: 
      if (i_r < 0) return 0.0;
      return veltheta[i];
    case lfac_: 
      if (i_r < 0) return 1.0;
      aux = (1.0 - velr[i] * velr[i] - veltheta[i] * veltheta[i]);
      if (aux < 0.0) aux = 1e-20;
      x = sqrt(1.0 / aux);
      return x;
    default:
      printf("c_box :: get_var ERROR. unknown variable\n");
      fflush(stdout);
      abort();
  }
  
  return -1.0;
}

////////////////////////////////////////////////////////////////////////////////

#if BOOST_ == ENABLED_

  double c_box :: get_var_ctr(int vartype, int i_r, int i_theta, int i_t)
  {
    // same as c_box :: get_var, only for the case where the counter jet is
    // stored in separate BOX data.
    double aux, x;
    double Sn;
    int i = i_t + i_r * thetares * tres + i_theta * tres; // flattened index
  
    if (k == 0)  
      Sn = n / n_actual;
    else // assuming wind as only alternative
      Sn = pow(n / n_actual, 3) / pow(E / E_actual, 2);
  
    switch( vartype )
    {
      case eint_: 
        if (i_r < 0) return 1e-10 * n * m_p * v_light * v_light;
        return eint_ctr[i] * Sn;
      case rho_: 
        if (i_r < 0) return n * m_p;
        return dens_ctr[i] * Sn;
      case p_:
        if (i_r < 0) return 1e-10 * n * m_p * v_light * v_light;
        return pres_ctr[i] * Sn;    
      case v_x_: 
        if (i_r < 0) return 0.0;
        return velr_ctr[i];
      case v_y_: 
        if (i_r < 0) return 0.0;
        return veltheta_ctr[i];
      case lfac_: 
        if (i_r < 0) return 1.0;
        aux = (1. - velr_ctr[i]*velr_ctr[i] - veltheta_ctr[i]*veltheta_ctr[i]);
        if (aux < 0.0) aux = 1e-20;
          x = sqrt(1.0 / aux);
          return x;
        default:
          printf("c_box_ctr :: get_var ERROR. unknown variable\n");
          fflush(stdout);
          exit(1);
    }
  
    return -1.0;
  }

#endif // BOOST_ == ENABLED_

////////////////////////////////////////////////////////////////////////////////

void c_box :: find_cell_mix(double a_r, double a_theta, int &cr, int &ctheta)
{
  double a_r_base = a_r / Sr; // radial coordinate scaled to baseline value
  double R_max_cur; // time-interpolated R_max, different for each theta
  bool found; // true if correct radial index number found
  int cL, cM, cR; // auxiliary variables for computing radial index number
  double r0, r1; // lower and upper interpolated radial cell boundaries
  int i_cur, i_next; // flattened indices
  
  // find theta index first
  if (a_theta >= theta_max_cur)
  {
    cr = -1;
    ctheta = -1;
    return;
  }
  ctheta = min((int) (a_theta / dtheta_cur), thetares - 1);

  // check if not outside of blast wave at this angle
  R_max_cur = R_max[ctheta][i_t_cur + 1] * fract + 
    (1. - fract) * R_max[ctheta][i_t_cur];
  
  if (a_r_base >= R_max_cur)
  {
    cr = -1;
    ctheta = -1;
    return;
  }

  // determine appropriate r entry, first check boundaries
  found = false;
  cL = 0; cR = Rres - 1;
  
  i_cur = i_t_cur + cL * thetares * tres + ctheta * tres;
  i_next = i_t_cur + 1 + cL * thetares * tres + ctheta * tres;
  
  r1 = (r[i_cur] + dr[i_cur]) * (1. - fract) + (r[i_next] + dr[i_next]) * fract;
  if (r1 > a_r_base)
  {
    cr = 0; found = true;
  }
  else
  {
    i_cur = i_t_cur + cR * thetares * tres + ctheta * tres;
    i_next = i_t_cur + 1 + cR * thetares * tres + ctheta * tres;
    r0 = r[i_cur] * (1. - fract) + 
      r[i_next] * fract;
      
    if (r0 <= a_r_base)
    {
      cr = cR; found = true;
    }
  }

  while (!found)
  {
    cM = (cL + cR) / 2;

    i_cur = i_t_cur + cM * thetares * tres + ctheta * tres;
    i_next = i_t_cur + 1 + cM * thetares * tres + ctheta * tres;

    r0 = r[i_cur] * (1. - fract) + r[i_next] * fract;

    r1 = (r[i_cur] + dr[i_cur]) * (1.- fract) + 
      (r[i_next] + dr[i_next]) * fract;

    if (r1 > a_r_base and r0 <= a_r_base)
    {
      found = true;
      cr = cM;
    }
    else
    {
      if (r0 <= a_r_base)
        cL = cM;
      else
        cR = cM;
    }
  } 
}

////////////////////////////////////////////////////////////////////////////////

#if BOOST_ == ENABLED_

  void c_box :: find_cell_mix_ctr(double a_r, double a_theta, int &cr,
    int &ctheta)
  {
    double a_r_base = a_r / Sr; // radial coordinate scaled to baseline value
    double R_max_cur_ctr; // time-interpolated R_max, different for each theta
    bool found; // true if correct radial index number found
    int cL, cM, cR; // auxiliary variables for computing radial index number
    double r0, r1; // lower and upper interpolated radial cell boundaries
    int i_cur, i_next; // flattened indices
  
    // flip angle to counterjet convention
    a_theta = M_PI - a_theta;
  
    // find theta index first
    if (a_theta >= theta_max_cur_ctr)
    {
      cr = -1;
      ctheta = -1;
      return;
    }
    ctheta = min((int) (a_theta / dtheta_cur_ctr), thetares - 1);
  
    // check if not outside of blast wave at this angle
    R_max_cur_ctr = R_max_ctr[ctheta][i_t_cur + 1] * fract + 
      (1. - fract) * R_max_ctr[ctheta][i_t_cur];
  
    if (a_r_base >= R_max_cur_ctr)
    {
      cr = -1;
      ctheta = -1;
      return;
    }
  
    // determine appropriate r entry, first check boundaries
    found = false;
    cL = 0; cR = Rres - 1;

    i_cur = i_t_cur + cL * thetares * tres + ctheta * tres;
    i_next = i_t_cur + 1 + cL * thetares * tres + ctheta * tres;
  
    r1 = (r_ctr[i_cur] + dr_ctr[i_cur]) * 
      (1. - fract) + 
      (r_ctr[i_next] + dr_ctr[i_next]) * 
      fract;
    if (r1 > a_r_base)
    {
      cr = 0; found = true;
    }
    else
    {
      i_cur = i_t_cur + cR * thetares * tres + ctheta * tres;
      i_next = i_t_cur + 1 + cR * thetares * tres + ctheta * tres;
   
      r0 = r_ctr[i_cur] * (1. - fract) + 
        r_ctr[i_next] * fract;
      
      if (r0 <= a_r_base)
      {
        cr = cR; found = true;
      }
    }

    while (!found)
    {
      cM = (cL + cR) / 2;

      i_cur = i_t_cur + cM * thetares * tres + ctheta * tres;
      i_next = i_t_cur + 1 + cM * thetares * tres + ctheta * tres;
    
      r0 = r_ctr[i_cur] * (1. - fract) + 
        r_ctr[i_next] * fract;

      r1 = (r_ctr[i_cur] + dr_ctr[i_cur]) * 
        (1. - fract) + 
        (r_ctr[i_next] + dr_ctr[i_next]) *
        fract;
    
      if (r1 > a_r_base and r0 <= a_r_base)
      {
        found = true;
        cr = cM;
      }
      else
      {
        if (r0 <= a_r_base)
          cL = cM;
        else
          cR = cM;
      }
    } 
  }

#endif

////////////////////////////////////////////////////////////////////////////////

void c_box :: set_local_empty(s_coordinates &cor, double *a_local)
{
  a_local[rho_] = n * m_p * pow(cor.r / 1e17, BM.k);
  a_local[eint_] = 1e-10 * a_local[rho_] * v_light * v_light;
  a_local[v_x_] = 0.0;
  a_local[v_y_] = 0.0;
  a_local[lfac_] = 1.0;
  a_local[N_] = a_local[rho_] * invm_p;
  a_local[D_] = a_local[lfac_] * a_local[rho_];
  a_local[v_] = sqrt(a_local[v_x_] * a_local[v_x_] + a_local[v_y_] * 
    a_local[v_y_]);  
}

////////////////////////////////////////////////////////////////////////////////

void c_box :: set_cell_coords(int ir, int ith, double &a_r, double &a_th)
{
  // normal case
  double dtheta;

  int i; // flattened index
  i = i_t_cur + ir * thetares * tres + ith * tres;
  a_r = (r[i] + 0.5 * dr[i]) * Sr;
  dtheta = theta_max[i_t_cur] / (double) thetares;
  a_th = ((double) ith + 0.5) * dtheta;
  return;
}

////////////////////////////////////////////////////////////////////////////////

#if BOOST_ == ENABLED_

  void c_box :: set_cell_coords_ctr(int ir, int ith, double &a_r, double &a_th)
  {
    // normal case
    double dtheta;

    int i; // flattened index
    i = i_t_cur + ir * thetares * tres + ith * tres;

    a_r = (r_ctr[i] + 0.5 * dr_ctr[i]) * Sr;
    dtheta = theta_max[i_t_cur] / (double) thetares;
    a_th = M_PI - ((double) ith + 0.5) * dtheta;
      // note that counterjet theta is zero pointing along counterjet axis
    return;
  }

#endif // BOOST_ == ENABLED_

////////////////////////////////////////////////////////////////////////////////

void c_box :: set_local(s_coordinates cor)
// This is the one set_local routine that can be called from outside (i.e. not
// protected). At this point, it is not clear whether the coordinates lie in the
// jet or counter-jet (if enabled)
{
  int i_r, i_theta; // box indices
  int i_cur; // flattened index lower bracketing time
  int i_next; // flattened index one time step up
  
  bool use_BM = false;
  bool use_BOX = false;
  bool use_empty = false; // empty space result if outside of jet in BM solution

  //----------------------------------------------------------------------------

  // determine local fluid state, using time interpolation for BOX snapshots

  // determine relevant time indices, or alternatively provide the BM solution
  // in case a time before the first box-covered time is requested
  int offset = 2; // 2
  if (BOX_enabled and t_e > t[offset] * Sr and t_e <= t[tres - 1] * Sr) 
    use_BOX = true;
  if (BM_enabled and !BOX_enabled) use_BM = true;
  if (BM_enabled and t_e <= t[offset] * Sr) use_BM = true;

  //----------------------------------------------------------------------------

  if (use_BM)
  {
    //printf("# [useBM]\n"); fflush(stdout);
    // time too early to be covered by box, or BOX disabled. 
    // Return BM value instead

    // return empty space in case a value outside of the jet(s) is probed
    if (cor.theta > theta_0 and jet == FORWARD_)
      use_empty = true;
    else if (M_PI - cor.theta > theta_0 and jet == RECEDING_)
      use_empty = true;
    else if (cor.theta > theta_0 and M_PI - cor.theta > theta_0 and 
      jet == BOTH_)
      use_empty = true;

    if (use_empty)
    {
      set_local_empty(cor, local);
      return;
    }

    // In the case of a boosted frame, there is no single lab frame time
    // across the grid, so we first need to determine the local lab frame time
    #if BOOST_ == ENABLED_
      BM.t = cor.t;
      BM.set_global_from_time(cor.t);
    #endif

    // definitely return empty once BM no longer applicable
    if (BM.lfac_shock < 1.2)
    {
      set_local_empty(cor, local);
      return;
    }

    BM.set_local(cor.r); // set local coordinates, includes empty space for
      // values outside shock radius.

    local[lfac_] = BM.lfac;
    local[rho_] = BM.rho;
    local[eint_] = BM.e_therm;  

    local[v_x_] = BM.beta; // velocity in the radial direction, units of c
    local[v_y_] = 0.0;
    local[v_] = local[v_x_]; // absolute value of velocity in units of c

    local[N_] = local[rho_] * invm_p;
    local[D_] = local[lfac_] * local[rho_];

    return;
  }

  //----------------------------------------------------------------------------

  if (use_BOX)
  {
    #if BOOST_ == DISABLED_
    
      // check for counter jet angles
      if (cor.theta > M_PI * 0.5)
      {
        if (jet != FORWARD_) // i.e. counter jet enabled
        {
          cor.theta = M_PI - cor.theta; // flip angle
        }
        else
        {
          set_local_empty(cor, local);
        }
      }
    
    #else
    
      if (cor.theta > M_PI * 0.5)
      {
        if (jet != FORWARD_ and counterjet) // i.e. counter jet enabled
        {
          
          // find counter jet cell indices
          find_cell_mix_ctr(cor.r, cor.theta, i_r, i_theta);
          if (i_r < 0)
          {
            set_local_empty(cor, local);
            return;
          }

          i_cur = i_t_cur + i_r * thetares * tres + i_theta * tres;
          i_next = i_t_cur + 1 + i_r * thetares * tres + i_theta * tres;

          // set local variables, interpolated values drawn from snapshots
          local[rho_] = (1. - fract) * dens_ctr[i_cur] +
            fract * dens_ctr[i_next];
          local[rho_] *= Sn;

          local[eint_] = (1. - fract) * eint_ctr[i_cur] +
            fract * eint_ctr[i_next];
          local[eint_] *= Sn;

          local[v_x_] = (1. - fract) * velr_ctr[i_cur] +
            fract * velr_ctr[i_next];

          local[v_y_] = (1. - fract) * veltheta_ctr[i_cur] +
            fract * veltheta_ctr[i_next];

          // set local variables, derived values
          local[lfac_] = sqrt(max(1.0, 1.0 / (1.0 - local[v_x_] * local[v_x_] -
            local[v_y_] * local[v_y_])));

          local[N_] = local[rho_] * invm_p;
          local[D_] = local[lfac_] * local[rho_];
          local[v_] = sqrt(local[v_x_] * local[v_x_] + 
            local[v_y_] * local[v_y_]);

          return;
          
        }
        else
        {
          set_local_empty(cor, local);
        }
      }
      else // i.e. cor.theta <= m_PI * 0.5
      {
        if (jet == RECEDING_)
        {
          set_local_empty(cor, local);
          return;
        }
      }

    #endif
    
    // At this point we are probing either the forward jet in boosted or non-
    // boosted frame, or the counter jet in a non-boosted frame

    find_cell_mix(cor.r, cor.theta, i_r, i_theta);
    if (i_r < 0)
    {
      set_local_empty(cor, local);
      return;
    }

    i_cur = i_t_cur + i_r * thetares * tres + i_theta * tres;
    i_next = i_t_cur + 1 + i_r * thetares * tres + i_theta * tres;
    
    // set local variables, interpolated values drawn from snapshots
    local[rho_] = (1. - fract) * dens[i_cur] + fract * dens[i_next];
    local[rho_] *= Sn;

    local[eint_] = (1. - fract) * eint[i_cur] + fract * eint[i_next];
    local[eint_] *= Sn;

    local[v_x_] = (1. - fract) * velr[i_cur] + fract * velr[i_next];

    local[v_y_] = (1. - fract) * veltheta[i_cur] + fract * veltheta[i_next];

    // set local variables, derived values
    local[lfac_] = sqrt(max(1.0, 1.0 / (1.0 - local[v_x_] * local[v_x_] -
      local[v_y_] * local[v_y_])));

    local[N_] = local[rho_] * invm_p;
    local[D_] = local[lfac_] * local[rho_];
    local[v_] = sqrt(local[v_x_] * local[v_x_] + local[v_y_] * local[v_y_]);

    return;
  }

  //----------------------------------------------------------------------------

  if (!use_BM and !use_BOX) // uncovered emission time
  {
    set_local_empty(cor, local);
  }
}

////////////////////////////////////////////////////////////////////////////////

double c_box :: get_var(int vartype, s_coordinates cor)
{
  // NOT YET IMPLEMENTED
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

void c_box :: set_scale_factors()
{
  // scaling settings
  if (k == 0)
  {
    Sr = pow((E * n_actual) / (E_actual * n), 1. / 3.);
    Sn = n / n_actual;
    BM.k = 0.;
  }
  else // assuming wind-type environment is the only alternative
  {
    Sr = (E * n_actual) / (E_actual * n);
    Sn = pow(n / n_actual, 3) / pow(E / E_actual, 2);
    BM.k = 2.;
  }
}

////////////////////////////////////////////////////////////////////////////////

void c_box :: set_global()
{
  set_scale_factors();
  
  // analytical solution settings
  BM.E = E;
  BM.n_ext = n;

  // if a Lorentz boost is applied to the grid, there is no single lab frame 
  // time, so setting BM globablly makes no sense
  #if BOOST_ != ENABLED_
    BM.set_global_from_time(t_e);
  #endif

  for (i_t_cur = 1; i_t_cur < tres - 1; i_t_cur++)
  {
    if (t[i_t_cur] * Sr > t_e) break;
  }
  i_t_cur--;

  // interpolated theta_max, due to time interpolation. Note that, in the
  // boosted frame case, time interpolation uses boosted times, not lab times
  fract = (t_e - t[i_t_cur] * Sr) / 
    (t[i_t_cur + 1] * Sr - t[i_t_cur] * Sr);

  theta_max_cur = fract * theta_max[i_t_cur + 1] + 
    (1.0 - fract) * theta_max[i_t_cur];
  dtheta_cur = theta_max_cur / (double) thetares;

  #if BOOST_ == ENABLED_
    if (counterjet)
    {
      theta_max_cur_ctr = fract * theta_max_ctr[i_t_cur + 1] + 
        (1.0 - fract) * theta_max_ctr[i_t_cur];
      dtheta_cur_ctr = theta_max_cur_ctr / (double) thetares;
    }
  #endif

  // at early times or for BM only use theta_0 instead of BOX provided theta_max
  if (t_e < t[0] * Sr or !BOX_enabled)
  {
    theta_max_cur = theta_0;
    
    #if BOOST_ == ENABLED_
      if (counterjet) theta_max_cur_ctr = theta_0;
    #endif
  }
}

////////////////////////////////////////////////////////////////////////////////

c_multibox :: c_multibox()
{
  no_boxes = -1; // default setting, as of yet unspecified number of boxes
  initialized = false;
}

////////////////////////////////////////////////////////////////////////////////

void c_multibox :: initialize()
{
  if (!initialized and no_boxes > 0)
  {
    box = new c_box[no_boxes];
    initialized = true;
  }
  else
  {
    if (initialized)
      printf("ERROR c_multibox: trying to initialize twice.\n");
    if (no_boxes <= 0)
      printf("ERROR c_multibox: no boxes to initialize\n");
    fflush(stdout);
    exit(1);
  }
}

////////////////////////////////////////////////////////////////////////////////

c_multibox :: ~c_multibox()
{
  if (initialized) delete[] box;
}

////////////////////////////////////////////////////////////////////////////////
// returns scaled max radius of boxes, so assumes that E, n are currently set

double c_multibox :: get_largest_radius()
{
  int i_th0;
  double R = 0;
  double Sr; // scale factor
  
  for (i_th0 = 0; i_th0 < no_boxes; i_th0++)
  {
    if (box[0].k == 0)
      Sr = pow((E * box[i_th0].n_actual) / (box[i_th0].E_actual / n), 1. / 3.);
    else // assuming wind the only alternative
      Sr = (E * box[i_th0].n_actual) / (box[i_th0].E_actual / n);
    
    R = max(R, box[i_th0].R_max[0][box[i_th0].tres - 1] * Sr);
  }
  
  return R;
}

////////////////////////////////////////////////////////////////////////////////
// Set up source for fluid information for all boxes: BM, BOX or BOTH

void c_multibox :: set_BOXBM(bool BOX_setting, bool BM_setting)
{
  int i_th;

  for (i_th = 0; i_th < no_boxes; i_th++)
  {
    box[i_th].BOX_enabled = BOX_setting;
    box[i_th].BM_enabled = BM_setting;
  }
}

////////////////////////////////////////////////////////////////////////////////
// prepare multibox for probing a given initial opening angle.

void c_multibox :: set_theta_0(double a_th0)
{
  th0 = a_th0;
  // stay on safe end of roundoff errors, which might keep us inside BOX range:
  if (a_th0 > box[no_boxes - 1].theta_0) th0 = th0 - 1e-10;
  if (a_th0 < box[0].theta_0) th0 = th0 + 1e-10;

  // determine boxes w opening angles surrounding the requested angle
  for (i_th0 = 1; i_th0 < no_boxes - 1; i_th0++)
    if (box[i_th0].theta_0 > th0) break;
  i_th0--;

  // determine theta_0 interpolation fraction
  fractheta0 = (th0 - box[i_th0].theta_0)
    / (box[i_th0 + 1].theta_0 - box[i_th0].theta_0);

  if (th0 > box[no_boxes - 1].theta_0)
    i_th0 = -1; // not covered
  
  // pass along energy and density settings to relevant boxes
  box[i_th0].E = E;
  box[i_th0 + 1].E = E;
  box[i_th0].n = n;
  box[i_th0 + 1].n = n;
}

////////////////////////////////////////////////////////////////////////////////
