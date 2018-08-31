////////////////////////////////////////////////////////////////////////////////
//
// dump_box.cpp
//
// Created Jun 15, 2011 by HJvE
// Last modified Feb 19, 2016 by HJvE
//
// dump_box creates an hdf5 dump from a grid snapshot file
// SYNTAX: dump_box <filename> <index>
//
////////////////////////////////////////////////////////////////////////////////

#include "dump_box.h"

////////////////////////////////////////////////////////////////////////////////

c_dump_box :: c_dump_box()
{
  varname = new char[1]; // dummy declaration
}

c_dump_box :: ~c_dump_box()
{
  // memory management;
  
  delete[] varname;
}

////////////////////////////////////////////////////////////////////////////////

int c_dump_box :: dump_2D()
{
  wmax = 0.0, wmin = 1e300; // fluid quantity extrema on requested grid

  w = array_2d<double>(dumpxres, dumpyres);
  x = array_2d<double>(dumpxres, dumpyres);
  y = array_2d<double>(dumpxres, dumpyres);

  if (!quiet)
  {
    printf("# output filename = %s\n", outputfilename);
    printf("# image range: x = %e .. %e, y = %e .. %e\n", dumpxmin, dumpxmax, 
      dumpymin, dumpymax);
    printf("# resolution of dump: %dx%d\n", dumpxres, dumpyres);
    printf("# fluid variable number: %d\n", fluidvar);
    printf("#--------------------------------------------------------------\n");
    fflush(stdout);
  }
  
  for (cx = 0; cx < dumpxres; cx++)
  {
    for (cy = 0; cy < dumpyres; cy++)
    {
      x[cx][cy] = 
        dumpxmin + (cx + 0.5) * (dumpxmax - dumpxmin) / (double) dumpxres;
      y[cx][cy] = 
        dumpymin + (cy + 0.5) * (dumpymax - dumpymin) / (double) dumpyres;

      r = sqrt(x[cx][cy] * x[cx][cy] + y[cx][cy] * y[cx][cy]);
      theta = atan(x[cx][cy] / y[cx][cy]);

      cor.x = x[cx][cy];
      cor.y = y[cx][cy];
      cor.r = r;
      cor.theta = theta;

      box.set_local(cor);
      w[cx][cy] = box.local[fluidvar];

      if (w[cx][cy] > wmax) wmax = w[cx][cy];
      if (w[cx][cy] < wmin) wmin = w[cx][cy];
    } // cy
  } // cx

  printf("#wmin = %e\n", wmin);
  printf("#wmax = %e\n", wmax);

  printf("#xmin = %e\n", dumpxmin);
  printf("#xmax = %e\n", dumpxmax);
  printf("#ymin = %e\n", dumpymin);
  printf("#ymax = %e\n", dumpymax);

  printf("# t = %e\n", cor.t);

    hid_t outfile, dataset, dataspace;
    hsize_t dim_1d, dim_2d[2];
  
    outfile = H5Fcreate(outputfilename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  
    dim_1d = 1;
    dataspace = H5Screate_simple(1, &dim_1d, NULL);
  
    dataset = H5Dcreate1(outfile, "wmax", H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
      &wmax);
    H5Dclose(dataset);
    dataset = H5Dcreate1(outfile, "wmin", H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
      &wmin);
    H5Dclose(dataset);
    dataset = H5Dcreate1(outfile, "xmax", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
      &dumpxmax);
    H5Dclose(dataset);
    dataset = H5Dcreate1(outfile, "xmin", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
      &dumpxmin);
    H5Dclose(dataset);
    dataset = H5Dcreate1(outfile, "ymax", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
      &dumpymax);
    H5Dclose(dataset);
    dataset = H5Dcreate1(outfile, "ymin", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
      &dumpymin);
    H5Dclose(dataset);
    dataset = H5Dcreate1(outfile, "t", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
      &(cor.t));
    H5Dclose(dataset);
  
    H5Sclose(dataspace);
  
    dim_2d[0] = dumpxres;
    dim_2d[1] = dumpyres;
    dataspace = H5Screate_simple(2, &dim_2d[0], NULL);
  
    dataset = H5Dcreate1(outfile, "w", H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
      &w[0][0]);
    H5Dclose(dataset);

    dataset = H5Dcreate1(outfile, "x", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
      &x[0][0]);
    H5Dclose(dataset);

    dataset = H5Dcreate1(outfile, "y", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT, 
      &y[0][0]);
    H5Dclose(dataset);

    H5Sclose(dataspace);
  
    H5Fclose(outfile);


    delete_array_2d(w);
    delete_array_2d(x);
    delete_array_2d(y);

  return 0; // success
}

////////////////////////////////////////////////////////////////////////////////

int c_dump_box :: dump_slice()
{
  wmax = 0.0, wmin = 1e300; // fluid quantity extrema on requested grid

  ws = array_1d<double>(dumprres);
  rs = array_1d<double>(dumprres);

  cor.theta = theta;
    
  for (cr = 0; cr < dumprres; cr++)
  {
    rs[cr] = dumprmin + (cr + .5) * (dumprmax - dumprmin) / (double) dumprres;
    cor.r = rs[cr];
     
    box.set_local(cor);
    ws[cr] = box.local[fluidvar];
    
    if (ws[cr] > wmax) wmax = ws[cr];
    if (ws[cr] < wmin) wmin = ws[cr];
  }

  printf("# wmin = %e\n", wmin);
  printf("# wmax = %e\n", wmax);

  printf("# rmin = %e\n", dumprmin);
  printf("# rmax = %e\n", dumprmax);
  printf("# theta = %e\n", theta);

  printf("# t = %e\n", cor.t);

  for (cr = 0; cr < dumprres; cr++)
    printf("%e, %e\n", rs[cr], ws[cr]);

  fflush(stdout);

    delete_array_1d(ws);
    delete_array_1d(rs);

  return 0; // success
}

////////////////////////////////////////////////////////////////////////////////

int c_dump_box :: read_parameters(int argc, char* argv[])
{
  double time_provided = false;
  
  // get a time
  if (parse_double("-t=", argdouble, argc, argv)) 
  {
    cor.t = argdouble;
    time_provided = true;
  } 

  if (parse_int("-snapshot_time=", argint, argc, argv))
  {
    if (argint < 0 or argint >= box.tres)
    {
      printf("ERROR: illegal box snapshot number provided for time.\n");
      printf("Provided %d, should have been between 0 and %d\n", argint, 
        box.tres);
      fflush(stdout);
      return 1;
    }
    
    cor.t = box.t[argint];
    time_provided = true;
  }

  if (!time_provided)
  {
    printf("ERROR: no time provided. Use -t=<time> with time in seconds,"
      "or -snapshot_time=<no> with <no> a snapshot number.\n");
    fflush(stdout);
    return 1;
  }

  // parse user command line options
  if (parse_double("-xmin=", argdouble, argc, argv))
    dumpxmin = argdouble;
  else
    dumpxmin = 0.;
    
  if (parse_double("-xmax=", argdouble, argc, argv))
    dumpxmax = argdouble;
  else
    dumpymax = cor.t * v_light;
    
  if (parse_double("-ymin=", argdouble, argc, argv))
    dumpymin = argdouble;
  else
    dumpymin = 0.;

  if (parse_double("-ymax=", argdouble, argc, argv))
    dumpymax = argdouble;
  else
    dumpymax = cor.t * v_light;

  if (parse_int("-xres=", argint, argc, &argv[0]))
    dumpxres = argint;
  else
    dumpxres = 512;

  if (parse_int("-yres=", argint, argc, &argv[0]))
    dumpyres = argint;
  else
    dumpyres = 512;

  slice = false; // default, assume 2D dump requested
  
  // user command line options for slice (these overrule 2D options)
  if (parse("-slice", argc, argv) or parse("-s", argc, argv)) slice = true;

  if (parse_double("-rmin=", argdouble, argc, argv)) 
  {
    dumprmin = argdouble;
    slice = true;
  }
  else
    dumprmin = 0.;

  if (parse_double("-rmax=", argdouble, argc, argv))
  {
    dumprmax = argdouble;
    slice = true;
  }
  else
    dumprmax = cor.t * v_light;
  
  if (parse_int("-snapshot_rmax=", argint, argc, argv))
  {
    if (argint < 0 or argint >= box.Rres)
    {
      printf("ERROR: illegal box snapshot number provided for Rmax.\n");
      printf("Provided %d, should have been between 0 and %d\n", argint, 
        box.tres);
      fflush(stdout);
      return 1;
    }

    slice = true;
    dumprmax = box.R_max[0][argint];
  }

  if (parse_double("-fractional_rmin=", argdouble, argc, argv)) 
  {
    dumprmin = argdouble * dumprmax;
    slice = true;
  }

  if (parse_int("-rres=", argint, argc, &argv[0]))
  {
    dumprres = argint;
    slice = true;
  }
  else
    dumprres = 100;
  
  if (parse_double("-slice_theta=", argdouble, argc, argv))
  {
    theta = argdouble;
    slice = true;
  }
  else
    theta = 0.0;
  
  if (parse_int("-var=", argint, argc, &argv[0]))
    fluidvar = argint;
  else
    fluidvar = rho_;
  
  if (parse_string("-varname=", varname, argc, argv))
    fluidvar = sscanvar(varname);
  
  if (fluidvar < 0)
  {
    printf("# WARNING: wrong fluid name provided: %s. Defaulting to rho\n",
      varname);
    fluidvar = rho_;
  }

  if (fluidvar == p_)
  {
    printf("# WARNING pressure stored in BOX. Plotting eint instead.\n");
    fluidvar = eint_;
  }

  if (fluidvar == tau_)
  {
    printf("# WARNING tau stored in BOX. Plotting eint instead.\n");
    fluidvar = eint_;
  }

  if (parse("-quiet", argc, argv) or parse("-q", argc, &argv[0])) 
    quiet = true;
  else
    quiet = false;

  if (!parse_string("-output=", outputfilename, argc, &argv[0]))
    strcpy(outputfilename, "dump.h5");
  
  return 0; // return success
}

////////////////////////////////////////////////////////////////////////////////

int c_dump_box :: main(int argc, char* argv[])
{
  int result;

  //----------------------------------------------------------------------------

  outputfilename = new char[512];

  // check if filename passed as argument and if so setup the filename
  if (argc < 3 or parse("-help", argc, argv))
  {
    if (argc < 3)
      printf("Not enough arguments.\n");
      
    printf("use: \"dump_box <filename> -t=[float]\" or "
      "\"dump_box <filename> -snapshot_time=[int]\"\n");
    printf("Overview of command-line settings:\n");
    printf("  -t=[float]               BOX time in seconds, comoving frame\n");
    printf("  -snapshot_time=[int]     BOX time, taken from stored BOX times, "
      "typically ranging from 0 - 100 (inclusive)\n");
    printf("  -xmin=[float]            lower boundary Cartesian, lab frame x "
      "(cm), default is 0\n");
    printf("  -xmax=[float]            upper boundary Cartesian, lab frame x "
      "(cm), default is c * t\n");
    printf("  -ymin=[float]            lower boundary Cartesian, lab frame y "
      "(cm) ,default is 0\n");
    printf("  -ymax=[float]            upper boundary Cartesian, lab frame y "
      "(cm), default is c * t\n");
    printf("  -xres=[int]              resolution of BOX image x direction, "
      "default is 512\n");
    printf("  -yres=[int]              resolution of BOX image y direction, "
      "default is 512\n");
    printf("  -var=[int]               fluid variable to output, by number, "
      "default is 0 (rho)\n");
    printf("  -varname=[string]        fluid variable to output, by name, "
      "default is rho\n");
    printf("  -help                    print this message, and quit\n");
    printf("  -output=[string]         name of output hdf5 file. Is set to "
      "dump.h5 by default\n");
    printf("  -quiet, -q               avoid output to stdout. Disabled by"
      " default\n");
    printf("\n");
    
    printf("The following command-line settings trigger a single slice dump at"
      " fixed theta and directly to stdout\n");
    printf("  -slice, -s               Switch to 1D slice at fixed theta\n");
    printf("  -rmin=[float]            lower boundary radial r (cm), default 0"
      "\n");
    printf("  -rmax=[float]            upper boundary radial r (cm), default "
      "c*t\n");
    printf("  -snapshot_rmax=[int]     upper boundary radial r, taken from "
      "stored BOX peak radii\n");
    printf("  -fractional_rmin=[float] lower boundary radial r, as fraction of "
      "upper boundary\n");
    printf("  -rres                    resolution r direction, default 100\n");
    printf("  -slice_theta             angle theta of slice (rad), default 0"
      "\n");
    
    printf("\n");
    printf("Fluid variables by name and number:\n");
    printf("  0   rho    comoving density (gram cm^-3)\n");
    printf("  1   eint   internal energy density (erg cm^-3)\n");
    printf("  2   v_x    velocity component in generalized x-direction "
      "(i.e. r-direction), lab frame, units of c \n");
    printf("  3   v_y    velocity component in generalized y-direction "
      "(i.e. theta-direction), lab frame, units of c\n");
    printf("  4   v      magnitude of velocity, lab frame, units of c\n");
    printf("  5   lfac   Lorentz factor, lab frame\n");
    printf("  6   D      lab frame density (gram cm^-3)\n");
    printf("  9   N      comoving number density (cm^-3)\n");

    if (argc < 3)
      return 1;
    else
      return 0;
  }

  // load the box file into memory and set up global state for requested time
  box.load(argv[1]);

  result = read_parameters(argc, argv);

  if (result == 0)
  {
    box.t_e = cor.t;  
    box.set_global();

    if (slice)
      result = dump_slice();
    else
      result = dump_2D();
  }

  box.close();

  //----------------------------------------------------------------------------
  // Memory management and exit program
  
  delete[] outputfilename;

  return result;
}


////////////////////////////////////////////////////////////////////////////////
// MAIN program
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
  c_dump_box dump_box; // the main class
  
  return dump_box.main(argc, argv);
}

