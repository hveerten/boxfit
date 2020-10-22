////////////////////////////////////////////////////////////////////////////////
//
// numerical.cpp
//
// Created August 31 2011, by HJvE
// Last Modified August 31 2011, by HJvE
//
// General numerical routines, references can be found in Press et al. '86, '92
//
////////////////////////////////////////////////////////////////////////////////

#include "numerical.h"

////////////////////////////////////////////////////////////////////////////////

c_random :: c_random()
{
  initialized = false;

  // Park & Miller suggestions for initial values Lehmer generator
  a = 16807; // multiplier
  m = 2147483647; // modules, a large prime number
  q = 127773; // q = m / a
  r = 2836; // r = m % a
  scale = 4.656612875245797e-10; // is 1.0 / m

  BDlength = 32;
  BDscale = 1 + (m - 1) / BDlength; 
}

////////////////////////////////////////////////////////////////////////////////

c_random :: ~c_random()
{
  if (initialized) delete[] BaysDurham;
}

////////////////////////////////////////////////////////////////////////////////

void c_random :: initialize(long a_seed)
{
  int i;
  long hi, lo; // auxiliary values, using Schrage's method to avoid int overflow

  if (!initialized)
  {
    seed = a_seed; // for reference, store the initial seed value
    currentseed = a_seed; // the current seed keeps changing
  
    // initialize Bays-Durham shuffle table.
    BaysDurham = new long[BDlength];
  
    // pre-run
    for (i = 0; i < 8; i++)
    {
      hi = currentseed / q;
      lo = currentseed - hi * q; // equivalent to currentseed % q
      currentseed = a * lo - r * hi;
      if (currentseed < 0) currentseed += m;
    }
  
    // fill the BaysDurham shuffle table.
    for (i = 0; i < BDlength; i++)
    {
      hi = currentseed / q;
      lo = currentseed - hi * q; // faster equivalent to currentseed % q
      currentseed = a * lo - r * hi;
      if (currentseed < 0) currentseed += m;
      BaysDurham[BDlength - i - 1] = currentseed; // we count down to allow for
        // direct comparison against random number routines from Press et al '86
    }
    
    entryseed = BaysDurham[0];
    initialized = true;
  }
}

////////////////////////////////////////////////////////////////////////////////

double c_random :: local()
{
  int i; // Bays-Durham table entry
  double result; // random number between 0..1 (excluding boundaries)
  long hi, lo;
  
  // scale entry seed random number down from 0 .. m (excluding boundaries) to
  // 0 .. BDlength - 1 (including boundaries), to get random table entry
  i = entryseed / BDscale;
  
  // update the entry seed
  entryseed = BaysDurham[i];
  
  // calculate the result from the updated entry seed
  result = scale * entryseed;
  if (result > 1.0 - 1.0e-14) result = 1.0 - 1e-14;  // exclude upper boundary

  // restock the table entry with an updated current seed value
  hi = currentseed / q;
  lo = currentseed - hi * q; // equivalent to currentseed % q, but faster
  currentseed = a * lo - r * hi;
  if (currentseed < 0) currentseed += m;
  BaysDurham[i] = currentseed;
  
  return result;  
}

////////////////////////////////////////////////////////////////////////////////

double c_random :: localGaussian()
{
  double a1, a2, asqr, a3; // auxiliary variables
  static double spare;          // routine calculates two random numbers, so
  static bool spareset = false; // we keep one for next function call
  
  if (!spareset)
  {
    // get two random numbers in the unit circle
    do
    {
      a1 = 2.0 * local() - 1.0;
      a2 = 2.0 * local() - 1.0;
  
      asqr = a1 * a1 + a2 * a2;
    } while (asqr >= 1.0 or asqr <= 0.0);
    a3 = sqrt(-2.0 * log(asqr) / asqr);
    
    spare = a1 * a3;
    spareset = true;
    return a2 * a3;
  }
  else
  {
    spareset = false;
    return spare;
  }
}  

////////////////////////////////////////////////////////////////////////////////

double c_random :: global()
{
  #if OPEN_MPI_ == DISABLED_
    return local();
  #endif // OPEN_MPI_ == DISABLED_
  
  #if OPEN_MPI_ == ENABLED_
  
    double x;
    
    if (myid == 0) x = local();
          
    MPI_Bcast(&x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    return x;

  #endif // PARALELL_ == ENABLED_  
}

////////////////////////////////////////////////////////////////////////////////

double c_random :: globalGaussian()
{
  #if OPEN_MPI_ == DISABLED_
    return localGaussian();
  #endif // OPEN_MPI_ == DISABLED_
  
  #if OPEN_MPI_ == ENABLED_
  
    double x;
    
    if (myid == 0) x = localGaussian();
          
    MPI_Bcast(&x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    return x;

  #endif // PARALELL_ == ENABLED_  
}

////////////////////////////////////////////////////////////////////////////////

c_simplex :: c_simplex()
{
  initialized = false;
  max_iter = 500;
  alpha = 1.0;
  beta = 0.5;
  gamma = 2.0;
  temp = 100.0;
}

////////////////////////////////////////////////////////////////////////////////

c_simplex :: ~c_simplex()
{
  if (initialized) clear_memory();
}

////////////////////////////////////////////////////////////////////////////////

void c_simplex :: initialize(int a_n)
{
  if (!initialized)
  {
    n = a_n;
    invn = 1.0 / (double) n;
    invnplus1 = 1.0 / (double) (n + 1);

    x = array_2d<double>(a_n + 1, a_n);
    x_centroid = array_1d<double>(a_n);
    x_trial = array_2d<double>(2, a_n);

    x_max = array_1d<double>(a_n);
    x_min = array_1d<double>(a_n);

    y_trial = array_1d<double>(2);
    y_trial_pert = array_1d<double>(2);

    y = array_1d<double>(a_n + 1);
    y_pert = array_1d<double>(a_n + 1);
    
    initialized = true;
  }
  else
  {
    printf("c_simplex :: initialize ERROR: double initialization!\n");
    fflush(stdout);
    abort();
  }
}

////////////////////////////////////////////////////////////////////////////////

void c_simplex :: clear_memory()
{
  if (initialized)
  {
    delete_array_2d(x);
    delete_array_1d(x_centroid);
    delete_array_2d(x_trial);
    delete_array_1d(y_trial);
    delete_array_1d(y_trial_pert);
    delete_array_1d(y_pert);
    delete_array_1d(x_max);
    delete_array_1d(x_min);
    delete_array_1d(y);
    initialized = false;
  }
}

////////////////////////////////////////////////////////////////////////////////

void c_simplex :: evaluate(int i)
{
  y[i] = (*p_f)(x[i]);
}

////////////////////////////////////////////////////////////////////////////////

void c_simplex :: set_centroid()
{
  int i_var, i_sim;
  
  // determine centroid
  for (i_var = 0; i_var < n; i_var++) x_centroid[i_var] = 0.0;
   
  for (i_sim = 0; i_sim <= n; i_sim++)
  {
    if (i_sim != i_h) // exlude highest value from centroid
    {
      for (i_var = 0; i_var < n; i_var++)
      {
        x_centroid[i_var] += x[i_sim][i_var];
      }
    }
  }
  
  // normalize
  for (i_var = 0; i_var < n; i_var++)
    x_centroid[i_var] = x_centroid[i_var] * invn;
}

////////////////////////////////////////////////////////////////////////////////

void c_simplex :: reflect()
{
  int i_var;
  
  y_trial[0] = 0.0;
  for (i_var = 0; i_var < n; i_var++)
  {
    x_trial[0][i_var] = (1 + alpha) * x_centroid[i_var] - alpha * x[i_h][i_var];
    if (x_trial[0][i_var] >= x_max[i_var]) y_trial[0] = 1e300; // illegal move
    if (x_trial[0][i_var] <= x_min[i_var]) y_trial[0] = 1e300;
  }  
  if (y_trial[0] > 0) { y_trial_pert[0] = 1e300; return; } // we're done
  y_trial[0] = (*p_f)(x_trial[0]);
  y_trial_pert[0] = y_trial[0] + temp * log(p_random->global());
  
}

////////////////////////////////////////////////////////////////////////////////

void c_simplex :: expand()
{
  int i_var;
  
  y_trial[1] = 0.0;
  for (i_var = 0; i_var < n; i_var++)
  {
    x_trial[1][i_var] = (1.0 + gamma) * x_trial[0][i_var] - 
      gamma * x_centroid[i_var];
    
    if (x_trial[1][i_var] >= x_max[i_var]) y_trial[1] = 1e300;
    if (x_trial[1][i_var] <= x_min[i_var]) y_trial[1] = 1e300;
  }
  if (y_trial[1] > 0.0) { y_trial_pert[1] = 1e300; return; }
  y_trial[1] = (*p_f)(x_trial[1]);
  y_trial_pert[1] = y_trial[1] + temp * log(p_random->global());
}

////////////////////////////////////////////////////////////////////////////////

void c_simplex :: replace(int i_h, int i)
{
  int i_var;
  
  for (i_var = 0; i_var < n; i_var++) // replace worst value by 2nd trial
    x[i_h][i_var] = x_trial[i][i_var];
  y[i_h] = y_trial[i];
}

////////////////////////////////////////////////////////////////////////////////

void c_simplex :: contract()
{
  int i_var;
  
  y_trial[1] = 0.0;
  for (i_var = 0; i_var < n; i_var++)
  {
    x_trial[1][i_var] = beta * x[i_h][i_var] + (1.0 - beta) * x_centroid[i_var];
    
    if (x_trial[1][i_var] >= x_max[i_var]) y_trial[1] = 1e300;
    
    if (x_trial[1][i_var] <= x_min[i_var]) y_trial[1] = 1e300;
  }

  if (y_trial[1] > 0.0) { y_trial_pert[1] = 1e300; return; }
  y_trial[1] = (*p_f)(x_trial[1]);
  y_trial_pert[1] = y_trial[1] + temp * log(p_random->global());
}

////////////////////////////////////////////////////////////////////////////////

void c_simplex :: contract_all()
{
  int i_sim, i_var;
  
  for (i_sim = 0; i_sim <= n; i_sim++) // loop over x's
  {
    if (i_sim != i_l)
    {
      for (i_var = 0; i_var < n; i_var++)
        x[i_sim][i_var] = 0.5 * (x[i_sim][i_var] + x[i_l][i_var]);
      y[i_sim] = (*p_f)(x[i_sim]);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void c_simplex :: minimize()
{
  // iterate up to max_iter times in order to minimize the chi^2 values in the
  // simplex.

  bool done = false;
  int i = 0;
  int i_sim, i_var; // dummy indices summation over simplex, fit variables
  bool worst; // auxiliary flag to check if 1st trial value still worst value
  double y_avg; // average y value for simplex. Used in stopping criterion
  double sigma; 
  double epsilon = 1e-8; // convergence criterion

  //----------------------------------------------------------------------------

  i_h = 0; i_l = 0;

  // assumes initial simplex is set up. 

  // start minimalization
  do
  {
    // set new perturbations
    for (i_sim = 0; i_sim <= n; i_sim++)
      y_pert[i_sim] = y[i_sim] - temp * log(p_random->global());    
    
    // find highest and lowest
    for (i_sim = 0; i_sim <= n; i_sim++)
    {
      if (y_pert[i_sim] > y_pert[i_h]) i_h = i_sim;
      if (y_pert[i_sim] < y_pert[i_l]) i_l = i_sim;
    }

    set_centroid();
    
    // apply reflection, calculate first trial
    reflect();
    
    if (y_trial_pert[0] < y_pert[i_l]) 
    {  // trial results in improvement beyond old best value
      // after succesful reflection, go for expansion
      expand();
      
      if (y_trial_pert[1] < y_pert[i_l]) // still better than old best value?
      {
        replace(i_h, 1);
      }
      else // 2nd trial not beter than old best value
      {
        replace(i_h, 0);
      }
    }
    else // first trial did not result in improvement over old best value
    {
      // check if the new trial is still the worst of the bunch
      worst = true;
      for (i_var = 0; i_var < n; i_var++)
        if (i_var != i_h and y_trial_pert[0] <= y_pert[i_var]) worst = false;
      
      if (worst) // first trial still the worst
      {
        // 1st trial at least better than previous worst?
        if (y_trial_pert[0] <= y_pert[i_h]) 
        { 
          replace(i_h, 0);
          y_pert[i_h] = y_trial_pert[0];
        }
        
        contract();
        
        // 2nd trial worse than worst / 1st trial?
        if (y_trial_pert[1] > y_pert[i_h]) 
        {
          contract_all();
        }
        else // 2nd trial better than worst / 1st trial
        {
          replace(i_h, 1);
        }
      }
      else // first trial not new best value, but beats at least one other
      { 
        replace(i_h, 0);
      }
    }

    // test stopping criteria, using unperturbed values
    i++; 
    if (i >= max_iter) done = true;

    y_avg = 0.0;
    for (i_sim = 0; i_sim <= n; i_sim++)
      y_avg += y[i_sim];
    y_avg = y_avg * invnplus1;

    sigma = 0.0;
    for (i_sim = 0; i_sim <= n; i_sim++)
      sigma += (y[i_sim] - y_avg) * (y[i_sim] - y_avg);
    sigma = sigma * invnplus1;

    if (sigma < epsilon * epsilon) done = true;

  } while(!done);

  // find best and worst x
  for (i_sim = 0; i_sim <= n; i_sim++)
  {
    if (y[i_sim] > y[i_h]) i_h = i_sim;
    if (y[i_sim] < y[i_l]) i_l = i_sim;
  }
  
}

////////////////////////////////////////////////////////////////////////////////

void c_simplex :: load(const char* filename)
{
  // the load routine is set up to mirror the save simplex routine. Because
  // additional information is saved along with the simplex coordinates, some
  // redundant information is loaded from disc: the annealing temperature and
  // the chi^2 values for the simplex are not used.

  FILE *p_file;
  int i_sim, i_var;
  int err;

  p_file = fopen(filename, "r");
  err = fscanf(p_file, "# current annealing temperature = %le\n", &temp); // not used
  if (err == 0)
    { printf("READ ERROR c_simplex_load()\n"); fflush(stdout); abort(); }
  err = fscanf(p_file, "# number of fit parameters = %d\n", &n);
  if (err == 0)
    { printf("READ ERROR c_simplex_load()\n"); fflush(stdout); abort(); }

  for (i_sim = 0; i_sim <= n; i_sim++)
  {
    for (i_var = 0; i_var < n; i_var++)
    {
      err = fscanf(p_file, "%le, ", &x[i_sim][i_var]);
      if (err == 0)
       { printf("READ ERROR c_simplex_load()\n"); fflush(stdout); abort(); }
    }
    err = fscanf(p_file, "%le\n", &y[i_sim]); // not used
    if (err == 0)
      { printf("READ ERROR c_simplex_load()\n"); fflush(stdout); abort(); }
  }

  fclose(p_file);
}

////////////////////////////////////////////////////////////////////////////////

void c_simplex :: save(const char* filename)
{
  FILE *p_file;
  int i_sim, i_var;
  
  p_file = fopen(filename, "w");
  fprintf(p_file, "# current annealing temperature = %e\n", temp);
  fprintf(p_file, "# number of fit parameters = %d\n", n);
  for (i_sim = 0; i_sim <= n; i_sim++)
  {
    for (i_var = 0; i_var < n; i_var++)
      fprintf(p_file, "%e, ", x[i_sim][i_var]);
    fprintf(p_file, "%e\n", y[i_sim]);
  }

  fclose(p_file);
}
