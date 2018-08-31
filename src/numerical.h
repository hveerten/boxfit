#ifndef NUMERICAL_H_
#define NUMERICAL_H_

////////////////////////////////////////////////////////////////////////////////
//
// numerical.h
//
// Created August 31 2011, by HJvE
// Last Modified August 31 2011, by HJvE
//
// General numerical routines, references can be found in Press et al.
//
////////////////////////////////////////////////////////////////////////////////

#include "environment.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "extramath.h"
#include "arraytools.h"
#if OPEN_MPI_ == ENABLED_
#include "mpi.h"
#endif // OPEN_MPI_ == ENABLED_

////////////////////////////////////////////////////////////////////////////////

class c_random
{
  public:
    c_random();
    ~c_random();
    void initialize(long a_seed);
    
    // random number from uniform distribution. Returns double between 0 and 1,
    // excluding the boundaries
    double local(); // get random number for local core.
    double global(); // get same random number on all cores.
    
    // random number from Gaussian distribution with width 1.0 and mean 0.0
    double localGaussian();  
    double globalGaussian();

    long seed; // initial seed value provided on initialization
    long currentseed; // current seed value
    long entryseed; // current entry number
    
  protected:
  
    bool initialized;
  
    // Park & Miller suggestions for initial values Lehmer generator
    long a, m, q, r;
    double scale;

    // We use Bays-Durham shuffle table to break serial correlations
    long* BaysDurham;
    int BDlength, BDscale;
};

////////////////////////////////////////////////////////////////////////////////

class c_simplex
{
  public:
  
    int n; // number of variables
    double (*p_f)(double *); // pointer to function to be minimized
    double **x; // set of n + 1 simplex coordinates
    double *x_centroid, **x_trial;
    double *x_max, *x_min; // edges of parameter space
    int max_iter; // single temperature maximum number of iterations
    double temp; // annealing temperature temperature
    double alpha, beta, gamma; // control parameters for downhill simplex
    int i_h, i_l; // entry labels highest and lowest y values

    double *y; // function evaluation results at simplex coordinates
    double *y_trial;
    double *y_pert, *y_trial_pert; // temperature perturbed y, y_trial

    c_random *p_random; // pointer to external random number generator

    c_simplex();
    ~c_simplex();
    void initialize(int a_n);
    void clear_memory();  
    void evaluate(int i); // calculate y for single entry x_i in simplex list
    void minimize(); // perform simplex function mimimalization
    
    void load(const char* filename);
    void save(const char* filename);

  protected:
    
    bool initialized; // initialization flag, for memory management
    double invn; // equal to 1.0 / n
    double invnplus1;

    void reflect();
    void expand();
    void contract();
    void contract_all();
    void replace(int i_h, int i);
    void set_centroid();
};

#endif // NUMERICAL_H_
