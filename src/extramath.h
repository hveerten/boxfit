////////////////////////////////////////////////////////////////////////////////
//
// extramath.h
//
// additional mathematical constants and routines
// Created pre December 18, 2008 by HJvE
// Last modified June 8, 2012 by HJvE
//
////////////////////////////////////////////////////////////////////////////////

#ifndef EXTRAMATH_H_
#define EXTRAMATH_H_

#include <stdlib.h>
#include <stdio.h>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
// global variables
////////////////////////////////////////////////////////////////////////////////

extern double smallest_number;
extern double largest_number;
extern double zero;
extern double inf;

////////////////////////////////////////////////////////////////////////////////
// mathematical constants
////////////////////////////////////////////////////////////////////////////////

#define PI            3.14159265358979323846
#define invPI         0.318309886
#define inv2PI        0.159154943
#define inv4PI        0.079577471
#define onethird      0.33333333333333333333
#define onesixth      0.16666666666666666666
#define onefifteenth  0.06666666666666666666
#define twothirds     0.66666666666666666666
#define rad2deg       57.29577951
#define deg2rad       0.017453292

////////////////////////////////////////////////////////////////////////////////
// basic routines
////////////////////////////////////////////////////////////////////////////////

#define SIGN(a, b) ((b) >= 0.0 ? abs(a) : -abs(a) )

// NOTE: SWITCH TO fmax from <cmath> !
double dmax(const double &a, const double &b);
double dmin(const double &a, const double &b);
double dmin(const double &a, const double &b, const double &c );
double dmax(const double &a, const double &b, const double &c );

////////////////////////////////////////////////////////////////////////////////
// BEGIN CODE NOT PUBLIC

// check whether two variables are approximately equal
bool approx(const double a, const double b, const double epsilon);
double chebval(double x, double *c, int n); // c version of numpy function

extern unsigned long MC_bifurcation_sample;
extern unsigned long MC_sample_minimum;
extern double MC_sample_increase;

int sign(double x);
int sign(int x);
bool even(int x);
bool odd(int x);
double mod(double a, double b);
double rpow(double x, int y); // fast recursive ('r') power algorithms for non  MAKE ARGUMENTS INTO CONSTANTS!!!!!!!!!
double rpow(double x, char y); // floating point numbers only
int rpow(int x, int y);
int rpow( int x, unsigned int y );
int rpow(int x, char y);
long rpow(long x, int y);
long rpow(long x, char y);
int imin(const int &a, const int &b); // return lowest of the two arguments
int imax(const int &a, const int &b);
long lmax(const long &a, const long &b);
long lmin(const long &a, const long &b);
void swap(double &a, double &b);
double SQR(const double &x); // squares a variable

void intersect_ray_sphere( double c[3], double r, double p[3], double d[3], 
  double w[3], double &t, bool &intersection );
// calculates the point of intersection of a ray originating at p, moving in the
// direction d with a sphere centered at o and with radius r. The resulting 
// point of intersection is returned in w. If there is no intersection
// then intersection is set to false. All in Cartesian coordinates. t is the
// distance between the point p and the intersection w

void intersect_ray_cone( double mu, double p[3], double d[3], double w[3],
  double &t, bool &intersection);
// calculates the point of intersection of a ray originating at p, moving in the
// direction d (should be normalized!) with a cone along the z-axis, starting
// at the origin. The point of intersection is returned in w. All is in
// Cartesian coordinates. mu is the cosine of the maximum angle of the cone,
// which is half the opening angle. t is the distance between the point p and
// the intersection w

////////////////////////////////////////////////////////////////////////////////
// Standard integration routines.
// In 1D: Bulirsch-Stoer with Richardson extrapolation 
// in 2D: Specialized Monte-Carlo with Sobol' sequence random numbers
////////////////////////////////////////////////////////////////////////////////

template <class T> double integrate_1D( T *a_base, int integrand_ID, double x_0, 
  double x_1, double res, double *err, double h_min );
  // will call Bulirsch-Stoer with Richardson extrapolation.
  // This is built to work on class member functions. The arguments are:
  // *a_base         pointer to instance of class of which the integrand is a
  //                 member function.
  // integrand_ID    entry in the array of function pointers provided in class
  //                 a_base.
  // x_0, x_1        Integration domain.
  // res             resolution, or the maximum stepsize.
  // err             Fractional error ( abs value of error divided by result )
  //                 On output err will contain the ABSOLUTE ERROR
  // h_min           minimum stepsize allowed
  // Note: the function integrate_1D expects relay function 'relay_integrand_1D'
  // within class, which in turns calls the correct member function via a
  // pointer to class member function for function that returns double and
  // takes a single (type double) argument.

template <class T> double integrate_1D( T *a_base, int integrand_ID, double x_0, 
  double x_1, int res_ID, double *err, double h_min );
  // We have replaced the fixed value 'res' with a pointer ID to a function in
  // a_base that provides a smart, coordinate dependent, upper limit on the
  // stepsize. All again works via the relay function

template <class T> double integrate_1D_abserr( T *a_base, int integrand_ID, 
  double x_0, double x_1, double res, double *err, double h_min );
  // same as previous, only now with absolute error instead of fractional

template <class T, class I> double integrate_1D( T *a_base, int integrand_ID,
  long S_0, long S_1, I *a_ivar, int ivar_ID, int depth, long level, int 
  max_depth, double f_err_in, double *abs_err_out, long *no_set );

////////////////////////////////////////////////////////////////////////////////

template < class I, class T > class integral_1D_dd;

////////////////////////////////////////////////////////////////////////////////

template <class T> double integrate_2D_miser( T *a_base, int integrand_ID, 
  double *err, unsigned long npts_start );
  // err          on input fractional error allowed. On output absolute error
  // functions expects relay function 'relay_integrand_2D' in class a_class
  // relay_integrand_2D will refer to correct member function, using
  // integrand_ID
  // domain is taken to be 0-1 in both x and y dir

template <class T > double integrate_2D_MC( T *a_base, int integrand_ID, 
  double *err);
  // domain is taken to be 0-1 in both x and y dir
  
template <class T, class I0, class I1> double integrate_2D_MC( T *a_base, 
  int integrand_ID, double x_0, double x_1, I0 *y0_base, int y0_ID, I1 *y1_base,
  int y1_ID, double *err );

template <class T, class I0, class I1> double integrate_2D_miser( T *a_base, 
  int integrand_ID, double x_0, double x_1, I0 *y0_base, int y0_ID, I1 *y1_base,
  int y1_ID, double *err, unsigned long npts_start );

template <class T, class I0, class I1> double integrate_2D_RE( T *a_base, int
  integrand_ID, double x_0, double x_1, double x_max, double x_min, I0 *y0_base,
  int y0_ID, I1 *y1_base, int y1_1D, double y_max, double y_min, double *err );
  // x_max maximum stepsize       x direction
  // x_min minimum stepsize       x direction
    
////////////////////////////////////////////////////////////////////////////////
// random numbers
////////////////////////////////////////////////////////////////////////////////

double rnd(double max );

///////////////////////////////////////////////////////////////////////////////
// Special Functions
///////////////////////////////////////////////////////////////////////////////

double log_gamma_function(double x);

///////////////////////////////////////////////////////////////////////////////
// polynomial-solving
///////////////////////////////////////////////////////////////////////////////

double polynomial_Newton_Rhapson(double a4, double a3, double a2,
            double a1, double a0, double start, double epsilon);

///////////////////////////////////////////////////////////////////////////////
// Monte-Carlo integration
///////////////////////////////////////////////////////////////////////////////

#define MAXBIT 30

void sobseq2D(double x[], const bool& init_call); // TO DO: rewrite routine to properly start counting at 0 instead of 1
double sobseq_integrate(double (*f)(double, double), double x_b, double x_e,
  double y_b, double y_e, double acc, int min_points); 
  // min_points is minimum number of points before starting to check the error

////////////////////////////////////////////////////////////////////////////////
// integral of 2nd degree polynomial
////////////////////////////////////////////////////////////////////////////////

double integrate_polynome(double x1, double x2, double x3, double y1, double y2,
  double y3 );
  
////////////////////////////////////////////////////////////////////////////////
// polynomial interpolation
////////////////////////////////////////////////////////////////////////////////

void pol_int(double xa[], double ya[], int n, double x, double *y, double *dy);
// dy is the error estimate (absolute value, not fractional)
double pol_int( double xa[], double ya[], int n, double x, double *dy );
  // same as previous routine, only result returned in different way
double pol_int( double xa[], double ya[], int n, double *dy );
  // if x not given, extrapolation to x=0.0 is assumed

///////////////////////////////////////////////////////////////////////////////
// Runge-Kutta algorithms
///////////////////////////////////////////////////////////////////////////////

double RungeKutta_fixed_step(double (*dydx)(double,double), double y_b, double 
  x_b, double x_e, long n_init);
// solves 1st order differential equation, starting from (x_b,y_b), using 4th
// order Runge-Kutta algorithm with fixed stepsize. This routine is build for 
// speed, hence there is no error checking.

double RungeKutta_adaptive_step(double (*dydx)(double), double y_b, 
  double x_b, double x_e, long n_init, int levels, double error);

double RungeKutta_adaptive_step(double (*dydx)(double,double), double y_b, 
  double x_b, double x_e, long n_init, int levels, double error);
// solves 1st order differential equation, starting from (x_b,y_b), using 4th
// order Runge-Kutta algorithm with adaptive stepsize. The routine is actually
// fifth order, because of local interpolation, but with up to 4th order error 
// check. The parameter 'error' is the allowed error on the final result. It 
// gets scaled to an error per step. If this error is exceeded, the stepsize is 
// halved. This is done at most 'levels' times. If the error is still exceeded 
// after this, the program logs an error message.

double RungeKutta_adaptive_step(double (*dydx)(double, double, double),
  double y_b, double x_b, double x_e, double z, long n_init, int levels, 
  double error);
// Same as previous routine, only allows for functions of three parameters, the
// third parameter being an independent parameter assumed constant over the
// integration domain.

void mmid( double y, double dydx, double xs, double htot, int nstep, double 
  &yout, double (*deriv)( double ) );
////////////////////////////////////////////////////////////////////////////////
// Special Functions
////////////////////////////////////////////////////////////////////////////////

double gammafunc(double x);
// Returns gamma function of x

void bessik(double x, double xnu, double *ri, double *rk, double *rip, double
  *rkp);
// returns modified Bessel function, ri = I_nu, rk = K_nu and their derivatives
// with rip being I'_nu amd rkp being K'_nu

double chebev(double a, double b, double c[], int m, double x);
void beschb(double x, double *gam1, double *gam2, double *gamp1, double *gammi);

////////////////////////////////////////////////////////////////////////////////

// the following templates are not currently in use, since all integrations are
// now done using gsl routines. The templates are commented temporarily

/*
template <class T> double modmid(T *a_base, int integrand_ID, double x_0, 
  double x_1, double dydx0, double dydx1, int n_steps );

template <class T> double integrate_1D_step( T *a_base, int integrand_ID, double
  x, double *stepsize, double steperr, double *abserr, double dydx0, 
  double dydx1, double h_min );

template <class T> double integrate_1D( T *a_base, int integrand_ID, double x_0, 
  double x_1, double res, double *err, double h_min );

template <class T> double integrate_1D( T *a_base, int integrand_ID, double x_0, 
  double x_1, int res_ID, double *err, double h_min );

template <class T> double integrate_1D_step_abserr( T *a_base, int integrand_ID, 
  double x, double *stepsize, double steperr, double *abserr, double dydx0, 
  double dydx1, double h_min );

template <class T> double integrate_1D_abserr( T *a_base, int integrand_ID, 
  double x_0, double x_1, double res, double *err, double h_min );

template < class I, class T > class integral_1D_dd;

// there are more templates in extramath_templates.h that have no prior 
// declaration included above. The prior declarations are not necessary, but
// omitting them will cause the intel c compiler to complain unneccessarily 

#include "extramath_templates.h" */

// END CODE NOT PUBLIC

#endif /*EXTRAMATH_H_*/
