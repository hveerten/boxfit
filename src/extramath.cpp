////////////////////////////////////////////////////////////////////////////////
//
// extramath.h
//
// additional mathematical constants and routines
// Created pre December 18, 2008 by HJvE
// Last modified September 7, 2011 by HJvE
//
////////////////////////////////////////////////////////////////////////////////

#include "extramath.h"

////////////////////////////////////////////////////////////////////////////////
// global variables
////////////////////////////////////////////////////////////////////////////////

double smallest_number = 1.0e-300;
double largest_number = 1.0e300;
double zero = 0.0;
double inf = 1.0 / zero;

////////////////////////////////////////////////////////////////////////////////

double dmax(const double &a, const double &b)
{ if (a > b) return a; else return b; }

////////////////////////////////////////////////////////////////////////////////

double dmin(const double &a, const double &b)
{ if (a <= b) return a; else return b; }

////////////////////////////////////////////////////////////////////////////////

double dmin(const double &a, const double &b, const double &c )
{
  double d = dmin(a, b);
  return dmin(d, c);
}

////////////////////////////////////////////////////////////////////////////////

double dmax(const double &a, const double &b, const double &c )
{
  double d = dmax(a, b);
  return dmax(d, c);
}

////////////////////////////////////////////////////////////////////////////////
// BEGIN CODE NOT PUBLIC

bool approx(const double a, const double b, const double epsilon)
{
  // check if two variables a, b are equal within accuracy epsilon
  if ((a - b) * (a - b) > fabs(epsilon * epsilon * a * b)) return false;
  return true;
}

////////////////////////////////////////////////////////////////////////////////

double chebval(double x, double *c, int n) // c version of numpy function
{
  double c0, c1, x2, tmp; // aux variables for calculating return
  int i;
  
  if (n == 1)
  {
    c0 = c[0];
    c1 = 0.;
  }
  else if (n == 2)
  {
    c0 = c[0];
    c1 = c[1];
  }
  else
  {
    x2 = 2. * x;
    c0 = c[n - 2];
    c1 = c[n - 1];
    for (i = 3; i <= n; i++)
    {
      tmp = c0;
      c0 = c[n - i] - c1;
      c1 = tmp + c1 * x2;
    }
  }
  printf("FUNCTION chebval (extramath.cpp) NOT YET TESTED\n");
  return c0 + c1 * x;
}

////////////////////////////////////////////////////////////////////////////////


int sign(double x)
{
  if (x<0.0) return -1; else return 1;
}

int sign(int x)
{
  if (x<0) return -1; else return 1;
}

double mod( double a, double b)
{
  return a / b - (double) ((int)(a/b));
}

bool even(int x)
{	return !(x & 1); }

bool odd(int x)
{ return (x & 1); }

double rpow(double x, int y)
{  if (y==0) return 1; else return x * rpow(x, y-1); }

double rpow(double x, char y)
{  if (y==0) return 1; else return x * rpow(x, y-1); }

int rpow(int x, int y)
{  if (y==0) return 1; else return x * rpow(x, y-1); }

int rpow(int x, unsigned int y)
{  if (y==0) return 1; else return x * rpow(x, y-1); }

int rpow(int x, char y)
{  if (y==0) return 1; else return x * rpow(x, y-1); }

long rpow(long x, int y)
{  if (y==0) return 1; else return x * rpow(x, y-1); }

long rpow(long x, char y)
{  if (y==0) return 1; else return x * rpow(x, y-1); }

int imin(const int &a, const int &b)
{ if (a <= b) return a; else return b; }

int imax(const int &a, const int &b)
{ if (a > b) return a; else return b; }

long lmax(const long &a, const long &b )
{ if (a > b) return a; else return b; }

long lmin(const long &a, const long &b)
{ if (a <= b) return a; else return b; }

void swap(double &a, double &b)
{ double c; c = a; a = b; b = c; }

double SQR(const double &x)
{ return x*x; }

////////////////////////////////////////////////////////////////////////////////

void intersect_ray_sphere( double c[3], double r, double p[3], double d[3], 
  double w[3], double &t, bool &intersection )
// calculates to point of intersection of a ray originating at p, moving in the
// direction d with a sphere centered at c and with radius r. The resulting 
// point of intersection is returned in w. If there is no intersection
// then intersection is set to false. All in Cartesian coordinates
{
  double oc[3]; // difference vector between origin sphere and ray start
  double loc_sqrd;
  double tca; // the point p + tca * d is the closest to c
  
  // assuming |d| to be unity
  oc[0] = p[0] - c[0]; oc[1] = p[1] - c[1]; oc[2] = p[2] - c[2];
  double D, doc;
  doc = d[0] * oc[0] + d[1] * oc[1] + d[2] * oc[2];
  loc_sqrd = oc[0]*oc[0]+oc[1]*oc[1]+oc[2]*oc[2];
  D = doc*doc - loc_sqrd + r*r;

  if (D == 0) // ray touches sphere, doesn't cross surface 
  {
		tca = - doc;
		if ( tca < 0 )
    { 
			w[0] = 0.0; w[1] = 0.0; w[2] = 0.0; t =-1.0; 
      intersection = false; return; 
    }
    w[0] = p[0] + tca*d[0]; w[1] = p[1] + tca*d[1]; w[2] = p[2] + tca*d[2];
    t = tca;
    intersection = true;
    return;
  }
  
  if (D < 0) // no intersection at all
  { 
  	w[0] = 0.0; w[1] = 0.0; w[2] = 0.0; t =-1.0; 
    intersection = false; return; 
  }
  
  // final option D > 0
  D = sqrt(D); // we need the square root of discriminant
  tca = -doc - D;
  if ( tca < 0 ) tca = -doc + D; // make sure we get closest intersection 
  if ( tca < 0 ) // both intersection point in wrong direction
  {
  	w[0] = 0.0; w[1] = 0.0; w[2] = 0.0; t =-1.0; 
    intersection = false; return; 
  }
  w[0] = p[0] + tca*d[0]; w[1] = p[1] + tca*d[1]; w[2] = p[2] + tca*d[2];
  intersection = true;
  t = tca;
  return;
}

////////////////////////////////////////////////////////////////////////////////

void intersect_ray_cone( double mu, double p[3], double d[3], double w[3],
  double &t, bool &intersection)
{
	double A, B, C, D; // abc-formula components
	double t1, t2;
	t = 1e300;///////////////// TEMP?//////////////////////////////////////////////////////////
	
	
	A = mu*mu - d[2]*d[2];
	B = ( p[0]*d[0] + p[1]*d[1] + p[2]*d[2] ) * mu * mu - p[2] * d[2];
	C = ( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] ) * mu * mu - p[2] * p[2];
	D = B*B - A*C; // discriminant
	
	if ( D < 0.0 ) // no intersection at all
	{
		w[0] = 0.0; w[1] = 0.0; w[2] = 0.0; t = -1.0;
	  intersection = false; return; 
  }
  
  D = sqrt(D);
	
	t1 = ( -B + D ) / A; t2 = ( -B - D ) / A;
	
	if ( t1 < 0 and t2 < 0 ) // no intersection in this ray direction
	{	
		w[0] = 0.0; w[1] = 0.0; w[2] = 0.0; t = -1.0;
	  intersection = false; return; 
	}
  
	if ( t1 < 0 and t2 >= 0 ) t = t2;
	if ( t2 < 0 and t1 >= 0 ) t = t1;
	if ( t1 >= 0 and t2 >= 0 ) t = dmin( t1, t2 ); 
	  // further one intersects the counterjet instead
	  
	if ( t > 1e200 )
	{
		printf("wtf?\n");
		printf("A = %e, B = %e, C= %e, D = %e\n", A, B, C, D );
		fflush(stdout); abort();
  }
		
	intersection = true;
	w[0] = p[0] + t * d[0]; w[1] = p[1] + t * d[1]; w[2] = p[2] + t * d[2];
	
	// important: at the moment 'true' is returned also when only the counterjet
	// is intersected. This could be checked by comparing 't' to the distance
	// to the origin.
}

////////////////////////////////////////////////////////////////////////////////
// Standard integration routines
////////////////////////////////////////////////////////////////////////////////

void mmid( double y, double dydx, double xs, double htot, int nstep, double 
  &yout, double (*deriv)( double ) )
{
  int n;
  double x, swap, h2, h, ym, yn;
  
  // set stepsize
  h = htot / nstep;
  
  // take first step
  ym = y;
  yn = y + h * dydx;
  x = xs + h;
  yout = (*deriv)( x ); // temporary store derivative in yout
  
  // loop over steps between first and last step
  h2 = 2.0 * h;
  for (n=2; n<= nstep; n++)
  {
    swap = ym + h2 * yout;
    ym = yn;
    yn = swap;
    x += h;
    yout = (*deriv)( x );
  }
  
  // final step, yout now contain the output
  yout = 0.5 * (ym + yn + h * yout );
}

///////////////////////////////////////////////////////////////////////////////
// random numbers
///////////////////////////////////////////////////////////////////////////////

double rnd(double max )
{
  return ((double) (rand() % 1000 ) + (double (rand() % 1000) ) / 1000.0 ) 
    * max * 1e-3;
}

///////////////////////////////////////////////////////////////////////////////
// Special Functions
///////////////////////////////////////////////////////////////////////////////

double log_gamma_function(double x)
{
	double z,y,tmp,ser;
	static double cof[6]={76.18009172947146, -86.50532032941677,
		24.01409824083091, -1.231739572450155, 0.1208650973866179e-2,
		-0.5395239384953e-5 };
	int j;
	
	y=z=x;
	tmp=z+5.5;
	tmp -= (z+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

////////////////////////////////////////////////////////////////////////////////
// Root finding
////////////////////////////////////////////////////////////////////////////////

double polynomial_Newton_Rhapson(double a4, double a3, double a2,
            double a1, double a0, double start, double epsilon)
// warning: the routine assumes that f'(x) = 0 doesn't occur!			 
{
  double r,ru; // 'result' and 'result update'
  double f, fprime;
  int n=0;

  r = start;
  ru = 2*start;
  
  while((fabs(ru - r) > epsilon) && (n<15)) 
    // 15 is a pretty low number, for speed's sake
  {
    r = ru;
    f = (a4*r*r*r*r + a3*r*r*r + a2*r*r + a1*r + a0);
    fprime = (4.0*a4*r*r*r + 3.0*a3*r*r + 2.0*a2*r + a1);
    ru = r - f / fprime;
    n++;	    
  }
  
//  if (n==15) 
//    printf("[Error in Newton-Rhapson convergence, too many iterations]\n");

  return ru;
}		

///////////////////////////////////////////////////////////////////////////////
// Monte-Carlo integration
///////////////////////////////////////////////////////////////////////////////

void sobseq2D(double x[], const bool& init_call)
{
  int j,k,l;
  unsigned int m;
  unsigned long i, im, ipp;
  static double fac;
  static unsigned long in, ix[3], *iu[MAXBIT+1];
  static unsigned long mdeg[3]={0,1,2};
  static unsigned long ip[3]={0,0,1};
  static unsigned long iv[2*MAXBIT+1]=
    {0,1,1,1,1,1,1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9};
    
  if(init_call)
  {
    for (k=1;k<=2;k++) ix[k]=0;
    in=0;
    if (iv[1] != 1) return;
    fac=1.0/(1L << MAXBIT);
    for (j=1, k=0; j<=MAXBIT; j++,k+=2) iu[j] = &iv[k];
    for (k=1;k<=2;k++)
    {
      for(m=1;m<=mdeg[k];m++) iu[m][k] <<= (MAXBIT-m);
      for (j=mdeg[k]+1;j<=MAXBIT;j++)
      {
        ipp=ip[k];
        i=iu[j-mdeg[k]][k];
        i ^= (i >> mdeg[k]);
        for (l=mdeg[k]-1;l>=1;l--)
        {
          if (ipp & 1) i ^= iu[j-1][k];
          ipp >>= 1;
        }
        iu[j][k]=i;
      }
    }
  } else
  {
    im=in++;
    for(j=1;j<=MAXBIT;j++)
    {
      if (!(im & 1)) break;
      im >>= 1;
    }
    if (j > MAXBIT) printf("MAXBIT to small in sobseq2D\n");
    im = (j-1)*2;
    for(k=1;k<=2;k++)
    {
      ix[k] ^= iv[im+k];
      x[k]=ix[k]*fac;
    }
  }
}

double sobseq_integrate(double (*f)(double, double), double x_b, double x_e,
  double y_b, double y_e, double acc, int min_points)
{
  double switch_factor = 1; // correct for switching boundaries
  double vol;
  double x[3]; // when we have rewritten sobseq2D to properly start at 0,
    // we need of course to change all this as well!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  unsigned long n = 0; // number of points
  double result = 0.0;
  double var = 0.0;
  double error;
  double dist_x;
  double dist_y;
  double f_curr;
  double result_over_n;
  bool sobseqswitch;

  // swap boundaries if necessary and calculate volume and distances
  if (y_e < y_b) { swap(y_b, y_e); switch_factor = switch_factor * (-1.0); }
  if (x_e < x_b) { swap(x_b, x_e); switch_factor = switch_factor * (-1.0); }
  dist_x = x_e - x_b;
  dist_y = y_e - y_b;
  vol = dist_x * dist_y;

  // initialize sobseq sequence
  sobseqswitch = true;
  sobseq2D(x, sobseqswitch);
  
  // first the minimum number of points
  for(n=0;n<(unsigned int) min_points;n++)
  {
    sobseqswitch = false;
    sobseq2D(x, sobseqswitch);
    f_curr = f( x[1] * dist_x, x[2]* dist_y );
    result += f_curr;
    var += f_curr*f_curr; // needed for error estimation later on
  }
  
  // then onward to the desired accuracy
  do
  {
    n++;
    sobseqswitch = false;
    sobseq2D(x, sobseqswitch);
    f_curr = f( x[1] * dist_x, x[2]* dist_y );
    result += f_curr;
    var += SQR(f_curr); // needed for error estimation
    result_over_n = result / (double) n;
    error = sqrt((var / (double) n - SQR(result_over_n)) / (double) n) 
      / (result_over_n);
  } while( error > acc );

  return switch_factor * vol * result / (double) n;
}

////////////////////////////////////////////////////////////////////////////////
// integral of 2nd degree polynomial
////////////////////////////////////////////////////////////////////////////////

double integrate_polynome(double x1, double x2, double x3, double y1, double y2,
  double y3 )
{
  double Dx31 = x3 - x1;
  double invDx21 = 1.0 / ( x2 - x1);
  double invDx31 = 1.0 / Dx31;
  double invDx32 = 1.0 / ( x3 - x2);
//  printf("invDx21 = %e, invDx31 = %e, invDx32 = %e\n", invDx21, invDx31, invDx32);
  return (y1 * invDx21 * invDx31 - y2 * invDx21 * invDx32 + y3 * invDx31 *
    invDx32 ) * (x3*x3*x3-x1*x1*x1 ) * onethird + ( - y1 * (x3+x2) * invDx21 *
    invDx31 + y2 * (x1+x3) * invDx21 * invDx32 - y3 * (x1+x2) * invDx31 *
    invDx32 ) * 0.5 * (x3*x3 - x1 * x1) + ( x3*x2*y1 * invDx21 * invDx31 - 
    y2*x1*x3 * invDx21 * invDx32 + y3*x1*x2 * invDx31*invDx32 ) * Dx31;
}

////////////////////////////////////////////////////////////////////////////////
// polynomial interpolation
////////////////////////////////////////////////////////////////////////////////

double pol_int( double xa[], double ya[], int n, double x, double *dy )
{
  int i, m, ns = 0;
  double den, dif, dift, ho, hp, w;
  double c[n]; double d[n];
  double y;
  
  dif = fabs(x - xa[0]);
  
  // find xa entry closest to x and initialise c, d
  for (i = 0; i < n; i++)
  {
    dift = fabs(x - xa[i]);
    if (dift < dif) 
    { ns = i; dif = dift; }
    c[i] = ya[i];
    d[i] = ya[i];
  }
  
  // initial guess for y
  y = ya[ns];
  ns--;
  
  for(m = 1; m < n; m++) // loop over nobase columns in tableau
  {
    for (i = 0; i < n - m; i++) // loop over rows in tableau
    {
      ho = xa[i] - x;
      hp = xa[i + m] - x;
      w = c[i + 1] - d[i];
      den = w / (ho - hp);
      d[i] = hp * den;
      c[i] = ho * den;
    }
    y += (*dy = (2*(ns+1) < (n-m) ? c[ns+1] : d[ns--] ));
  }

  return y;
  
}

double pol_int( double xa[], double ya[], int n, double *dy )
{
  // NB the whole point was to write a smarter routine for x=0.0.....so do that at some point!
  return pol_int( xa, ya, n, 0.0, dy );
}


void pol_int(double xa[], double ya[], int n, double x, double *y, double *dy)
{
  *y = pol_int(xa, ya, n, x, dy);
}

///////////////////////////////////////////////////////////////////////////////
// Runge-Kutta algorithms
///////////////////////////////////////////////////////////////////////////////

double RungeKutta_fixed_step(double (*dydx)(double,double), double y_b, double 
  x_b, double x_e, long n_init)
{
	double h, half_h; // stepsize
	long int n; // number of steps
	double x_n, y_n; // intermediate values for x and y
	double k_1, k_2, k_3, k_4; // midpoints
	
	y_n = y_b; x_n = x_b; // start at provided boundary
	
	h = (x_e - x_b) / (double) n_init; // calculate stepsize
	half_h = 0.5*h;
	
	for (n=0; n<n_init; n++)
	{
		k_1 = h * dydx(x_n, y_n);
		k_2 = h * dydx(x_n + half_h, y_n + 0.5*k_1);
		k_3 = h * dydx(x_n + half_h, y_n + 0.5*k_2);
		k_4 = h * dydx(x_n + h, y_n + k_3);
		y_n += onesixth * (k_1 + k_4) + onethird * (k_2 + k_3);
		x_n += h;
	}
	
	return y_n; // sign h to correct for case x_b > x_e
}

double RungeKutta_adaptive_step(double (*dydx)(double), double y_b, 
  double x_b, double x_e, long n_init, int levels, double error)
{
	double h, h_half, h_quarter, h_threequarter;
	double x_n, y_n; // intermediate position and result
	int current_level = 0; // counts how often stepsize decreased
	double k1, k2, k3, k4; // intermediate evaluations big step
	double l1, l2, l3, l4; // intermediate evaluations 1st small step
	double m1, m2, m3, m4; // intermediate evaluations 2nd small step
	double y_big; // new y according to big step
	double y_small; // new y according to first and then both small steps
	double current_error; // current error
	bool error_found; // true if unacceptable error and step has to be redone
	double error_per_h; // allowed error per step
	long i; // counts how many steps still to go to b_e
	
	if (n_init < 2) n_init = 2;
	h = (x_e - x_b) / (double) (n_init / 2); // h is stepsize for the big steps
	h_half = h / 2.0;
	h_quarter = h_half / 2.0;
	h_threequarter = h_half + h_quarter;
	error_per_h = error / (double) (n_init / 2);
	i = n_init / 2;
	
	x_n = x_b; y_n = y_b;
	
	do
	{
		// one big step
		k1 = h * dydx(x_n);
		k2 = h * dydx(x_n + h_half);
		k3 = h * dydx(x_n + h_half);
		k4 = h * dydx(x_n + h);
		y_big = y_n + onesixth * (k1 + k4) + onethird * (k2 + k3);
		
		// 1st small step
		l1 = 0.5 * k1;
		l2 = h_half * dydx(x_n + h_quarter);
		l3 = h_half * dydx(x_n + h_quarter);
		l4 = h_half * dydx(x_n + h_half);
		y_small = y_n + onesixth * (l1 + l4) + onethird * (l2 + l3);
		
		// 2nd small step
    m1 = h_half * dydx(x_n + h_half);
    m2 = h_half * dydx(x_n + h_threequarter);
    m3 = h_half * dydx(x_n + h_threequarter);
    m4 = h_half * dydx(x_n + h);
    y_small += onesixth * (m1 + m4) + onethird * (m2 + m3);
		
		// compare errors
		current_error = y_small - y_big;
		if (fabs(current_error) > fabs(error_per_h))
		{
			if (current_level < levels)
			{
				error_found = true;
				h = 0.5*h;
				h_half = 0.5*h_half;
				h_quarter = 0.5*h_quarter;
				h_threequarter = 0.5*h_threequarter;
				i = 2*i;
				current_level++;
				error_per_h = 0.5*error_per_h;
			} else
			{
				error_found = false;
				x_n += h; i--;
			}
		} else 
		{
			error_found = false;
			x_n += h; i--;
			if(current_level > 0 && even(i))
			{
				h = 2.0*h;
				i = i / 2;
				h_half = 2.0*h_half;
				h_quarter = 2.0*h_quarter;
				h_threequarter = 2.0*h_threequarter;
				current_level--;
				error_per_h = 2.0*error_per_h;
			}
		}
		
		// store intermediate result
		if (!error_found)
		{
		  y_n = y_small + onefifteenth * current_error; // local extrapolation
		} 
	} while(i>0);
	
	return y_n;
}

double RungeKutta_adaptive_step(double (*dydx)(double,double), double y_b, 
  double x_b, double x_e, long n_init, int levels, double error)
{
	double h, h_half, h_quarter, h_threequarter;
	double x_n, y_n; // intermediate position and result
	int current_level = 0; // counts how often stepsize decreased
	double k1, k2, k3, k4; // intermediate evaluations big step
	double l1, l2, l3, l4; // intermediate evaluations 1st small step
	double m1, m2, m3, m4; // intermediate evaluations 2nd small step
	double y_big; // new y according to big step
	double y_small; // new y according to first and then both small steps
	double current_error; // current error
	bool error_found; // true if unacceptable error and step has to be redone
	double error_per_h; // allowed error per step
	long i; // counts how many steps still to go to b_e
	
	if (n_init < 2) n_init = 2;
	h = (x_e - x_b) / (double) (n_init / 2); // h is stepsize for the big steps
	h_half = h / 2.0;
	h_quarter = h_half / 2.0;
	h_threequarter = h_half + h_quarter;
	error_per_h = error / (double) (n_init / 2);
	i = n_init / 2;
	
	x_n = x_b; y_n = y_b;
	
	do
	{
		// one big step
		k1 = h * dydx(x_n, y_n);
		k2 = h * dydx(x_n + h_half, y_n + 0.5*k1);
		k3 = h * dydx(x_n + h_half, y_n + 0.5*k2);
		k4 = h * dydx(x_n + h, y_n + k3);
		y_big = y_n + onesixth * (k1 + k4) + onethird * (k2 + k3);
		
		// 1st small step
		l1 = 0.5 * k1;
		l2 = h_half * dydx(x_n + h_quarter, y_n + 0.5*l1);
		l3 = h_half * dydx(x_n + h_quarter, y_n + 0.5*l2);
		l4 = h_half * dydx(x_n + h_half, y_n + l3);
		y_small = y_n + onesixth * (l1 + l4) + onethird * (l2 + l3);
		
		// 2nd small step
    m1 = h_half * dydx(x_n + h_half, y_small);
    m2 = h_half * dydx(x_n + h_threequarter, y_small + 0.5*m1);
    m3 = h_half * dydx(x_n + h_threequarter, y_small + 0.5*m2);
    m4 = h_half * dydx(x_n + h, y_small + m3);
    y_small += onesixth * (m1 + m4) + onethird * (m2 + m3);
		
		// compare errors
		current_error = y_small - y_big;
		if (fabs(current_error) > fabs(error_per_h))
		{
			if (current_level < levels)
			{
				error_found = true;
				h = 0.5*h;
				h_half = 0.5*h_half;
				h_quarter = 0.5*h_quarter;
				h_threequarter = 0.5*h_threequarter;
				i = 2*i;
				current_level++;
				error_per_h = 0.5*error_per_h;
			} else
			{
				error_found = false;
				x_n += h; i--;
			}
		} else 
		{
			error_found = false;
			x_n += h; i--;
			if(current_level > 0 && even(i))
			{
				h = 2.0*h;
				i = i / 2;
				h_half = 2.0*h_half;
				h_quarter = 2.0*h_quarter;
				h_threequarter = 2.0*h_threequarter;
				current_level--;
				error_per_h = 2.0*error_per_h;
			}
		}
		
		// store intermediate result
		if (!error_found)
		{
		  y_n = y_small + onefifteenth * current_error; // local extrapolation
		} 
	} while(i>0);
	
	return y_n;
}

double RungeKutta_adaptive_step(double (*dydx)(double, double, double),
  double y_b, double x_b, double x_e, double z, long n_init, int levels, 
  double error)
{
	double h, h_half, h_quarter, h_threequarter;
	double x_n, y_n; // intermediate position and result
	int current_level = 0; // counts how often stepsize decreased
	double k1, k2, k3, k4; // intermediate evaluations big step
	double l1, l2, l3, l4; // intermediate evaluations 1st small step
	double m1, m2, m3, m4; // intermediate evaluations 2nd small step
	double y_big; // new y according to big step
	double y_small; // new y according to first and then both small steps
	double current_error; // current error
	bool error_found; // true if unacceptable error and step has to be redone
	double error_per_h; // allowed error per step
	long i; // counts how many steps still to go to b_e
	
	if (n_init < 2) n_init = 2;
	h = (x_e - x_b) / (double) (n_init / 2); // h is stepsize for the big steps
	h_half = h / 2.0;
	h_quarter = h_half / 2.0;
	h_threequarter = h_half + h_quarter;
	error_per_h = error / (double) (n_init / 2);
	i = n_init / 2;
	
	x_n = x_b; y_n = y_b;
	
	do
	{
		// one big step
		k1 = h * dydx(x_n, y_n,z);
		k2 = h * dydx(x_n + h_half, y_n + 0.5*k1,z);
		k3 = h * dydx(x_n + h_half, y_n + 0.5*k2,z);
		k4 = h * dydx(x_n + h, y_n + k3,z);
		y_big = y_n + onesixth * (k1 + k4) + onethird * (k2 + k3);
		
		// 1st small step
		l1 = 0.5 * k1;
		l2 = h_half * dydx(x_n + h_quarter, y_n + 0.5*l1,z);
		l3 = h_half * dydx(x_n + h_quarter, y_n + 0.5*l2,z);
		l4 = h_half * dydx(x_n + h_half, y_n + l3,z);
		y_small = y_n + onesixth * (l1 + l4) + onethird * (l2 + l3);
		
		// 2nd small step
    m1 = h_half * dydx(x_n + h_half, y_small,z);
    m2 = h_half * dydx(x_n + h_threequarter, y_small + 0.5*m1,z);
    m3 = h_half * dydx(x_n + h_threequarter, y_small + 0.5*m2,z);
    m4 = h_half * dydx(x_n + h, y_small + m3,z);
    y_small += onesixth * (m1 + m4) + onethird * (m2 + m3);
		
		// compare errors
		current_error = y_small - y_big;
		if (fabs(current_error) > fabs(error_per_h))
		{
			if (current_level < levels)
			{
				error_found = true;
				h = 0.5*h;
				h_half = 0.5*h_half;
				h_quarter = 0.5*h_quarter;
				h_threequarter = 0.5*h_threequarter;
				i = 2*i;
				current_level++;
				error_per_h = 0.5*error_per_h;
			} else
			{
				error_found = false;
				x_n += h; i--;
			}
		} else 
		{
			error_found = false;
			x_n += h; i--;
			if(current_level > 0 && even(i))
			{
				h = 2.0*h;
				i = i / 2;
				h_half = 2.0*h_half;
				h_quarter = 2.0*h_quarter;
				h_threequarter = 2.0*h_threequarter;
				current_level--;
				error_per_h = 2.0*error_per_h;
			}
		}
		
		// store intermediate result
		if (!error_found)
		{
		  y_n = y_small + onefifteenth * current_error; // local extrapolation
		} 
	} while(i>0);
	
	return y_n;
}

////////////////////////////////////////////////////////////////////////////////
// Special Functions
////////////////////////////////////////////////////////////////////////////////

double gammafunc(double x)
{
  static double result, series;
  static int j;
  static const double cof[6] = {76.18009172947126, -86.50532032941677,
    24.01409824083091, -1.231739572450155, 0.1208650973866179e-2,
    -0.5395239384953e-5};
  result = pow( x + 4.5, x - 0.5 );
  result *= exp( -(x + 4.5 ) );
  series = 1.000000000190015;
  for (j=0; j<6; j++) series += cof[j]/x++;
  return result * 2.5066282746310005 * series;
}

////////////////////////////////////////////////////////////////////////////////

// auxiliray function chebev
double chebev(double a, double b, double c[], int m, double x)
{
  double d=0.0, dd=0.0, sv, y, y2;
  int j;
  
  if ((x-a)*(x-b) > 0.0) { printf("x not in range chebev\n"); abort(); }
  y2 = 2.0 * (y=(2.0*x-a-b)/(b-a));
  for (j=m-1; j>=1; j-- )
  {
    sv=d;
    d=y2*d-dd+c[j];
    dd=sv;
  }
  return y*d-dd+0.5*c[0];
}

////////////////////////////////////////////////////////////////////////////////

// auxiliary function 'beschb'
void beschb(double x, double *gam1, double *gam2, double *gamp1, double *gammi)
{
  static const int NUSE1 = 5;
  static const int NUSE2 = 5;
  
  double xx;
  static double c1[] = { -1.142022680371168e0, 6.5165112670737e-3,
    3.087090173086e-4, -3.4706269649e-6, 6.9437664e-9, 3.67795e-11, -1.356e-13};
  static double c2[] = { 1.843740587300905e0, -7.68528408447867e-2,
    1.2719271366546e-3, -4.9717367042e-6, -3.31261198e-8, 2.423096e-10,
    -1.702e-13, -1.49e-15 };
    
  xx=8.0*x*x-1.0;
  *gam1 = chebev(-1.0, 1.0, c1, NUSE1, xx );
  *gam2 = chebev(-1.0, 1.0, c2, NUSE2, xx);
  *gamp1 = *gam2-x*(*gam1);
  *gammi = *gam2+x*(*gam1);
}

////////////////////////////////////////////////////////////////////////////////

void bessik(double x, double xnu, double *ri, double *rk, double *rip, double
  *rkp)
{
  static const double EPS = 1.0e-16;
  static const double FPMIN = 1.0e-30;
  static const int MAXIT = 10000;
  static const double XMIN = 2.0;
  
  int i, l, nl;  
  double a, a1, b, c, d, del, del1, delh, dels, e, f, fact, fact2, ff, gam1,
    gam2, gammi, gampl, h, p, pimu, q, q1, q2, qnew, ril, ril1, rimu, rip1,
    ripl, ritemp, rk1, rkmu, rkmup, rktemp, s, sum, sum1, x2, xi, xi2, xmu,
    xmu2;
  
  if (x <= 0.0 || xnu < 0.0 ) { printf("bad arg in bessik\n"); abort(); };
  nl = (int) (xnu+0.5);
  xmu=xnu-nl;
  xmu2 = xmu * xmu;
  xi = 1.0 / x;
  xi2 = 2.0 * xi;
  h = xnu * xi;
  if ( h < FPMIN) h = FPMIN;
  b = xi2 * xnu;
  d= 0.0;
  c = h;
  for ( i=1; i<=MAXIT; i++)
  {
    b+= xi2;
    d=1.0/(b+d);
    c=b+1.0/c;
    del=c*d;
    h=del*h;
    if (fabs(del-1.0) < EPS) break;
  }
  
  if (i > MAXIT) { printf("x too large in bessik; try asymptotic expansion.\n");
    abort(); }
  ril=FPMIN;
  ripl=h*ril;
  ril1=ril;
  rip1=ripl;
  fact=xnu*xi;
  for (l=nl; l>=1; l-- )
  {
    ritemp = fact*ril + ripl;
    fact -= xi;
    ripl = fact*ritemp+ril;
    ril = ritemp;
  }  
  f = ripl / ril;
  if (x < XMIN)
  {
    x2 = 0.5*x;
    pimu = PI*xmu;
    fact = (fabs(pimu) < EPS ? 1.0 : pimu/sin(pimu));
    d = -log(x2);
    e=xmu*d;
    fact2 = (fabs(e) < EPS ? 1.0 : sinh(e)/e);
    beschb( xmu, &gam1, &gam2, &gampl, &gammi );
    ff=fact*(gam1*cosh(e)+gam2*fact2*d);
    sum = ff;
    e=exp(e);
    p=0.5*e/gampl;
    q=0.5/(e*gammi);
    c=1.0;
    d=x2*x2;
    sum1=p;
    for ( i=1; i<= MAXIT; i++ ) {
      ff = (i*ff+p+q)/(i*i - xmu2);
      c *= (d/i);
      p /= (i-xmu);
      q /= (i+xmu);
      del = c*ff;
      sum += del;
      del1 = c*(p-i*ff);
      sum1 += del1;
      if (fabs(del) < fabs(sum)*EPS ) break;
    }
    if (i > MAXIT) { printf("bessk series failed to converge\n"); abort(); }
    rkmu=sum;
    rk1=sum1*xi2;
  } else
  {
    b = 2.0*(1.0+x);
    d=1.0/b;
    h=delh=d;
    q1=0.0;
    q2=1.0;
    a1=0.25-xmu2;
    q=c=a1;
    a = -a1;
    s=1.0+q*delh;
    for (i=2; i<=MAXIT; i++)
    {
      a -= 2*(i-1);
      c = -a*c/i;
      qnew = (q1-b*q2) / a;
      q1 = q2;
      q2 = qnew;
      q += c*qnew;
      b += 2.0;
      d = 1.0/(b+a*d);
      delh=(b*d-1.0)*delh;
      h += delh;
      dels = q*delh;
      s += dels;
      if ( fabs(dels/s) < EPS ) break;
    }
    if (i > MAXIT) { printf("bessik failed to converge in cf2\n"); abort(); }
    h = a1*h;
    rkmu = sqrt(PI/(2.0*x))*exp(-x)/s;
    rk1=rkmu*(xmu+x+0.5-h)*xi;
  }
  rkmup=xmu*xi*rkmu-rk1;
  rimu=xi/(f*rkmu-rkmup);
  *ri=(rimu*ril1)/ril;
  *rip=(rimu*rip1)/ril;
  for (i=1; i<=nl; i++)
  {
    rktemp=(xmu+1)*xi2*rk1+rkmu;
    rkmu=rk1;
    rk1=rktemp;
  }
  *rk=rkmu;
  *rkp=xnu*xi*rkmu-rk1;
}

// END CODE NOT PUBLIC
