template <class T> double modmid(T *a_base, int integrand_ID, double x_0, 
  double x_1, double dydx0, double dydx1, int n_steps )
// dydx is included as argument because we expect modmid routines to be
// performed in adjacent domains. We don't want to do double work when the
// new domain's left boundary is equal to the previous domain's right boundary
// and an earlier modmid routine already gave us the integrand at the boundary.
{
  int n;
  double ym, yn, x, swap, h2, h, yout;
  
  // determine stepsize
  h = (x_1 - x_0) / n_steps;
  
  // first step
  ym = 0.0; yn = h * dydx0;
  
  x = x_0 + h;
  yout = (*a_base).relay_integrand_1D( integrand_ID, x );
    // yout functions as temporary storage of derivative
  
  // general step
  h2 = 2.0 * h;
  for ( n=2; n < n_steps; n++ )
  {
    swap = ym + h2 * yout;
    ym = yn;
    yn = swap;
    x += h;
    yout = (*a_base).relay_integrand_1D( integrand_ID, x );
  }
  
  // now take advantage of known dydx1
  swap = ym + h2 * yout;
  ym = yn;
  yn = swap;
  x+= h;
  yout = dydx1;
  
  // last step and return
  return 0.5 * (ym + yn + h * yout );
}

////////////////////////////////////////////////////////////////////////////////

template <class T> double integrate_1D_step( T *a_base, int integrand_ID, double
  x, double *stepsize, double steperr, double *abserr, double dydx0, 
  double dydx1, double h_min )
{
  int m; // no of steps is 2*(m+1)
  double ReP; // result of Richardson extrapolation and error
  double errin = abs(steperr); // errin stores requested maximum error
  double y[10]; // contains result for modmid with n = 2, 4, 6, ...
  double h[10];

  // get started with n=2. We need at least this much accuracy
  y[0] = modmid( a_base, integrand_ID, x, x+(*stepsize), dydx0, dydx1, 2 );
  h[0] = 0.5; // for ReP we go to limit h -> 0
  m = 1;
//  printf("fractional error limit: %e\n", steperr );
  
  // try at this stepsize
  while( m < 6 )
  {
    y[m] = modmid( a_base, integrand_ID, x, x+(*stepsize),dydx0, dydx1,2*(m+1));
    h[m] = 1.0 / ( 2*(m+1) );
    ReP = pol_int( h, y, m+1, abserr );
    *abserr = abs(*abserr);
    if( ReP == 0.0 )
    {
      if (y[m] == 0.0) { return ReP; }
      else steperr = abs((*abserr) / y[m] );
    } else
    steperr = abs((*abserr) / ReP); // change from abs to frac err
    if ( steperr < errin ) return ReP; 
    if ( abs( y[m] - y[m-1] ) < errin * y[m] ) 
    {
      *abserr = abs( y[m] - y[m-1] );
      return y[m];
    }
    m++;
  }

  // check if we still can decrease stepsize
  if (*stepsize <= h_min) 
  {
  // TURN BACK ON EVENTUALLY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//    printf("when integration over ID = %d\n", integrand_ID);
//    printf("steperr = %e\n", steperr );
//    printf("minimum stepsize (%e) reached!\n", h_min); fflush(stdout);
    
    return y[m-1];
  }
  
  // we can and have to decrease stepsize
  *stepsize = (*stepsize) / 2.0;
  dydx1 = (*a_base).relay_integrand_1D( integrand_ID, x+*stepsize );
  return integrate_1D_step( a_base, integrand_ID, x, stepsize, errin, abserr, 
    dydx0, dydx1, h_min );
}

////////////////////////////////////////////////////////////////////////////////

template <class T> double integrate_1D( T *a_base, int integrand_ID, double x_0, 
  double x_1, double res, double *err, double h_min )
{
//  printf("starting integration over ID %d\n", integrand_ID );
//  printf("domain is %e - %e, res = %e, err %e, h_min=%e\n", x_0, x_1, res, *err,
//    h_min );
    
  double y;
  double x;
  double H, H_max; // stepsize, maximum stepsize
  double steperr, sumabserr; // error, weighed sum of squared errors
  double abserr;
  double domain, inv_domain; // integration domain
  double dydx0, dydx1; // integrand values at the bounds of subdomain H
  
  // check domain size nonzero
  if (x_1 == x_0)
  {
    *err = 0.0;
    return 0.0;
  }
  
  // make sure domain boundaries are in correct order
  if (x_1 < x_0) 
    return -integrate_1D( a_base, integrand_ID, x_1, x_0, res, err, h_min);
  
  // if res and / or h_min provided equal to zero, it means that the routine
  // has to make its own guess
  if ( res == 0.0 ) res = *err * ( x_1 - x_0 );
  if ( h_min == 0.0 ) h_min = 1e-2 * res;
  
  // if h_min of same order as domain, then this integral will not contribute
  // much of significance. We give just a 3 points modmid result
  if ( h_min >= (x_1 - x_0) )
  {
    dydx0 = (*a_base).relay_integrand_1D( integrand_ID, x_0 );
    dydx1 = (*a_base).relay_integrand_1D( integrand_ID, x_1 );
    double y2, y3; // two calls to modmid to have an error estimate as well
    y2 = modmid( a_base, integrand_ID, x_0, x_1, dydx0, dydx1, 2 );
    y3 = modmid( a_base, integrand_ID, x_0, x_1, dydx0, dydx1, 3 );
    *err = abs(y2 - y3);
    printf("integral returns (small domain): %e\n", y3);
    return y3;
  }
  
  // declare variables
  H_max = H = res; // set stepsizes
  x = x_0;
  sumabserr = 0; // cumulative error
  domain = x_1 - x_0;
  inv_domain = 1.0 / domain;
  y = 0.0;
  dydx0 = (*a_base).relay_integrand_1D( integrand_ID, x );
  
  while (x < x_1)
  {
    if (x + H > x_1) H = x_1 - x; // correct for overshooting

    steperr = (*err); //reset frac error (it will be modified to obtained error)
    dydx1 = (*a_base).relay_integrand_1D( integrand_ID, x+H );
    y += integrate_1D_step( a_base, integrand_ID, x, &H, steperr, &abserr, dydx0, 
      dydx1, h_min ); // steperr is now fractional error
    sumabserr += abserr;
    x += H;
    dydx0 = dydx1; // use right boundary value of the integrand of this step as 
      // left boundary value of the integrand for the next step
    H *= 2.0; // last attempt was succesful, so try a bigger step
    if (H > H_max) H = H_max; // but not too big
  }
  
  *err = sumabserr;
  
//  printf("integral returns: %e\n",y);
  return y;
}

////////////////////////////////////////////////////////////////////////////////

template <class T> double integrate_1D( T *a_base, int integrand_ID, double x_0, 
  double x_1, int res_ID, double *err, double h_min )
{
//  printf("starting integration over ID %d\n", integrand_ID );
//  printf("domain is %e - %e, res = %e, err %e, h_min=%e\n", x_0, x_1, res, *err,
//    h_min );
    
  double y;
  double x;
  double H, H_max; // stepsize, maximum stepsize
  double steperr, sumabserr; // error, weighed sum of squared errors
  double abserr;
  double domain, inv_domain; // integration domain
  double dydx0, dydx1; // integrand values at the bounds of subdomain H
  
  // check domain size nonzero
  if (x_1 == x_0)
  {
    *err = 0.0;
    return 0.0;
  }
  
  // make sure domain boundaries are in correct order
  if (x_1 < x_0) 
    return -integrate_1D( a_base, integrand_ID, x_1, x_0, res_ID, err, h_min);
  
  // if res and / or h_min provided equal to zero, it means that the routine
  // has to make its own guess
  H_max = a_base->relay_integrand_1D( res_ID, x_0 ); // use smart stepsize max
  if ( h_min == 0.0 ) h_min = 1e-4 * H_max;
  
  // if h_min of same order as domain, then this integral will not contribute
  // much of significance. We give just a 3 points modmid result
  if ( h_min >= (x_1 - x_0) )
  {
    dydx0 = (*a_base).relay_integrand_1D( integrand_ID, x_0 );
    dydx1 = (*a_base).relay_integrand_1D( integrand_ID, x_1 );
    double y2, y3; // two calls to modmid to have an error estimate as well
    y2 = modmid( a_base, integrand_ID, x_0, x_1, dydx0, dydx1, 2 );
    y3 = modmid( a_base, integrand_ID, x_0, x_1, dydx0, dydx1, 3 );
    *err = abs(y2 - y3);
    printf("integral returns (small domain): %e\n", y3);
    return y3;
  }
  
  // declare variables
  H_max = H = a_base->relay_integrand_1D( res_ID, x_0 ); // set stepsizes
  x = x_0;
  sumabserr = 0; // cumulative error
  domain = x_1 - x_0;
  inv_domain = 1.0 / domain;
  y = 0.0;
  dydx0 = (*a_base).relay_integrand_1D( integrand_ID, x );
  
  while (x < x_1)
  {
    if (x + H > x_1) H = x_1 - x; // correct for overshooting

    steperr = (*err); //reset frac error (it will be modified to obtained error)
    dydx1 = (*a_base).relay_integrand_1D( integrand_ID, x+H );
    y += integrate_1D_step( a_base, integrand_ID, x, &H, steperr, &abserr, dydx0, 
      dydx1, h_min ); // steperr is now fractional error
    sumabserr += abserr;
    x += H;
    dydx0 = dydx1; // use right boundary value of the integrand of this step as 
      // left boundary value of the integrand for the next step
    H *= 2.0; // last attempt was succesful, so try a bigger step
    H_max = a_base->relay_integrand_1D( res_ID, x ); // use smart stepsize max
    if (H > H_max) H = H_max; // but not too big
  }
  
  *err = sumabserr;
  
//  printf("integral returns: %e\n",y);
  return y;
}

////////////////////////////////////////////////////////////////////////////////

template <class T> double integrate_1D_step_abserr( T *a_base, int integrand_ID, 
  double x, double *stepsize, double steperr, double *abserr, double dydx0, 
  double dydx1, double h_min )
{
  int m; // no of steps is 2*(m+1)
  double ReP; // result of Richardson extrapolation and error
  double errin = abs(steperr); // errin stores requested maximum error
  double y[10]; // contains result for modmid with n = 2, 4, 6, ...
  double h[10];

  // get started with n=2. We need at least this much accuracy
  y[0] = modmid( a_base, integrand_ID, x, x+(*stepsize), dydx0, dydx1, 2 );
  h[0] = 0.5; // for ReP we go to limit h -> 0
  m = 1;
//  printf("fractional error limit: %e\n", steperr );
  
  // try at this stepsize
  while( m < 6 )
  {
    y[m] = modmid( a_base, integrand_ID, x, x+(*stepsize),dydx0, dydx1,2*(m+1));
    h[m] = 1.0 / ( 2*(m+1) );
    ReP = pol_int( h, y, m+1, abserr );
    *abserr = abs(*abserr);
    if( ReP == 0.0 )
    {
      if (y[m] == 0.0) { return ReP; }
      else steperr = abs((*abserr) / y[m] );
    } else
    if ( steperr < errin ) return ReP;
    if ( abs( y[m] - y[m-1] ) < errin * y[m] ) 
    {
      *abserr = abs( y[m] - y[m-1] );
      return y[m];
    }
    m++;
  }

  // check if we still can decrease stepsize
  if (*stepsize <= h_min) 
  {
  // TURN BACK ON EVENTUALLY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//    printf("when integration over ID = %d\n", integrand_ID);
//    printf("steperr = %e\n", steperr );
//    printf("minimum stepsize (%e) reached!\n", h_min); fflush(stdout);
//    abort();
    return y[m-1];
  }
  
  // we can and have to decrease stepsize
  *stepsize = (*stepsize) / 2.0;
  dydx1 = (*a_base).relay_integrand_1D( integrand_ID, x+*stepsize );
  return integrate_1D_step( a_base, integrand_ID, x, stepsize, errin, abserr, 
    dydx0, dydx1, h_min );
}

////////////////////////////////////////////////////////////////////////////////

template <class T> double integrate_1D_abserr( T *a_base, int integrand_ID, 
  double x_0, double x_1, double res, double *err, double h_min )
{
//  printf("starting integration over ID %d\n", integrand_ID );
//  printf("domain is %e - %e, res = %e, err %e, h_min=%e\n", x_0, x_1, res, *err,
//    h_min );
    
  double y;
  double x;
  double H, H_max; // stepsize, maximum stepsize
  double steperr, sumabserr; // error, weighed sum of squared errors
  double abserr;
  double domain, inv_domain; // integration domain
  double dydx0, dydx1; // integrand values at the bounds of subdomain H
  
  // check domain size nonzero
  if (x_1 == x_0)
  {
    *err = 0.0;
    return 0.0;
  }
  
  // make sure domain boundaries are in correct order
  if (x_1 < x_0) return 
    -integrate_1D_abserr( a_base, integrand_ID, x_1, x_0, res, err, h_min);
  
  // if res and / or h_min provided equal to zero, it means that the routine
  // has to make its own guess
  if ( res == 0.0 ) res = *err * ( x_1 - x_0 );
  if ( h_min == 0.0 ) h_min = 1e-2 * res;
  
  // if h_min of same order as domain, then this integral will not contribute
  // much of significance. We give just a 3 points modmid result
  if ( h_min >= (x_1 - x_0) )
  {
    dydx0 = (*a_base).relay_integrand_1D( integrand_ID, x_0 );
    dydx1 = (*a_base).relay_integrand_1D( integrand_ID, x_1 );
    double y2, y3; // two calls to modmid to have an error estimate as well
    y2 = modmid( a_base, integrand_ID, x_0, x_1, dydx0, dydx1, 2 );
    y3 = modmid( a_base, integrand_ID, x_0, x_1, dydx0, dydx1, 3 );
    *err = abs(y2 - y3);
    printf("integral returns (small domain): %e\n", y3);
    return y3;
  }
  
  // declare variables
  H_max = H = res; // set stepsizes
  x = x_0;
  sumabserr = 0; // cumulative error
  domain = x_1 - x_0;
  inv_domain = 1.0 / domain;
  y = 0.0;
  dydx0 = (*a_base).relay_integrand_1D( integrand_ID, x );
  
  while (x < x_1)
  {
    if (x + H > x_1) H = x_1 - x; // correct for overshooting

    steperr = (*err); //reset frac error (it will be modified to obtained error)
    dydx1 = (*a_base).relay_integrand_1D( integrand_ID, x+H );
    y += integrate_1D_step_abserr( a_base, integrand_ID, x, &H, steperr, 
      &abserr, dydx0, dydx1, h_min ); // steperr is now fractional error
    sumabserr += abserr;
    x += H;
    dydx0 = dydx1; // use right boundary value of the integrand of this step as 
      // left boundary value of the integrand for the next step
    H *= 2.0; // last attempt was succesful, so try a bigger step
    if (H > H_max) H = H_max; // but not too big
  }
  
  *err = sumabserr;
  
//  printf("integral returns: %e\n",y);
  return y;
}

////////////////////////////////////////////////////////////////////////////////
// 1D integration routines with implicit discrete domain
////////////////////////////////////////////////////////////////////////////////

template < class I, class T > class integral_1D_dd
{
  public:
    // input variables, mathematical
    double frac_err; // fractional error
    long S0, S1; // domain

    // input variables, pointers
    I *I_base; // pointer to class containing the integrand
    int integrand_ID; // number given to relay function within I_base to
      // identify the integrand function
    T *T_base; // pointer to class containing transformation function from S to
      // real (continuous) domain
    int transform_ID; // number given to relay function within T_base to
      // identify the transformation function

    // input variables, procedural settings for integration algorithm
    long dS_max; // initial resolution
    int level_max;
   
    // output variables
    double abs_err; // absolute error
    double y; // the integral
    void solve(); // calculates y
    long no_set; // stores number of integrand evaluation calls;
    double* dydx;
    bool* dydx_set;
      
  protected:
    long S, dS; // current position, stepsize
    void improve_domain(); // replaces dS with a dS smaller or equal than it
      // that has advantageous divisor properties with respect to the integers
      // up to and including level_max
    void set_dydx( long i );
    double integrate( int level ); // integration via Riemann sum
    double modmid( int level ); // integration via modified midpoint
    double dxds( long s );
    long f( long level);
};

// note on template notation: first line says we're dealing with a template
// class member function. Seconds <I,T> says that the template variables I and
// T are actually declared within the class
template < class I, class T >
void integral_1D_dd<I,T> :: set_dydx( long i )
{
  if ( !dydx_set[i - S0] )
  {
    dydx_set[i - S0] = true;
    no_set++;
    dydx[ i - S0 ] = (*I_base).relay_integrand_1D_dd( integrand_ID, i );
  }
}

template < class I, class T>
double integral_1D_dd<I,T> :: dxds( long s )
{
  if ( s == S0) return (*T_base).relay_integrand_1D_dd( transform_ID, S0+1) -
   (*T_base).relay_integrand_1D_dd( transform_ID, S0);
  if ( s == S1 ) return (*T_base).relay_integrand_1D_dd( transform_ID, S1) -
   (*T_base).relay_integrand_1D_dd( transform_ID, S1 - 1);

  // we're not at the edges
  return 0.5 * ( (*T_base).relay_integrand_1D_dd( transform_ID, s+1) -
   (*T_base).relay_integrand_1D_dd( transform_ID, s-1) );
}

template < class I, class T >
void integral_1D_dd<I,T> :: improve_domain()
{
  int cur_max_level; // level max may be out of reach due to small domain size

  // first, ensure that cur_max_level requires smaller (or equal) domain than dS
  cur_max_level = level_max;
  while ( f( cur_max_level) > dS )
  { cur_max_level--; }
  
  // now refine dS to number divisible by all levels up to cur_max_level
  while ( dS % f( cur_max_level) != 0 ) { dS--; }
}

template < class I, class T > 
long integral_1D_dd<I,T> :: f( long level)
{
  switch( level )
  {
    case 1: return 1; break;
    case 2: return 2; break;
    case 3: return 6; break;
    case 4: return 12; break;
    case 5: return 60; break;
    case 6: return 120; break;
    case 7: return 840; break;
    case 8: return 1680; break;
    case 9: return 1890; break;
    case 10: default: return 3780; break;
  }
}

template < class I, class T >
double integral_1D_dd<I,T> :: integrate( int level )
{ // this routine is called only if dS % level == 0
  int i;
  double dx;
  double dy;
  long h;

  // calculate integral, using polynome with subdepth points
  // for now a simple approx, just the mean
  dy = 0;
  h = dS / level;

  dx = 0.5 * ( (*T_base).relay_integrand_1D_dd( transform_ID, S + h) -
       (*T_base).relay_integrand_1D_dd( transform_ID, S ) );
  set_dydx( S );
  dy = dx * dydx[ S - S0 ];
//  printf("  ( dy = %e, dx = %e )\n", dy, dx );
  for( i=1; i< level; i++ )
  {
    dx = 0.5 * ( (*T_base).relay_integrand_1D_dd( transform_ID, S + (i+1)*h) -
       (*T_base).relay_integrand_1D_dd( transform_ID, S + (i-1)*h ) );
    set_dydx( S + i*h );
    dy += dx * dydx[ S + i*h - S0 ];
//    printf("  ( dy = %e, dx = %e )\n", dy, dx );
  }
  dx = 0.5 * ( (*T_base).relay_integrand_1D_dd( transform_ID, S + i*h ) - 
      (*T_base).relay_integrand_1D_dd( transform_ID, S + (i-1) * h ) );
  set_dydx( S + i*h);
  dy += dx * dydx[ S + i*h - S0 ];
//  printf("  ( dy = %e, dx = %e )\n", dy, dx );
  
  printf("domain: [%ld .. %ld] [%e .. %e], result: %e, level %d\n", S, S + dS, 
    (*T_base).relay_integrand_1D_dd( transform_ID, S ),
    (*T_base).relay_integrand_1D_dd( transform_ID, S + dS ),
    dy, level);  
  return dy; 
}

template <class I, class T> 
double integral_1D_dd<I,T> :: modmid( int level )
{
  int n; // dummy index counting number of steps
  long S_sub, h, h2;
  double ym, yn, swap, yout;
  
  // determine stepsize
  h = dS / level;
  
  // first step
  ym = 0.0; 
  set_dydx( S );
  yn = h * dxds( S ) * dydx[ S - S0 ];
  S_sub = S + h;  
  set_dydx( S_sub );

  yout = dxds( S_sub) * dydx[ S_sub - S0 ];
    // yout functions as temporary storage of derivative
  
  // general step
  h2 = 2 * h;
  for ( n=2; n <= level; n++ )
  {
    swap = ym + h2 * yout;
    ym = yn;
    yn = swap;
    S_sub += h;
    set_dydx( S_sub);
    yout = dydx[ S_sub - S0 ] * dxds( S_sub );
  }
  
  // last step and return
  return 0.5 * (ym + yn + h * yout );
}

template < class I, class T > 
void integral_1D_dd<I,T> :: solve()
{
  long i; // dummy index variable
  int level;
  double dy0, dy1;
  double dx;
  double suberr;
  
  // sanity check on dS_max and order of boundaries
  if (S1 < S0)
  {
    i = S1; S1 = S0; S0 = i; // swap boundaries
    solve(); // call this function (which calculates y)
    y *= -1.0; // correct for interchanging order boundaries
  }
  if (dS_max > S1 - S0) dS_max = S1 - S0;
  
  // clean variables and declare memory
  dS = dS_max;
  improve_domain();
  dS_max = dS; // this new dS_max is more sensible
  level = 2;
  y = 0.0;
  dydx = new double[S1 - S0 + 1];
  dydx_set = new bool[ S1 - S0 + 1];
  for( i = 0; i<= S1 - S0; i++ ) dydx_set[i] = false;
  abs_err = 0.0;
  no_set = 0;
  S = S0;

  dy0 = modmid( 1 );
  do
  {
    if ( dS == 1)
    {
      y += dy0;
      // special error handling
      dx = (*T_base).relay_integrand_1D_dd( transform_ID, S + 1 ) -
        (*T_base).relay_integrand_1D_dd( transform_ID, S );
      abs_err += abs( 0.5 * ( dydx[ S + 1 - S0 ] - dydx[ S - S0 ] ) * dx );
      // continue to next part of domain, attempting larger stepsize
      S++;
      dS = 2;
      if (S+dS > S1) dS = S1 - S;
      if ( S < S1 ) dy0 = modmid( 1 );
    }
    else
    {
      if ( dS % level == 0)
      {
        dy1 = modmid( level );
        suberr = abs( dy1 - dy0 );
        if ( suberr <= abs( frac_err * dy1 ) )
        {
//          printf("value accepted at level %d\n", level);
          y += dy1;
          abs_err += suberr;
          S += dS;
          level = 2;
          dS *= 2;
          if ( dS > dS_max ) dS = dS_max;
          if ( S + dS > S1) dS = S1 - S;
          if ( dS > 0 ) improve_domain();
          if ( S < S1 ) dy0 = modmid( 1 );
        }
        else
        {
          if ( level < level_max )
          {
            level++;
            dy0 = dy1;
          }
          else
          {
            dS--;
            improve_domain(); 
            level = 2;
            dy0 = modmid( 1 );
          }
        }
      }
      else 
      {
        dS--;
        improve_domain();
        level = 2;
        dy0 = modmid( 1 );
      }
    }
    
  } while ( S < S1 );

  // clean memory
  delete[] dydx;
  delete[] dydx_set;
}
  
////////////////////////////////////////////////////////////////////////////////
// 2D integration routines
////////////////////////////////////////////////////////////////////////////////

template <class T> void miser( T *a_base, int integrand_ID, double regn[], 
  unsigned long npts, double *ave, double *var )
{
  double x_sobseq[3];
  double x, y;

  static unsigned long MNPT = MC_bifurcation_sample;
  static unsigned long MNBS = 4 * MC_bifurcation_sample;
  static double PFAC = 0.1;

  double regn_temp[2*2+1];
  unsigned long n, nptl, nptr, npre;
  int j, jb;
  double avel, varl;
  double fracl, fval;
  double rgl, rgm, rgr, sigl, siglb, sigr, sigrb;
  double sum, sumb, summ, summ2;

  double fmaxl[2+1];
  double fmaxr[2+1];
  double fminl[2+1];
  double fminr[2+1];
  double pt[2+1];
  double rmid[2+1];
  summ = summ2 =0.0;
  
  if (npts < MNBS ) // to few points to bisect. Straight MC follows
  {

    for (n=1; n<=npts; n++)
    {
      sobseq2D(x_sobseq, false);
      x = ( regn[3] - regn[1]) * x_sobseq[1] + regn[1];
      y = ( regn[4] - regn[2] ) * x_sobseq[2] + regn[2];

      fval = (*a_base).relay_integrand_2D( integrand_ID, x, y );
      summ += fval;
      summ2 += fval * fval;
    }
    *ave = summ / npts;
    *var = dmax ( smallest_number, (summ2 - summ * summ / npts) / (npts*npts) );
  }
  else // do preliminary uniform sampling
  {
    npre = lmax( (unsigned long) (npts * PFAC ), MNPT );
    for (j=1; j <= 2; j++)
    {  // intialize boundaries for each dimension
      rmid[j] = 0.5 * (regn[j] + regn[2+j] ); 
        // the formula here is: "middle" = 0.5 * ("end" - "begin") + "begin"
      fminl[j] = fminr[j] = largest_number; // so ANY following number will be
        // smaller than this number and dmin will take the following number
      fmaxl[j] = fmaxr[j] = -largest_number;
    }

    for (n=1; n <= npre; n++ ) // loop over the points in the sample
    {
      
      sobseq2D(x_sobseq, false); 
      pt[1] = x = (regn[3] - regn[1]) * x_sobseq[1] + regn[1]; // x component
      pt[2] = y = (regn[4] - regn[2]) * x_sobseq[2] + regn[2]; // y component

      fval = (*a_base).relay_integrand_2D( integrand_ID, x, y );
      
      // following 2 lines are new...
      summ += fval; summ2 += fval * fval;
      
      for( j=1; j <= 2; j++ )
      {
        if (pt[j] <= rmid[j]) // so sample point in left part
        {
          fminl[j] = dmin(fminl[j], fval );
          fmaxl[j] = dmax(fmaxl[j], fval );
        }
        else
        {
         fminr[j] = dmin( fminr[j], fval );
          fmaxr[j] = dmax( fmaxr[j], fval );
        }
      }
    }
      // now fminl[j] contains smallest fval and fmaxl largest fval
    
    // also new. Check if we do not accidentally hit accuracy
    *ave = summ / npre;
    *var = dmax ( smallest_number, (summ2 - summ * summ / npre) / (npre*npre) );
    double err = 5e-3;
    if ( *ave == 0.0 or abs(*var) <= abs(err * err * (*ave) * (*ave) ) )
      return;
    
    sumb = largest_number;  // choose which dimension to bisect
    jb=0;
    siglb=sigrb=1.0;
    for (j=1; j<=2; j++)
    {
      if (fmaxl[j] > fminl[j] && fmaxr[j] > fminr[j] )
      {
        sigl=dmax(smallest_number, pow(fmaxl[j] - fminl[j], 2.0/3.0 ));
        sigr=dmax(smallest_number, pow(fmaxr[j] - fminr[j], 2.0/3.0 ));
        sum = sigl + sigr;
        if (sum <= sumb )
        {
          sumb = sum;
          jb = j;
          siglb=sigl;
          sigrb=sigr;
        }
      }
    }
    if (!jb) 
    {
//      printf("just intersecting along x\n");
      if (x_sobseq[1] > x_sobseq[2] ) jb = 1;
        else jb = 2; // random choose direction if we don't know what to do
    }
//    else printf("NOT intersecting along x!\n");

    rgl=regn[jb];
    rgm=rmid[jb];
    rgr=regn[2+jb];
    fracl = abs( (rgm-rgl) / (rgr-rgl));
    nptl=(unsigned long) (MNPT+(npts-3*MNPT)*fracl*siglb
      /(fracl*siglb+(1.0-fracl)*sigrb) );
    nptr=npts-MNPT-nptl;
    for (j=1; j <= 2; j++)
    {
      regn_temp[j]=regn[j];
      regn_temp[2+j]=regn[2+j];
    }
    regn_temp[2+jb]=rmid[jb];
    
    miser( a_base, integrand_ID, regn_temp, nptl, &avel, &varl );
    
    regn_temp[jb] = rmid[jb];
    regn_temp[2+jb] = regn[2+jb];
    
    miser( a_base, integrand_ID, regn_temp, nptr, ave, var );
    
    *ave = fracl * avel + (1-fracl)*(*ave);
    *var = fracl * fracl * varl + (1-fracl)*(1-fracl)*(*var);
  }
  return;
}

////////////////////////////////////////////////////////////////////////////////

template <class T > double integrate_2D_MC( T *a_base, int integrand_ID, 
  double *err)
{
  double fval, summ, summ2;
  double x_sobseq[3];
  int npts=10;
  double ave, var, varnpts3;
  int n;
  
  summ = 0.0; summ2 = 0.0;
  sobseq2D( x_sobseq, true ); // initialize sobseq sequence

  // initial round
  for (n=1; n<=npts; n++)
  {
    sobseq2D(x_sobseq, false);
    
//    printf("(theta', phi') = (%e, %e) -> ", x_sobseq[1], x_sobseq[2] );

    fval = (*a_base).relay_integrand_2D( integrand_ID, x_sobseq[1],x_sobseq[2]);
    summ += fval;
    summ2 += fval * fval;
  }

  do
  {
    sobseq2D(x_sobseq, false);

    fval = (*a_base).relay_integrand_2D( integrand_ID, x_sobseq[1],x_sobseq[2]);
    summ += fval;
    summ2 += fval * fval;
    varnpts3 = dmax( smallest_number, summ2 * npts - summ * summ );
    npts++;
  } while ( abs(varnpts3) > abs( (*err) * (*err) * summ * npts*npts) ); 

  ave = summ / npts;
  var = dmax ( smallest_number, (summ2 - summ * summ / npts) / (npts*npts) );
 
  *err = sqrt(var);
  return ave;
}

////////////////////////////////////////////////////////////////////////////////

template <class L, class I0, class I1> class MC_aux
{
  public:
  MC_aux( L *a_base, int a_int_ID, I0 *a_y0_base, int a_y0_ID,
    I1 *a_y1_base, int a_y1_ID, double a_x0, double a_x1) 
  { 
    // establish pointers to old integrand and y boundaries
    base = a_base; integrand_ID = a_int_ID;
    y0_base = a_y0_base; y0_ID = a_y0_ID;
    y1_base = a_y1_base; y1_ID = a_y1_ID;
    
    x0 = a_x0;
    x1 = a_x1;
    Dx = x1 - x0;
  }

  ~MC_aux()
  {
  }

  double relay_integrand_2D( int dummy, double xp, double yp )
  {
    x = xp * Dx + x0;
    y0 = (*y0_base).relay_integrand_1D( y0_ID, x );
    y1 = (*y1_base).relay_integrand_1D( y1_ID, x );
    if (y1 == y0) return 0.0;
    Dy = y1 - y0;
    y = yp*Dy + y0;
//    printf("(theta, phi) = (%e,%e)\n", x, y );
//    printf("[at theta=%e, phi between %e and %e]\n", x, y0, y1);
//    printf("theta between %e - %e\n", x0, x1);
//    abort();
    return Dy * (*base).relay_integrand_2D( integrand_ID, x, y );
  }

  protected:
  double x1, x0, Dx;
  double y1, y0, Dy;
  double x, y;

  L *base;
  int integrand_ID;
  I0 *y0_base;
  int y0_ID;
  I1 *y1_base;
  int y1_ID;

};

////////////////////////////////////////////////////////////////////////////////

template <class L, class I0, class I1> double integrate_2D_MC( L *a_base, 
  int integrand_ID, double x_0, double x_1, I0 *y0_base, int y0_ID, I1 *y1_base,
  int y1_ID, double *err )
{
  if (x_0 == x_1) { *err = 0.0; return 0.0; }
  
  MC_aux <L, I0, I1> aux(a_base, integrand_ID, y0_base, y0_ID, y1_base, y1_ID, 
    x_0, x_1);
  
  return (x_1 - x_0) * integrate_2D_MC( &aux, 0, err);  
}

////////////////////////////////////////////////////////////////////////////////


template <class L, class I0, class I1> double integrate_2D_miser( L *a_base, 
  int integrand_ID, double x_0, double x_1, I0 *y0_base, int y0_ID, I1 *y1_base,
  int y1_ID, double *err, unsigned long npts_start )
{
  if (x_0 == x_1) { *err = 0.0; return 0.0; }

  MC_aux <L, I0, I1> aux(a_base, integrand_ID, y0_base, y0_ID, y1_base, y1_ID, 
    x_0, x_1);
  
//  printf("(x_0 = %e, x_1=%e)\n", x_0, x_1);  
  return (x_1 - x_0) * integrate_2D_miser( &aux, 0, err, npts_start);  
}
////////////////////////////////////////////////////////////////////////////////

template <class T> double integrate_2D_miser( T *a_base, int integrand_ID, 
  double *err, unsigned long npts_start )
{
  // declare variables
  double x[3]; // when we have rewritten sobseq2D to properly start at 0,
    // we need of course to change all this as well!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  unsigned long npts;
  double ave = 0.0;
  double avetemp = 0.0;
  double var = 0.0;
  double vartemp = 0.0;
  unsigned long total_pts;
  double regn[4+1] = { 0.0, 0.0, 0.0, 1.0, 1.0 };
  double z; // result

  // initialize sobseq sequence
  sobseq2D(x, true);

  // first guess
  npts = lmax( (unsigned long) MC_sample_minimum, npts_start );
  miser( a_base, integrand_ID, regn, npts, &avetemp, &vartemp );
  total_pts = npts;
  var = vartemp * npts * npts;
  ave = avetemp * npts;

  // later guesses

//  printf("MC_sample_minimum = %ld\n", MC_sample_minimum );
//  printf("first guess yields: var %e, ave %e, frac error = %e\n", var, ave,
//    sqrt(abs(var)) / ave);
  while ( abs( var ) >= (*err) * (*err) * ave * ave )
  {
    npts = (unsigned long) (MC_sample_increase * npts);
    miser( a_base, integrand_ID, regn, npts, &avetemp, &vartemp );
    var += vartemp * npts * npts; 
    ave += avetemp * npts;
    total_pts += npts;
    if (ave == 0.0) 
    {  
       printf("could not find nonzero part integration domain\n");// abort();
       *err = -1.0; // negative error tells us something went wrong
       return 0.0;
    }
//    printf("frac error %e\n", sqrt(abs(var))/ave);
  }

  *err = sqrt( var ) / total_pts;
  z = ave / total_pts;
  
//  printf("error: %e, result: %e\n", *err, z );

  return z; 
}

////////////////////////////////////////////////////////////////////////////////
// 2D integration via looped 1D integration
////////////////////////////////////////////////////////////////////////////////

template <class T> double modmid(T *a_base, int integrand_ID, double y_0, 
  double y_1, double ddzdxdy0, double ddzdxdy1, int n_steps, double x )
// dydx is included as argument because we expect modmid routines to be
// performed in adjacent domains. We don't want to do double work when the
// new domain's left boundary is equal to the previous domain's right boundary
// and an earlier modmid routine already gave us the integrand at the boundary.
// note the extra parameter z, as we integrate over a surface (say, xz) where  
// integrands over x can depend on the value of z
{
  int n;
  double dzdxm, dzdxn, y, swap, h2, h, dzdxout;
  
  
  // determine stepsize
  h = (y_1 - y_0) / n_steps;
  
  // first step
  dzdxm = 0.0; dzdxn = h * ddzdxdy0;
  
  y = y_0 + h;
  dzdxout = (*a_base).relay_integrand_2D( integrand_ID, x, y );
    // yout functions as temporary storage of derivative
  
  // general step
  h2 = 2.0 * h;
  for ( n=2; n < n_steps; n++ )
  {
    swap = dzdxm + h2 * dzdxout;
    dzdxm = dzdxn;
    dzdxn = swap;
    y += h;
    dzdxout = (*a_base).relay_integrand_2D( integrand_ID, x, y );
  }
  
  // now take advantage of known dydx1
  swap = dzdxm + h2 * dzdxout;
  dzdxm = dzdxn;
  dzdxn = swap;
  y += h;
  dzdxout = ddzdxdy1;
  
  // last step and return
  return 0.5 * (dzdxm + dzdxn + h * dzdxout );
}

////////////////////////////////////////////////////////////////////////////////

template <class T> double integrate_1D_step( T *a_base, int integrand_ID, double
  y, double *stepsize, double steperr, double *abserr, double ddzdxdy0, 
  double ddzdxdy1, double h_min, double x )
{
  int m; // no of steps is 2*(m+1)
  double ReP; // result of Richardson extrapolation and error
  double errin = abs(steperr); // errin stores requested maximum error
  double dzdx[10]; // contains result for modmid with n = 2, 4, 6, ...
  double h[10];

  // get started with n=2. We need at least this much accuracy
  dzdx[0] = modmid( a_base, integrand_ID, y, y+(*stepsize), ddzdxdy0, ddzdxdy1, 
    2, x );
  h[0] = 0.5; // for ReP we go to limit h -> 0
  m = 1;
//  printf("fractional error limit: %e\n", steperr );
  
  // try at this stepsize
  while( m < 6 )
  {
    dzdx[m] = modmid( a_base, integrand_ID, y, y+(*stepsize), ddzdxdy0, 
      ddzdxdy1, 2*(m+1), x );
    h[m] = 1.0 / ( 2*(m+1) );
    ReP = pol_int( h, dzdx, m+1, abserr );
    *abserr = abs(*abserr);
    if( ReP == 0.0 )
    {
      if (dzdx[m] == 0.0) { return ReP; }
      else steperr = abs((*abserr) / dzdx[m] );
    } else
    steperr = abs((*abserr) / ReP);
    if ( steperr < errin ) return ReP;
    if ( abs( dzdx[m] - dzdx[m-1] ) < errin * dzdx[m] ) 
    {
      *abserr = abs( dzdx[m] - dzdx[m-1] );
      return dzdx[m];
    }
    m++;
  }

  // check if we still can decrease stepsize
  if (*stepsize <= h_min) return dzdx[m-1];

  // we can and have to decrease stepsize
  *stepsize = (*stepsize) / 2.0;
  ddzdxdy1 = (*a_base).relay_integrand_2D( integrand_ID, x, y+*stepsize );
  return integrate_1D_step( a_base, integrand_ID, y, stepsize, errin, abserr, 
    ddzdxdy0, ddzdxdy1, h_min, x );
}

////////////////////////////////////////////////////////////////////////////////

template <class T> double integrate_1D( T *a_base, int integrand_ID, double y_0, 
  double y_1, double res, double *err, double h_min, double x )
{
//  printf("domain: %2.12e - %2.12e, res = %e\n", y_0, y_1, res );
  double dzdx;
  double y;
  double H, H_max; // stepsize, maximum stepsize
  double steperr, sumabserr; // error, weighed sum of squared errors
  double abserr;
  double domain, inv_domain; // integration domain
  double ddzdxdy0, ddzdxdy1; // integrand values at the bounds of subdomain H
  
  if (y_0 == y_1)
  {
    *err = 0.0;
    return 0.0;
  }
  
  // make sure domain boundaries are in correct order
  if (y_1 < y_0) 
    return -integrate_1D( a_base, integrand_ID, y_1, y_0, res, err, h_min, x);
  
  // if res and / or h_min given as zero, it means make your own guess
  if ( res == 0.0 ) res = *err * ( y_1 - y_0 );
  if ( h_min == 0.0 ) h_min = *err * 1e-4;

  // if h_min of same order as domain, then this integral will not contribute
  // much of significance. We give just a 3 points modmid result
  if ( h_min >= (y_1 - y_0) )
  {
    ddzdxdy0 = (*a_base).relay_integrand_2D( integrand_ID, x, y_0 );
    ddzdxdy1 = (*a_base).relay_integrand_2D( integrand_ID, x, y_1 );
    double dzdx2, dzdx3; // two calls to modmid to have an error estimate as well
    dzdx2 = modmid( a_base, integrand_ID, y_0, y_1, ddzdxdy0, ddzdxdy1, 2, x );
    dzdx3 = modmid( a_base, integrand_ID, y_0, y_1, ddzdxdy0, ddzdxdy1, 3, x );
    *err = abs(dzdx2 - dzdx3);
    return dzdx3;
  }
  
  // declare variables
  H_max = H = res; // set stepsizes
  if (H < 0 ) { printf("negative resolution requested...\n"); abort(); }
  if (H > 0.1*(y_1 - y_0) ) // if overly large stepsize given, make educated guess
    { H = 1e-2*(y_1 - y_0); h_min = 1e-3 * H; }
  if ( *err > 0.1*H ) { *err = 1e-2 * H; }
  y = y_0;
  sumabserr = 0; // cumulative error
  domain = y_1 - y_0;
  inv_domain = 1.0 / domain;
  dzdx = 0.0;
  ddzdxdy0 = (*a_base).relay_integrand_2D( integrand_ID, x, y );
  
  while (y < y_1 )
  {
    if (y + H > y_1) H = y_1 - y; // correct for overshooting

    steperr = (*err); // reset frac error (it is modified to obtained error)
    ddzdxdy1 = (*a_base).relay_integrand_2D( integrand_ID, x, y+H );
    dzdx += integrate_1D_step( a_base, integrand_ID, y, &H, steperr, &abserr, 
      ddzdxdy0, ddzdxdy1, h_min, x );
    sumabserr += abserr;
    y += H;
    ddzdxdy0 = ddzdxdy1; // use right boundary value of the integrand of this step as 
      // left boundary value of the integrand for the next step
    H *= 2.0; // last attempt was succesful, so try a bigger step
    if (H > H_max) H = H_max; // but not too big
  }
  
  *err = sumabserr;
  
//  printf("%2.12e, %2.12e\n", x, dzdx );
  return dzdx;
}

////////////////////////////////////////////////////////////////////////////////

template <class L, class I0, class I1> class c_aux
{
  public:
  typedef double (c_aux::*pt2integrand1D) ( double );

  void init( L *a_a_base, int a_integrand_ID, I0 *a_y0_base, int a_y0_ID,
    I1 *a_y1_base, int a_y1_ID, double a_y_max, double a_y_min, double *a_err )
  {
    a_base = a_a_base; integrand_ID = a_integrand_ID; 
    y0_base = a_y0_base; y0_ID = a_y0_ID;
    y1_base = a_y1_base; y1_ID = a_y1_ID;
    integrand_list[0] = &c_aux::integrand;
    err = a_err;
    y_max = a_y_max; y_min = a_y_min;
    if (a_y_max == 0.0 or a_y_min == 0.0 ) rescale_stepsizes = true;
      else rescale_stepsizes = false;
  }
  
  c_aux() 
  { 
    integrand_list = new pt2integrand1D[1]; 
  }

  ~c_aux()
  {
    delete[] integrand_list;
  }

  double relay_integrand_1D( int integrand_ID, double x )
  {
    return (*this.*integrand_list[ integrand_ID ])( x );
  }
          
  protected:
  L *a_base;
  int integrand_ID;
  I0 *y0_base;
  int y0_ID;
  I1 *y1_base;
  int y1_ID;
  double y_max;
  double *err;
  double y_min;
  pt2integrand1D *integrand_list;
  bool rescale_stepsizes; // true if y_max / y_min given as 0.0 in function call
  double integrand( double x )
  {
    double y_0, y_1; double suberr = 0.8*(*err); double y;
    y_0 = (*y0_base).relay_integrand_1D( y0_ID, x );
    y_1 = (*y1_base).relay_integrand_1D( y1_ID, x );
 //   printf("Now looking at theta = %e, phi = %e...%e\n", x, y_0, y_1 );
//    printf("with res = %e, suberr = %e, h_min = %e\n", res, suberr, h_min); fflush(stdout);

    // if y_max / y_min not given make own estimate for each value of x
    if ( rescale_stepsizes )
    {
      y_max = abs( suberr * ( y_1 - y_0 ) );
      y_min = suberr * suberr * y_max;
    }    

    // do the integral
    y = integrate_1D( a_base, integrand_ID, y_0, y_1, y_max, &suberr, y_min, x);
    
//    printf("at theta = %2.12e, phi = %e .. %e we get %2.12e (+- %e)\n", x, y_0, y_1, y, suberr );
    return y;
  }
};

////////////////////////////////////////////////////////////////////////////////

template <class T, class I0, class I1> double integrate_2D_RE( T *a_base, int
  integrand_ID, double x_0, double x_1, double x_max, double x_min, I0 *y0_base,
  int y0_ID, I1 *y1_base, int y1_ID, double y_max, double y_min, double *err )
{
  c_aux <T, I0, I1> aux;
  double z; // will contain the result
  
  // if x_max and x_min for x direction equal to zero, a guess is made
  if (x_max == 0.0 or x_min == 0.0)
  {
    x_max = *err * abs(x_1 - x_0);
    x_min = (*err) * (*err) * x_max;
  }
  
  aux.init( a_base, integrand_ID, y0_base, y0_ID, y1_base, y1_ID, y_max, y_min, 
    err );
    
  z = integrate_1D( &aux, 0, x_0, x_1, x_max, err, x_min );
  return z;
}

