////////////////////////////////////////////////////////////////////////////////
//
// radiation.cpp
//
// Created, Jun 27 2011, by HJvE
// Last Modified, NOv 24, 2014 by HJvE
//
// Routines to provide emission and absorption coefficients for a given
// fluid state.
//
// References:
//   van Eerten & Wijers 2009, MNRAS 394, 2164
//   van Eerten et al. 2010, MNRAS, 403, 300
//   van Eerten, Zhang & MacFadyen 2010, ApJ 722, 235
//
////////////////////////////////////////////////////////////////////////////////

#include "radiation.h"

////////////////////////////////////////////////////////////////////////////////

c_emab :: c_emab()
{
  overrule_theta_max = -1.0;
  //jet = BOTH_;
}

////////////////////////////////////////////////////////////////////////////////

c_emab :: ~c_emab()
{
}

////////////////////////////////////////////////////////////////////////////////

void c_emab :: get_emab(s_coordinates cor, double &em, double &ab)
{
  
  // local fluid quantities, radiation related
  double B, invlfac_m, invlfac_M;

  double mu, beaming; // cosine of beaming angle and beaming term
  // frequencies
  double nu_c, invnu_m = -1.0, invnu_M = -1.0, nu_cnu_m = -1.0;
  // resulting output
  double em_s; // emission coefficient for synchrotron radiation
  double em_s_local, ab_local = 0.0, ab_freq = 0.0;
  double em_s_freq = 0.0;
    // the emission and absorption coefficient are the product of a factor
    // dependending directly on the local fluid conditions and a dimensionless
    // function of frequency over critical frequency. If cooling is enabled
    // this functions has two dimensionless parameters, using two critical
    // frequencies.
  
  // some constants
  static double C_em = sqrt(3.0) / (4.0 * PI) * pow(e_e, 3.0) 
    / (m_e * v_light * v_light) * 0.5;
  static double C_nu = 3.0 * e_e / ( 4.0 * PI * m_e * v_light );
  static double C_ab_0 = sqrt(3.0) * pow(e_e, 3.0)
    / (16.0 * PI * m_e * m_e * v_light * v_light);
    // see eq. (A21) of van Eerten, Zhang, MacFadyen 2010

  //----------------------------------------------------------------------------

  sphericalfromcartesian(cor); // we're always going to need spherical
    // coordinates in order to calculate the beaming

  //----------------------------------------------------------------------------

  // set up local fluid conditions.
  p_fluid->set_local(cor);

  // set up local radiation related fluid conditions
  B = sqrt(epsilon_B * p_fluid->local[eint_] * 8.0 * PI);
  
  invlfac_m = ((p_synch - 1.0) * ksi_N * p_fluid->local[N_] * 
    m_e * v_light * v_light) / ((p_synch - 2.0) * epsilon_E * 
    p_fluid->local[eint_]);

  if (electron_cooling) 
    invlfac_M = sigma_T * B * B * cor.t / 
      (6.0 * PI * m_e * p_fluid->local[lfac_] * v_light);
  else
    invlfac_M = 0;

  // set beaming angle. Note that this angle is between the observer direction
  // and the flow velocity vector (rather than the fluid parcel position).
  // Only for spherically symmetric flow (i.e. 1D snapshots) the two are 
  // identical.
  if (p_fluid->local[v_] > 0)
  {
    mu = 
      (p_Obs->sintheta * cor.cos_phi * (p_fluid->local[v_x_] * cor.sin_theta
      + p_fluid->local[v_y_] * cor.cos_theta) + p_Obs->costheta * 
      (p_fluid->local[v_x_] * cor.cos_theta - p_fluid->local[v_y_] * 
      cor.sin_theta)) / p_fluid->local[v_];
  }
  else
  {
    mu = 1.0;
  }

  //printf("lfac = %e\n", p_fluid->local[lfac_]);

  if ( p_fluid->local[lfac_] > 1.00000001 ) 
    beaming = p_fluid->local[lfac_] - sqrt(p_fluid->local[lfac_] * 
      p_fluid->local[lfac_] - 1.0) * mu;
  else 
  {
    em = 0.0; ab = 0.0; return;
  }

  //----------------------------------------------------------------------------
  // calculate emissivity

  em_s_local = (p_synch - 1.0) * C_em * p_fluid->local[N_] * ksi_N * B;
  nu_c = p_Obs->nu * beaming * (1.0 + p_Obs->z); 
    // Doppler shift to comoving freq
  invnu_m = invlfac_m * invlfac_m / (B * C_nu);
  nu_cnu_m = nu_c * invnu_m;

  if (!electron_cooling)
  {
    em_s_freq = 9.6323 * powerlawsync(nu_c, invnu_m) / (3.0 * p_synch - 1.0);
  }
  else
  {
    invnu_M = invlfac_M * invlfac_M / (B * C_nu);
    em_s_freq = 9.6323 * powerlawsync(nu_c, invnu_m, invnu_M) / 
      (3.0 * p_synch - 1.0);
  }

  // numerical constant chosen to match Zhang & MacFadyen 2009, Granot 1999

  em_s = em_s_local * em_s_freq;
  em = em_s; // sum over all radiative processes (for now synchrotron only)
  em = em / (beaming * beaming); // transform em to lab frame

  //----------------------------------------------------------------------------
  // calculate absorption coefficient

  if (absorption)
  {
    // simplified power laws, global treatment of cooling
    // The effect of cooling on the absorption coefficient is ignored.
    // This is valid when the cooling break is sufficiently far from the
    // self-absorption break in the spectrum. If it isn't, the bigger
    // issue is the global treatment of cooling anyway.
    ab_local = C_ab_0 * ksi_N * p_fluid->local[N_] * B;
    ab_freq = (p_synch - 1.0) * (p_synch + 2.0) * invlfac_m / nu_c / nu_c; 
    
    if (nu_cnu_m > 1.0)
      ab_freq = ab_freq * pow(nu_cnu_m, -p_synch * 0.5);
    else
      ab_freq = ab_freq * pow(nu_cnu_m, onethird);

    ab = ab_local * ab_freq; // combine various aspects
    ab = ab * beaming; // transform absorption coefficient to the lab frame
  }
  else ab = 0.0;
  
  // in case of error
  if (em != em or ab != ab)
  //if (invnu_m > 1.)
  {
    printf("---------------------------------------------------------------\n");
    printf("fluid state:\n");
    printf("  N = %e, eint = %e, B = %e\n", p_fluid->local[N_], 
      p_fluid->local[eint_], B);
    printf("  lfac = %e, fluid.v = %e, fluid.v_x = %e, fluid.v_y = %e\n",
      p_fluid->local[lfac_], p_fluid->local[v_], p_fluid->local[v_x_], 
      p_fluid->local[v_y_]);

    printf("coordinates:\n");
    printf("  x = %e, y = %e, z = %e\n", cor.x, cor.y, cor.z);
    printf("  r = %e, theta = %e, phi = %e\n", cor.r, cor.theta, cor.phi);
    printf("  cos(phi) = %e\n", cor.cos_phi);
    printf("  cos(theta) = %e, sin(theta) = %e\n", cor.cos_theta, 
      cor.sin_theta);
    printf("  h = %e\n", cor.h);
    printf("radiation parameters\n");
    printf("  epsilon_E = %e, epsilon_B = %e\n", epsilon_E, epsilon_B);
    printf("  ksi_N = %e, p_synch = %e\n", ksi_N, p_synch);
    printf("emission parameters\n");
    printf("  em = %e, em_s_local = %e, em_s_freq = %e\n", em, em_s_local, 
      em_s_freq);
    printf("  ab = %e, ab_local = %e, ab_freq = %e\n", ab, ab_local, ab_freq);
    printf("  beaming = %e, mu = %e\n", beaming, mu);
    printf("  nu_c = %e, nu_cnu_m = %e\n", nu_c, nu_cnu_m);
    printf("  fluid.invlfac_m = %e, invnu_m = %e, ksi_N = %e\n", invlfac_m, 
      invnu_m, ksi_N);
    printf("observer position:\n");
    printf("  sin_theta_obs = %e\n", p_Obs->sintheta);
    printf("  cos_theta_obs = %e\n", p_Obs->costheta);
    printf("  theta_obs = %e\n", p_Obs->theta);
    fflush(stdout);  
    exit(1);
  }
  
}

////////////////////////////////////////////////////////////////////////////////

double c_emab :: powerlawsync(double nu_c, double invnu_m)
{
  double nu_cnu_m = nu_c * invnu_m;

  if (nu_cnu_m > 1.0) return pow(nu_cnu_m, (1.0 - p_synch) * 0.5);
  else return pow(nu_cnu_m, onethird);
}

////////////////////////////////////////////////////////////////////////////////

double c_emab :: powerlawsync(double nu_c, double invnu_m, double invnu_M)
{
  double nu_cnu_m = nu_c * invnu_m;
  double nu_cnu_M = nu_c * invnu_M;
  
  if (invnu_m < invnu_M) // fast cooling regime
  {
    if (nu_cnu_M < 1.0) return pow(nu_cnu_M, onethird);
    else
    {
      if (nu_cnu_m < 1.0) return pow(nu_cnu_M, -0.5);
      else return pow(invnu_M / invnu_m, -0.5) * 
        pow(nu_cnu_m, -0.5 * p_synch);
    }
  }
  else // slow cooling regime
  {
    if (nu_cnu_m < 1.0) return pow(nu_cnu_m, onethird);
    else
    {
      if (nu_cnu_M < 1.0) return pow(nu_cnu_m, 0.5 * (1.0 - p_synch));
      else return pow(invnu_m / invnu_M, 0.5 * (1.0 - p_synch)) *
        pow(nu_cnu_M, -0.5 * p_synch);
    }
  }
}
