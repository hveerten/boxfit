////////////////////////////////////////////////////////////////////////////////
//
// coords.h
//
// Created Jun 23, 2011 by HJvE
// Last modified July 29, 2016 by HJvE
//
// provides coord structs and functions
//
////////////////////////////////////////////////////////////////////////////////

#include "coords.h"

////////////////////////////////////////////////////////////////////////////////

void hfromcartesian(s_coordinates &cor)
{
  cor.h = sqrt(cor.x * cor.x + cor.y * cor.y);
}

////////////////////////////////////////////////////////////////////////////////

void sphericalfromcartesian(s_coordinates &cor)
// uses the Cartesian elements of the struct to calculate the spherical ones,
// taking care to avoid singularities in the transformation
{
  double xysqrd = cor.x * cor.x + cor.y * cor.y;
  cor.r = sqrt(xysqrd + cor.z * cor.z);
  
  // point approximately at the origin or along the z-axis
  if (cor.r < 1e-8 or xysqrd < 1e-8) 
  { 
    cor.phi = 0.0; 
    cor.theta = 0.0; 
    cor.sin_theta = sin(cor.theta); 
    cor.cos_theta = cos(cor.theta);
    cor.cos_phi = cos(cor.phi);
  }
  else
  {
    // r always larger than z and origin already taken care of, so theta is ok    
    cor.theta = acos(cor.z / cor.r); // acos returns value between 0, pi
    cor.sin_theta = sin(cor.theta); 
    cor.cos_theta = cos(cor.theta);
    
    // deal with first quadrant, including y = 0
    if (cor.x > 0 and cor.y >= 0)
    {
      if (cor.y < 1e-8) cor.phi = 0.0;
      else cor.phi = atan(cor.y / cor.x);
    }
    else // 2nd quadrant
    if (cor.x <= 0 and cor.y > 0)
    {
      if (fabs(cor.x) < 1e-8) cor.phi = 0.5 * PI;
      else cor.phi = PI - atan(-cor.y / cor.x);
    }
    else
    if (cor.x < 0 and cor.y <= 0)
    {
      if (fabs(cor.y) < 1e-8) cor.phi = PI;
      else cor.phi = PI + atan(cor.y / cor.x);
    }
    if (cor.x >= 0 and cor.y < 0)
    {
      if (cor.x < 1e-8) cor.phi = 1.5 * PI;
      else cor.phi = 2.0 * PI - atan(-cor.y / cor.x);
    }

    cor.cos_phi = cos(cor.phi);
  }
}

////////////////////////////////////////////////////////////////////////////////

void cylindricalfromspherical(s_coordinates &cor)
{
  cor.z = cor.r * cos(cor.theta);
  cor.h = cor.r * sin(cor.theta);
}

