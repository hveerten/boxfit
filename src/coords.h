////////////////////////////////////////////////////////////////////////////////
//
// coords.h
//
// Created Jun 23, 2011 by HJvE
// Last modified April 24, 2015 by HJvE
//
// provides coord structs and functions
//
////////////////////////////////////////////////////////////////////////////////

#ifndef COORDS_H_
#define COORDS_H_

#include "extramath.h"

struct s_coordinates
{
  // these coordinates and time will always be in the LAB frame
  
  double x, y, z; // z also cylindrical coordinate
  double r, theta, phi; // phi also cylindrical coordinate
  double h; // cylindrial coordinate perpendicular to jet axis
  double sin_theta, cos_theta;
  double cos_phi;
  
  double t; // local lab frame time
};

////////////////////////////////////////////////////////////////////////////////

void sphericalfromcartesian(s_coordinates &cor);
// uses the Cartesian elements of the struct to calculate the spherical ones,
// taking care to avoid singularities in the transformation

void hfromcartesian(s_coordinates &cor);
// only calculate the cylindrical h coordinate

void cylindricalfromspherical(s_coordinates &cor);

#endif // COORDS_H_
