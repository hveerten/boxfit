////////////////////////////////////////////////////////////////////////////////
//
// observer.cpp
//
// Created July 2, 2014 by HJvE
// Last modified July 2, 2014 by HJvE
//
// provides observer class, which contains all information on the observer
// position
//
////////////////////////////////////////////////////////////////////////////////

#include "observer.h"

void c_observer :: set_sincos()
{
  sintheta = sin(theta);
  costheta = cos(theta);
}
