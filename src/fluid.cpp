////////////////////////////////////////////////////////////////////////////////
//
// fluid.cpp
// 
// Created July 21, 2011 by HJvE
// Last modified December 29, 2011 by HJvE
//
// empty destructor, print string name from number
//
////////////////////////////////////////////////////////////////////////////////

#include "fluid.h"

void sprintvar(char* varstring, int varno)
{
  switch (varno)
  {
    case rho_: sprintf(varstring, "rho"); break;
    case eint_: sprintf(varstring, "eint"); break;
    case v_x_: sprintf(varstring, "v_x"); break;
    case v_y_: sprintf(varstring, "v_y"); break;
    case v_: sprintf(varstring, "v"); break;
    case lfac_: sprintf(varstring, "lfac"); break;
    case D_: sprintf(varstring, "D"); break;
    case p_: sprintf(varstring, "p"); break;
    case tau_: sprintf(varstring, "tau"); break;
    case N_: sprintf(varstring, "N"); break;
    case lfacbeta_: sprintf(varstring, "lfac x beta"); break;
    case dpdr_: sprintf(varstring, "d p / d r"); break;
    case drhodr_: sprintf(varstring, "d rho / d r"); break;
    case c_s_: sprintf(varstring, "c_s"); break;
    default:
      printf("sprintfvar ERROR: undefined variable %d\n", varno);
      fflush(stdout);
      abort();
  }
}

////////////////////////////////////////////////////////////////////////////////

int sscanvar(char* varstring)
{
  if (strcmp(varstring, "rho") == 0) return rho_;
  if (strcmp(varstring, "eint") == 0) return eint_;
  if (strcmp(varstring, "v_x") == 0) return v_x_;
  if (strcmp(varstring, "v_y") == 0) return v_y_;
  if (strcmp(varstring, "v") == 0) return v_;
  if (strcmp(varstring, "lfac") == 0) return lfac_;
  if (strcmp(varstring, "D") == 0) return D_;
  if (strcmp(varstring, "p") == 0) return p_;
  if (strcmp(varstring, "tau") == 0) return tau_;
  if (strcmp(varstring, "N") == 0) return N_;
  
  return -1;
}

////////////////////////////////////////////////////////////////////////////////

c_fluid :: ~c_fluid()
{}
