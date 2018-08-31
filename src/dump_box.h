////////////////////////////////////////////////////////////////////////////////
//
// dump_box.h
//
// Created Feb 22, 2016 by HJvE
// Last modified Feb 22, 2016 by HJvE
//
// dump_box creates an hdf5 dump from a grid snapshot file  ////////////////////////////////// UPDATE
// SYNTAX: dump_box <filename> <index>
//
////////////////////////////////////////////////////////////////////////////////

#include "environment.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "extramath.h"
#include "box.h"
#include "parse.h"
#include "physics.h"

////////////////////////////////////////////////////////////////////////////////

class c_dump_box
{
  public:

    c_dump_box();
    ~c_dump_box();

    int main(int argc, char* argv[]);

  //----------------------------------------------------------------------------
  
  protected:
  
    double dumpxmin, dumpxmax, dumpymin, dumpymax; // snapshot boundaries
    int dumpxres, dumpyres; // snapshot resultion
    int cx, cy; // dummy loop variables
    double r, theta; // spherical coordinates
    int fluidvar; // which fluid quantity to dump
    double **w; // 2D arrays containing fluid data
    double **x, **y;// 2D arrays containing coordinates for all array entries
    double wmax, wmin; // fluid quantity extrema on requested grid

    double argdouble; // dummy var used when parsing arguments
    int argint; // dummy var used when parsing arguments
    bool quiet; // if true don't print output on screen

    char *outputfilename;
    
    char* varname; // variable name, if provided by name as argument

    s_coordinates cor; // coordinates on box

    // slice related
    bool slice; // true is 1D slice requested rather than a 2D dump
    double dumprmin, dumprmax; // radial slice boundaries
    int dumprres; // slice resolution
    double *rs, *ws; // 1D arrays, coordinate and fluid data
    int cr; // dummy loop variable

    c_box box;
    
    //--------------------------------------------------------------------------
    
    int read_parameters(int argc, char* argv[]);
    int dump_slice();
    int dump_2D();
};
